"""Numerical views and reproducible CPU graph builders."""
import csv
import random
import shutil
import subprocess
from pathlib import Path

import numpy as np
from scipy import io, sparse
from scipy.sparse.linalg import LinearOperator, svds
from scipy.spatial import cKDTree
from sklearn.kernel_approximation import PolynomialCountSketch
from sklearn.preprocessing import normalize
from sklearn.random_projection import SparseRandomProjection

from .data import write_tsv
from .runtime import rscript_command


def projection(x, name, cfg, dest):
    inputs = normalize(x, norm="l2", axis=1)
    if name == "V_LINEAR":
        model = SparseRandomProjection(n_components=cfg["linear_dim"], density="auto",
                                       dense_output=True, random_state=cfg["seed"])
    else:
        model = PolynomialCountSketch(degree=2, gamma=1.0, coef0=0.0,
                                      n_components=cfg["tensor_dim"], random_state=cfg["seed"])
    model.fit(inputs)
    out = np.lib.format.open_memmap(dest, mode="w+", dtype=np.float32,
                                    shape=(x.shape[0], model.n_components))
    for start in range(0, x.shape[0], cfg["chunk_cells"]):
        stop = min(start + cfg["chunk_cells"], x.shape[0])
        z = model.transform(inputs[start:stop])
        out[start:stop] = normalize(z, norm="l2", axis=1)
    out.flush()
    return out


def gini_features(x, gene_ids):
    """Positive-only sorted values, with the preceding implicit zeros included in ranks."""
    c = x.tocsc()
    n, p = c.shape
    g = np.zeros(p)
    peak = np.zeros(p)
    det = np.diff(c.indptr) / n
    for j in range(p):
        v = np.sort(c.data[c.indptr[j]:c.indptr[j+1]])
        if len(v):
            ranks = np.arange(n-len(v)+1, n+1, dtype=float)
            g[j] = np.dot(2*ranks-n-1, v) / (n * v.sum())
            peak[j] = v[-1]
    u = np.log2(peak + 0.1)
    mean = np.asarray(x.mean(axis=0)).ravel()
    sqmean = np.asarray(x.power(2).mean(axis=0)).ravel()
    valid = (sqmean - mean**2) > 1e-12 * np.maximum(1, sqmean)

    def trend(mask):
        values = u[mask]
        if not len(values):
            return np.zeros_like(u)
        breaks = np.unique(np.quantile(values, np.linspace(0, 1, min(max(len(values)//8, 1), 64)+1)))
        idx = np.searchsorted(breaks[1:-1], values, side="right")
        knots = [(np.median(values[idx == i]), np.median(g[mask][idx == i])) for i in np.unique(idx)]
        xx, yy = zip(*knots)
        return np.interp(u, xx, yy)

    fit = trend(valid)
    if valid.any():
        residual = g - fit
        fit = trend(valid & (residual <= np.quantile(residual[valid], .75)))
    ordered = sorted(np.where(valid)[0], key=lambda j: (-(g[j]-fit[j]), gene_ids[j]))
    return ordered, g, det, fit


def gini_pca(ds, cfg, dest, feature_file):
    idx, g, det, fit = gini_features(ds.x, ds.gene_ids)
    sel = np.asarray(idx[:cfg["n_features"]], dtype=int)
    selected = set(sel)
    write_tsv(feature_file, ["gene_id", "gini", "detected", "fit", "residual", "selected"],
              (dict(gene_id=gene, gini=g[j], detected=det[j], fit=fit[j], residual=g[j]-fit[j],
                    selected=int(j in selected)) for j, gene in enumerate(ds.gene_ids)))
    active = np.flatnonzero(np.diff(ds.x.indptr))
    a = ds.x[active][:, sel].tocsr()
    if min(a.shape) < 2:
        raise ValueError("Gini PCA requires two variable features and two informative cells")
    n = a.shape[0]
    means = np.asarray(a.mean(axis=0)).ravel()
    scale = np.sqrt(np.maximum(0, (np.asarray(a.power(2).mean(axis=0)).ravel()-means**2)*n/(n-1)))
    scale[scale < 1e-8] = 1.
    scaled = (a @ sparse.diags(1. / scale)).tocsr()
    shift = means / scale

    def forward(v):
        return (scaled @ v) - (shift @ v)

    def backward(v):
        if v.ndim == 1:
            return scaled.T @ v - shift * v.sum()
        return scaled.T @ v - shift[:, None] * v.sum(axis=0)

    operator = LinearOperator(a.shape, matvec=forward, matmat=forward, rmatvec=backward,
                              rmatmat=backward, dtype=np.float64)
    d = min(cfg["n_pcs"], min(a.shape)-1)
    rng = np.random.RandomState(cfg["seed"])
    u, s, vt = svds(operator, k=d, v0=rng.normal(size=min(a.shape)))
    u, s, vt = u[:, ::-1], s[::-1], vt[::-1]
    # Deterministic PC sign convention.
    signs = np.sign(vt[np.arange(d), np.argmax(abs(vt), axis=1)])
    out = np.zeros((ds.x.shape[0], d), dtype=np.float32)
    out[active] = u * (s * signs)
    np.save(dest, out)
    np.savez(dest.parent / "gini_pca_transform.npz", features=sel, mean=means, scale=scale,
             loadings=(vt * signs[:, None]).T)
    return out


def neighbours(z, informative, cfg):
    norms = np.linalg.norm(z, axis=1)
    valid = np.flatnonzero(informative & (norms > 1e-10))
    unit = np.asarray(z[valid], dtype=np.float64)
    if len(unit):
        unit /= np.linalg.norm(unit, axis=1)[:, None]
    k = min(cfg["k"], len(valid)-1)
    if k <= 0:
        return valid, unit, np.empty((len(valid), 0), dtype=int), np.zeros(len(valid))
    tree = cKDTree(unit)
    nbrs, radius = np.empty((len(valid), k), dtype=int), np.zeros(len(valid))
    # Query k+1 first. A radius search includes every tie at the boundary.
    for start in range(0, len(valid), cfg["chunk_cells"]):
        stop = min(start + cfg["chunk_cells"], len(valid))
        dist, _ = tree.query(unit[start:stop], k=k+1)
        for i, r in zip(range(start, stop), dist[:, -1]):
            pool = np.asarray(tree.query_ball_point(unit[i], np.nextafter(r + 1e-12, np.inf)))
            pool = pool[pool != i]
            distances = np.linalg.norm(unit[pool] - unit[i], axis=1)
            select = np.lexsort((valid[pool], distances))[:k]
            nbrs[i] = pool[select]
            radius[i] = distances[select[-1]]
    return valid, unit, nbrs, radius


def python_graph(nbrs, cfg, dest):
    import igraph as ig
    n = len(nbrs)
    sets = [set(map(int, row)) | {i} for i, row in enumerate(nbrs)]
    edges, weights = [], []
    for i in range(n):
        for j in nbrs[i]:
            j = int(j)
            if j > i or i not in sets[j]:
                left, right = sorted((i, j))
                score = len(sets[i] & sets[j]) / len(sets[i] | sets[j])
                if score > 1/15:
                    edges.append((left, right))
                    weights.append(score)
    if edges:
        rows, cols = zip(*edges)
        upper = sparse.csr_matrix((weights, (rows, cols)), shape=(n,n))
        snn = upper + upper.T
    else:
        snn = sparse.csr_matrix((n, n))
    sparse.save_npz(dest, snn)
    if not edges:
        return {r: list(range(n)) for r in cfg["resolutions"]}
    graph = ig.Graph(n=n, edges=edges, directed=False)
    partitions = {}
    for res in cfg["resolutions"]:
        ig.set_random_number_generator(random.Random(cfg["seed"]))
        try:
            result = graph.community_leiden(objective_function="modularity", weights=weights,
                                            resolution=float(res), n_iterations=10)
        finally:
            ig.set_random_number_generator(None)
        partitions[res] = result.membership
    return partitions


def run_r(mode, work, cfg, *args):
    command = rscript_command(cfg)
    r = shutil.which(command)
    if not r:
        raise RuntimeError(f"Rscript executable not found: {command}")
    script = Path(__file__).resolve().parent / "r" / "seurat.R"
    command = [r, "--vanilla", str(script), mode, str(work), *map(str, args)]
    # R process sees its own configured environment, independent of discovery-provider secrets.
    result = subprocess.run(command, text=True, capture_output=True, timeout=cfg["r_timeout"])
    if result.returncode != 0:
        raise RuntimeError(f"Seurat {mode}: " + (result.stderr or result.stdout)[-1600:].strip())


def seurat_pca(ds, cfg, work, dest, output):
    if ds.raw is None:
        raise ValueError("Seurat VST requires raw counts; input declares log1p_normalized")
    active = np.flatnonzero(np.diff(ds.raw.indptr))
    if min(len(active), ds.raw.shape[1]) < 3:
        raise ValueError("Seurat PCA requires at least three informative cells and genes")
    io.mmwrite(str(work/"counts.mtx"), ds.raw[active].T.tocoo())
    run_r("pca", work, cfg, cfg["n_features"], cfg["n_pcs"], cfg["seed"], cfg["target_sum"],
          cfg["scale_max_bytes"])
    pca = np.loadtxt(work/"pca.tsv", ndmin=2)
    if pca.shape[0] != len(active) or not np.isfinite(pca).all():
        raise ValueError("Seurat returned invalid PCA rows")
    out = np.zeros((ds.x.shape[0], pca.shape[1]), dtype=np.float32)
    out[active] = pca
    np.save(dest, out)
    feature_idx = np.loadtxt(work/"hvg_indices.tsv", dtype=int, ndmin=1) - 1
    write_tsv(output/"hvg.tsv", ["gene_id", "symbol"],
              ({"gene_id": ds.genes[j]["gene_id"], "symbol": ds.genes[j].get("symbol","")} for j in feature_idx))
    shutil.copyfile(work/"pca_loadings.tsv", output/"pca_loadings.tsv")
    shutil.copyfile(work/"versions.tsv", output/"versions.tsv")
    return out


def seurat_graph(unit, cfg, work, dest):
    np.savetxt(work/"embedding.tsv", unit, delimiter="\t", fmt="%.16g")
    run_r("graph", work, cfg, cfg["k"], ",".join(map(str,cfg["resolutions"])), cfg["seed"])
    sparse.save_npz(dest, io.mmread(str(work/"snn.mtx")).tocsr())
    with (work/"groups.tsv").open() as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if len(rows) != len(unit) or any(len(row) != len(cfg["resolutions"]) for row in rows):
        raise ValueError("Seurat partition shape mismatch")
    return {r: [row[j] for row in rows] for j, r in enumerate(cfg["resolutions"])}
