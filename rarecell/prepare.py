"""Prepare immutable, annotation-free discovery snapshots with independently tracked stages."""
import csv
import fcntl
import importlib.metadata
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import sparse

from .data import canonical, digest, file_sha, read_dataset, write_tsv
from .views import (gini_pca, neighbours, projection, python_graph,
                    run_r, seurat_graph, seurat_pca)


NAMES = ("V_HVG", "V_GINI", "V_LINEAR", "V_TENSOR2")
DEFAULTS = dict(views=list(NAMES[:3]), graph_backend="seurat", seed=1009, resolutions=[0.4, 1.0],
                k=15, n_features=1000, n_pcs=30, linear_dim=128, tensor_dim=256, chunk_cells=256,
                target_sum=10000., markers=12, outlier_mad=3.5, local_cosine=.90, local_pairs=64,
                scale_max_bytes=536870912, rscript="Rscript", r_timeout=7200)
MARKER_FIELDS = ["resolution","group_id","gene_id","symbol","avg_log2FC","pct.1","pct.2","p_val","p_val_adj"]
FAILURES = (ValueError, RuntimeError, OSError, ImportError, subprocess.TimeoutExpired,
            ArithmeticError, np.linalg.LinAlgError)


def config(cfg):
    options = cfg.get("prepare", {})
    if set(options) - DEFAULTS.keys():
        raise ValueError(f"unknown prepare parameters: {sorted(set(options) - DEFAULTS.keys())}")
    p = {**DEFAULTS, **options}
    if not p["views"] or not set(p["views"]).issubset(NAMES) or len(set(p["views"])) != len(p["views"]):
        raise ValueError("select unique, recognized views")
    for key in ("seed","k","n_features","n_pcs","linear_dim","tensor_dim","chunk_cells","markers","scale_max_bytes","r_timeout"):
        if type(p[key]) != int or p[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    if p["seed"] >= 2**31 or p["graph_backend"] not in ("seurat","igraph"):
        raise ValueError("invalid seed or graph_backend")
    if not p["resolutions"] or len(set(p["resolutions"])) != len(p["resolutions"]) or not all(
        math.isfinite(r) and r>0 for r in p["resolutions"]):
        raise ValueError("resolutions must be unique, positive and finite")
    if not math.isfinite(p["outlier_mad"]) or p["outlier_mad"] < 0 or not 0 <= p["local_cosine"] <= 1:
        raise ValueError("invalid outlier/local pair threshold")
    if type(p["local_pairs"]) != int or p["local_pairs"] < 0:
        raise ValueError("local_pairs must be a nonnegative integer")
    p["resolutions"] = sorted(map(float, p["resolutions"]))
    return p


def blocked(exc):
    return {"status": "blocked", "reason": f"{type(exc).__name__}: {exc}"}


def marker_cards(ds, ids, background, total, all_detected, limit):
    if background is None and len(ids) == len(ds.cells):
        return []
    inside = ds.x[ids]
    si = np.asarray(inside.sum(axis=0)).ravel()
    di = inside.getnnz(axis=0)
    if background is None:
        n = len(ds.cells)-len(ids)
        so, do = total-si, all_detected-di
    elif len(background):
        other = ds.x[background]
        n, so, do = len(background), np.asarray(other.sum(axis=0)).ravel(), other.getnnz(axis=0)
    else:
        return []
    effect = si/len(ids) - so/n
    rank = np.lexsort((np.arange(len(effect)), -abs(effect)))
    return [dict(gene_id=ds.genes[j]["gene_id"], symbol=ds.genes[j].get("symbol",""),
                 effect=float(effect[j]), pct_in=di[j]/len(ids), pct_out=do[j]/n)
            for j in rank[:limit] if abs(effect[j])>1e-12]


class Pool:
    def __init__(self, cells, barcodes=None):
        self.cells = cells
        self.barcodes = dict(zip(cells, barcodes or cells))
        self.candidates = {}
        self.seeds = defaultdict(list)
        self.background = defaultdict(set)

    def add(self, idx, source):
        members = sorted(self.cells[i] for i in idx)
        if not members:
            raise ValueError("empty initial candidate")
        gid = "H" + digest(members)[:24]
        item = self.candidates.setdefault(gid, dict(id=gid, members=members, sources=[]))
        if item["members"] != members:
            raise RuntimeError("cell-set hash collision")
        item["sources"].append(source)
        if len(members) == 1:
            self.seeds[members[0]].append(source)
        return gid

    def rows(self, ds, limit):
        reverse = {c:i for i,c in enumerate(self.cells)}
        total = np.asarray(ds.x.sum(axis=0)).ravel()
        detection = ds.x.getnnz(axis=0)
        for gid,item in sorted(self.candidates.items()):
            idx = [reverse[c] for c in item["members"]]
            bg = None
            if len(idx)==1 and self.background.get(idx[0]):
                bg = sorted(self.background[idx[0]] - set(idx))
            item["marker_background"] = {"type": "local_neighbours", "members": [self.cells[i] for i in bg]} if bg else {"type":"complement"}
            item["marker_kind"] = "exploratory_single_cell" if len(idx)==1 else "group_effect"
            if ds.x[idx].nnz == 0:
                item["marker_status"] = "no_expression"
                item["markers"] = []
            else:
                item["markers"] = marker_cards(ds,idx,bg,total,detection,limit)
                item["marker_status"] = "ok" if item["markers"] else "no_contrast"
            yield item

    def seed_rows(self):
        for cell, sources in sorted(self.seeds.items()):
            yield {"id":"S"+digest(cell)[:24], "cell_id":cell, "candidate_id":"H"+digest([cell])[:24],
                   "barcode":self.barcodes[cell],"status":"pending","sources":sources}


def add_structure(pool, name, valid, unit, nbrs, radius, cfg):
    if not len(valid) or nbrs.shape[1]==0:
        return []
    med = np.median(radius)
    bound = float(med + cfg["outlier_mad"] * 1.4826 * np.median(abs(radius-med)) + 1e-8)
    flagged = radius > bound
    for i in np.flatnonzero(flagged):
        j = int(valid[i])
        pool.add([j], {"kind":"outlier","view":name,"kth_distance":float(radius[i]),"threshold":bound})
        pool.background[j].update(valid[nbrs[i]].tolist())
    pairs = []
    for i, row in enumerate(nbrs):
        j = row[0]
        if j > i and nbrs[j,0] == i:
            similarity = float(np.dot(unit[i], unit[j]))
            if similarity >= cfg["local_cosine"]:
                pairs.append((-similarity, int(valid[i]), int(valid[j])))
    for inverse_score,i,j in sorted(pairs)[:cfg["local_pairs"]]:
        pool.add([i,j], {"kind":"reciprocal_pair","view":name,"cosine":-inverse_score})
    return [{"cell_id":pool.cells[int(row)], "kth_distance":float(radius[i]),
             "threshold":bound, "outlier":int(flagged[i])} for i,row in enumerate(valid)]


def record_partitions(pool, name, valid, nbrs, partitions, data_version):
    rows = []
    excluded = sorted(set(range(len(pool.cells))) - set(valid.tolist()))
    for res, values in partitions.items():
        groups = defaultdict(list)
        for i, value in enumerate(values):
            # Seurat can use the *same* singleton label for unrelated cells.
            key = ("singleton", int(valid[i])) if str(value)=="singleton" else ("cluster", str(value))
            groups[key].append(int(valid[i]))
        for i in excluded:
            groups[("uninformative",i)] = [i]
        partition_id = "P" + digest([data_version,name,res,sorted(sorted(pool.cells[i] for i in g)
                                                    for g in groups.values())])[:24]
        aindex = {int(v):i for i,v in enumerate(valid)}
        for key, idx in groups.items():
            raw = key[0] if key[0]!="cluster" else key[1]
            src = {"kind":"partition","view":name,"resolution":res,"partition_id":partition_id,"raw_cluster":raw}
            if name=="V_HVG":
                src["marker_file"]="seurat/markers.tsv"
            gid = pool.add(idx,src)
            src["group_id"] = gid
            if len(idx)==1 and idx[0] in aindex:
                pool.background[idx[0]].update(valid[nbrs[aindex[idx[0]]]].tolist())
            rows.extend(dict(view=name,resolution=res,cell_id=pool.cells[i],group_id=gid) for i in idx)
    return rows


def structural_statistics(pool, view_dir, names, stages):
    """Precompute observed within-group SNN incidence from each frozen graph."""
    for record in pool.candidates.values():
        record["structure"] = {}
    for view in names:
        if stages.get(view+"_graph", {}).get("status") != "ok":
            for record in pool.candidates.values():
                record["structure"][view] = {"status": "graph_unavailable"}
            continue
        with (view_dir/(view+"_graph_rows.tsv")).open() as stream:
            ids = [row["cell_id"] for row in csv.DictReader(stream, delimiter="\t")]
        index = {cell:i for i,cell in enumerate(ids)}
        graph = sparse.load_npz(view_dir/(view+"_snn.npz")).tocsr()
        if graph.shape != (len(ids), len(ids)) or len(index) != len(ids):
            raise ValueError("Invalid SNN row identifiers for "+view)
        if not np.isfinite(graph.data).all() or np.any(graph.data < 0):
            raise ValueError("Invalid SNN edge weight for "+view)
        graph.sum_duplicates()
        graph.setdiag(0)
        graph.eliminate_zeros()
        degrees = np.asarray(graph.sum(axis=1)).ravel()
        stamps = np.full(len(ids), -1, dtype=np.int64)
        for stamp, record in enumerate(pool.candidates.values()):
            members = record["members"]
            represented = [index[c] for c in members if c in index]
            stamps[represented] = stamp
            within = 0.0
            for row in represented:
                start, stop = graph.indptr[row:row+2]
                idx = graph.indices[start:stop]
                within += graph.data[start:stop][stamps[idx] == stamp].sum()
            incident = float(degrees[represented].sum())
            record["structure"][view] = {
                "status": "ok", "n_graph_cells": len(represented),
                "graph_coverage": len(represented)/len(members),
                "n_zero_degree": int((degrees[represented] == 0).sum()),
                "within_weight": float(within), "incident_weight": incident,
                "within_fraction": float(within/incident) if incident > 0 else None,
            }


def write_jsonl(path, records):
    with path.open("w",encoding="utf8") as f:
        for item in records:
            f.write(canonical(item)+"\n")


def seurat_markers(ds, rows, cfg, work, dest):
    active = np.flatnonzero(np.diff(ds.x.indptr))
    lookup = {ds.cells[j]["cell_id"]:i+1 for i,j in enumerate(active)}
    write_tsv(work/"marker_partitions.tsv",["index","resolution","group_id"],
              ({"index":lookup[r["cell_id"]], "resolution":r["resolution"], "group_id":r["group_id"]}
               for r in rows if r["view"]=="V_HVG" and r["cell_id"] in lookup))
    run_r("markers", work, cfg)
    def converted():
        with (work/"markers_out.tsv").open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                gene = ds.genes[int(row.pop("gene_index"))-1]
                yield dict(**row,gene_id=gene["gene_id"],symbol=gene.get("symbol",""))
    write_tsv(dest,MARKER_FIELDS,converted())


def build(dataset, out, cfg, work):
    start = time.monotonic()
    start_cpu = os.times()
    ds = read_dataset(dataset, cfg["target_sum"])
    p = out/"views"
    p.mkdir()
    (out/"seurat").mkdir()
    sparse.save_npz(out/"data.npz", ds.x)
    if ds.raw is not None:
        sparse.save_npz(out/"raw.npz",ds.raw)
    write_tsv(out/"cells.tsv", ds.cells[0].keys(), ds.cells)
    write_tsv(out/"genes.tsv", ds.genes[0].keys(), ds.genes)
    base = ["cells.tsv","genes.tsv","data.npz"] + (["raw.npz"] if ds.raw is not None else [])
    version = digest([ds.info, {k:file_sha(out/k) for k in base}])
    valid_counts = np.diff(ds.x.indptr)>0
    pool = Pool(ds.cell_ids, [row.get("barcode") or row["cell_id"] for row in ds.cells])
    rows, done, stages = [], [], {}
    for i in np.flatnonzero(~valid_counts):
        pool.add([i], {"kind":"zero_expression","reason":"no detectable expression"})
    write_tsv(out/"cell_qc.tsv",["cell_id","zero_expression"],
              ({"cell_id":c,"zero_expression":int(not valid_counts[i])} for i,c in enumerate(pool.cells)))
    for name in NAMES:
        path = p/(name+".npy")
        if name not in cfg["views"]:
            stages[name]={"status":"disabled"}
            continue
        try:
            if name=="V_HVG":
                z = seurat_pca(ds,cfg,work,path,out/"seurat")
            elif name=="V_GINI":
                z = gini_pca(ds,cfg,path,p/"gini_features.tsv")
            else:
                z = projection(ds.x,name,cfg,path)
            done.append(name)
            stages[name]={"status":"ok","shape":list(z.shape)}
            if name == "V_HVG" and ds.info["layer"] == "linear_expression":
                stages[name]["selection_method"] = "dispersion"
        except FAILURES as exc:
            stages[name] = blocked(exc)
            continue
        valid, unit, nbrs, radius = neighbours(z,valid_counts,cfg)
        write_tsv(p/(name+"_graph_rows.tsv"),["cell_id"],({"cell_id":pool.cells[i]} for i in valid))
        write_tsv(p/(name+"_outliers.tsv"),["cell_id","kth_distance","threshold","outlier"],
                  add_structure(pool,name,valid,unit,nbrs,radius,cfg))
        try:
            dest = p/(name+"_snn.npz")
            if len(valid)<2:
                sparse.save_npz(dest,sparse.csr_matrix((len(valid),len(valid))))
                partitions={r:list(range(len(valid))) for r in cfg["resolutions"]}
            elif cfg["graph_backend"]=="igraph":
                partitions=python_graph(nbrs,cfg,dest)
            else:
                partitions=seurat_graph(unit,cfg,work,dest)
            rows.extend(record_partitions(pool,name,valid,nbrs,partitions,version))
            stages[name+"_graph"]={"status":"ok","informative":len(valid),"degenerate":len(valid)<2}
        except FAILURES as exc:
            stages[name+"_graph"]=blocked(exc)
    # If no graph was available, unknown membership remains individually traceable.
    covered = {c for v in pool.candidates.values() for c in v["members"]}
    for i,cell in enumerate(pool.cells):
        if cell not in covered:
            pool.add([i],{"kind":"unassigned_input","reason":"no available initial group"})
    structural_statistics(pool,p,cfg["views"],stages)
    markerfile = out/"seurat"/"markers.tsv"
    if any(r["view"]=="V_HVG" for r in rows):
        try:
            seurat_markers(ds,rows,cfg,work,markerfile)
            stages["seurat_markers"]={"status":"ok"}
        except FAILURES as exc:
            stages["seurat_markers"]=blocked(exc)
    if not markerfile.exists():
        write_tsv(markerfile,MARKER_FIELDS,[])
    write_tsv(out/"baselines.tsv",["view","resolution","cell_id","group_id"],rows)
    write_jsonl(out/"initial.jsonl", pool.rows(ds,cfg["markers"]))
    write_jsonl(out/"seeds.jsonl", pool.seed_rows())
    files = {f.relative_to(out).as_posix():file_sha(f) for f in sorted(out.rglob("*")) if f.is_file()}
    packages = {k:importlib.metadata.version(k) for k in ("numpy","scipy","scikit-learn")}
    if cfg["graph_backend"]=="igraph":
        try:
            packages["igraph"]=importlib.metadata.version("igraph")
        except importlib.metadata.PackageNotFoundError:
            pass
    # Costs describe this invocation. Immutable scientific files also have their
    # own hash, independent of timing, for equality checks across preparations.
    cpu = os.times()
    runtime = {"elapsed_seconds":time.monotonic()-start,
               "python_cpu_seconds":cpu.user+cpu.system-start_cpu.user-start_cpu.system,
               "subprocess_cpu_seconds":cpu.children_user+cpu.children_system-start_cpu.children_user-start_cpu.children_system}
    module = Path(__file__).parent
    sources = {key:file_sha(module/key) for key in ("prepare.py","data.py","views.py","r/seurat.R")}
    manifest={"schema":"rarecell.snapshot.v1","method_version":"prepare-1.0.0","data_version":version,
              "views":done,"n_cells":ds.x.shape[0],"n_genes":ds.x.shape[1],"nnz":ds.x.nnz,
              "zero_expression":int(sum(~valid_counts)),"layer":ds.info["layer"],
              "species":ds.info["species"],"tissue":ds.info["tissue"],"prepare":cfg,
              "packages":packages,"stages":stages,"files":files,"cost":runtime,"code_sha256":sources,
              "initial_evidence_hash":digest(files),
              "status":"partial" if any(s["status"]=="blocked" for s in stages.values()) else "complete"}
    manifest["snapshot_hash"]=digest(manifest)
    (out/"manifest.json").write_text(canonical(manifest)+"\n",encoding="utf8")


def prepare(dataset, out, cfg):
    """Atomically publish a new snapshot; a published directory is immutable."""
    out = Path(out).resolve()
    out.parent.mkdir(parents=True,exist_ok=True)
    settings=config(cfg)
    with (out.parent/("."+out.name+".prepare.lock")).open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if out.exists():
            raise FileExistsError(f"snapshot destination exists: {out}")
        with tempfile.TemporaryDirectory(dir=out.parent,prefix=".prepare-") as staging:
            stage=Path(staging)
            bundle=stage/"result"
            work=stage/"work"
            bundle.mkdir()
            work.mkdir()
            build(dataset,bundle,settings,work)
            os.rename(bundle,out)
    return out


def verify_snapshot(path):
    root=Path(path).resolve()
    record=json.loads((root/"manifest.json").read_text(encoding="utf8"))
    expected=record.pop("snapshot_hash")
    if digest(record)!=expected or record.get("schema")!="rarecell.snapshot.v1":
        raise ValueError("snapshot manifest hash or schema mismatch")
    for relative, value in record["files"].items():
        f=root/relative
        if f.is_symlink() or f.resolve()!=f or root not in f.resolve().parents or file_sha(f)!=value:
            raise ValueError(f"snapshot artifact integrity failure: {relative}")
    return dict(record,snapshot_hash=expected)
