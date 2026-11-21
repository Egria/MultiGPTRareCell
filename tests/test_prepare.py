import csv
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy import sparse
from sklearn.kernel_approximation import PolynomialCountSketch
from sklearn.preprocessing import normalize

from rarecell.data import read_dataset, write_tsv
from rarecell.prepare import NAMES, Pool, prepare, record_partitions, verify_snapshot
from rarecell.views import gini_features

EXAMPLE = Path(__file__).resolve().parents[1]/"examples"/"minimal_sparse"/"input"


def fixture(tmp_path, a, ids=None, genes=None, spec=None, meta=None):
    tmp_path.mkdir(parents=True)
    ids = ids or [f"C{i:02d}" for i in range(a.shape[0])]
    genes = genes or [f"G{i:02d}" for i in range(a.shape[1])]
    sparse.save_npz(tmp_path/"counts.npz", sparse.csr_matrix(a))
    meta = meta or [{} for _ in ids]
    rows = [dict(r, cell_id=c) for c,r in zip(ids,meta)]
    write_tsv(tmp_path/"cells.tsv",list(rows[0]),rows)
    write_tsv(tmp_path/"genes.tsv",["gene_id"],[dict(gene_id=g) for g in genes])
    cfg = dict(counts="counts.npz",cells="cells.tsv",genes="genes.tsv",orientation="cells_by_genes",
               layer="counts",species="human",tissue="test",allowed_metadata=[])
    (tmp_path/"dataset.yaml").write_text(yaml.safe_dump(dict(cfg,**(spec or {}))))
    return tmp_path/"dataset.yaml"


def make_input(spec, tmp_path):
    return fixture(tmp_path/"d", np.array([[1,3,0],[2,6,0],[0,0,0]]),spec=spec)


def test_depth_control_and_log_layer(tmp_path):
    path = make_input({},tmp_path)
    ds = read_dataset(path)
    assert np.issubdtype(ds.raw.dtype, np.integer)
    np.testing.assert_equal(ds.raw.toarray(), [[1,3,0],[2,6,0],[0,0,0]])
    np.testing.assert_allclose(ds.x.toarray()[0], np.log1p([2500,7500,0]))
    np.testing.assert_equal(ds.x.toarray()[0],ds.x.toarray()[1])
    assert ds.x[2].nnz==0
    p=fixture(tmp_path/"log",ds.x,spec=dict(layer="log1p_normalized"))
    log=read_dataset(p)
    assert log.raw is None
    np.testing.assert_equal(ds.x.toarray(),log.x.toarray())
    sn=prepare(path,tmp_path/"prepared",{"prepare":{"views":["V_LINEAR","V_TENSOR2"],
               "graph_backend":"igraph","k":1,"linear_dim":8,"tensor_dim":32}})
    empty = next(c for c in map(json.loads,(sn/"initial.jsonl").read_text().splitlines()) if c["members"]==["C02"])
    assert empty["markers"]==[] and empty["marker_status"]=="no_expression"
    rows=list(csv.DictReader((sn/"baselines.tsv").open(), delimiter="\t"))
    for view in ("V_LINEAR","V_TENSOR2"):
        z=np.load(sn/"views"/(view+".npy"))
        np.testing.assert_equal(z[0],z[1])
        for res in {r["resolution"] for r in rows}:
            group={r["cell_id"]:r["group_id"] for r in rows if r["view"]==view and r["resolution"]==res}
            assert group["C00"]==group["C01"] and group["C02"]!=group["C00"]


@pytest.mark.parametrize("spec",[
    {"orientation":"genes_by_cells"}, {"layer":"scaled"}, {"allowed_metadata":["GT"]},
    {"gt":"../gt.tsv"}, {"counts":"../external.npz"}
])
def test_reject_wrong_schema(tmp_path,spec):
    with pytest.raises(ValueError):
        read_dataset(make_input(spec,tmp_path))


@pytest.mark.parametrize("mutation",["shape","duplicate","infinite","negative","nan","csc"])
def test_sparse_input_validation(tmp_path,mutation):
    p=make_input({},tmp_path)
    a=sparse.load_npz(p.parent/"counts.npz").astype(float)
    if mutation=="shape":
        a=a[:2]
    elif mutation=="duplicate":
        (p.parent/"genes.tsv").write_text("gene_id\nA\nA\nC\n")
    elif mutation=="csc":
        a=a.tocsc()
    else:
        a.data[0]={"infinite":np.inf,"negative":-1.,"nan":np.nan}[mutation]
    sparse.save_npz(p.parent/"counts.npz",a)
    with pytest.raises(ValueError):
        read_dataset(p)


def test_sparse_gini_formula_and_detrending():
    rng=np.random.RandomState(47)
    a=rng.poisson(rng.uniform(.2,3.0, size=(25,180))).astype(float)
    a[:,0]=0
    a[:,1]=4
    a[3,2]=30
    ids=[f"G{i:04d}" for i in range(a.shape[1])]
    selected,g,det,fit=gini_features(sparse.csr_matrix(a),ids)
    for j in range(2, a.shape[1]):
        v=a[:,j]
        exact=np.abs(v[:,None]-v[None,:]).sum() / (2*len(v)*v.sum())
        assert g[j]==pytest.approx(exact,abs=1e-12)
    assert set(selected)==set(range(2,180))
    np.testing.assert_equal(det,(a>0).mean(axis=0))
    assert np.ptp(fit[2:])>.01
    assert selected==sorted(selected,key=lambda j:(-(g[j]-fit[j]),ids[j]))


def test_tensor_kernel_rank_and_finite_dimensions():
    rng=np.random.RandomState(61)
    rows=normalize(rng.uniform(0,1,(30,22)))
    s=rows@rows.T
    exact=s**2
    assert np.array_equal(s.argsort(),exact.argsort())
    assert np.array_equal(s.argsort(),(.25*s+.75*exact).argsort())
    errors={}
    for d in [64,2048]:
        errors[d]=[]
        for seed in [9,27,51]:
            z=PolynomialCountSketch(n_components=d,degree=2,coef0=0,gamma=1,random_state=seed).fit_transform(sparse.csr_matrix(rows))
            errors[d].append(np.sqrt(np.mean((z@z.T-exact)**2)))
    assert np.mean(errors[2048]) < np.mean(errors[64])
    assert max(errors[2048]) < .15
    print("Exact kernel RMSE by dimension and seeds 9,27,51:", errors)


def test_singleton_label_recovery_and_traceable_seeds():
    pool=Pool(["C-a","C-b","C-c","C-d"], ["AA","BA","CB","DB"])
    valid=np.array([0,1,2,3])
    nbr=np.array([[1],[0],[3],[2]])
    rows=record_partitions(pool,"V_HVG",valid,nbr,{.4:["singleton","singleton","2","2"]},"d")
    group={r["cell_id"]:r["group_id"] for r in rows}
    assert group["C-a"]!=group["C-b"] and group["C-c"]==group["C-d"]
    records=list(pool.seed_rows())
    assert len({r["id"] for r in records})==2
    assert {r["barcode"] for r in records}=={"AA","BA"}


def test_row_reordering_graph_candidates_markers_and_isolation(tmp_path):
    rng=np.random.RandomState(21)
    counts=rng.poisson(.45,(48,60))
    counts[:24,:8]+=2
    counts[24:,8:16]+=3
    counts[-1,:]=0
    ids=np.array([f"C{i:03d}" for i in range(48)])
    gs=np.array([f"G{i:03d}" for i in range(60)])
    cfg={"prepare":{"views":["V_GINI","V_LINEAR","V_TENSOR2"],"graph_backend":"igraph",
            "linear_dim":30,"tensor_dim":100,"k":4,"n_features":24,"n_pcs":6,"local_pairs":2}}
    original=None
    for i in range(2):
        order=rng.permutation(len(ids))
        genes=rng.permutation(len(gs))
        meta=[{"batch":"pre","gt_truth":f"canary_{i}_{j}"} for j in order]
        d=fixture(tmp_path/f"input{i}",counts[order][:,genes],ids[order].tolist(),gs[genes].tolist(),
                  spec={"allowed_metadata":["batch"]},meta=meta)
        s=prepare(d,tmp_path/f"snapshot{i}",cfg)
        m=verify_snapshot(s)
        assert m["status"]=="complete"
        assert "gt_truth" not in (s/"cells.tsv").read_text()
        assert "canary" not in (s/"initial.jsonl").read_text()
        records=[json.loads(x) for x in (s/"initial.jsonl").read_text().splitlines()]
        assert set(x for c in records for x in c["members"])==set(ids)
        assert len(records)==len({tuple(c["members"]) for c in records})
        assert any(len(c["sources"])>1 for c in records)
        partition=[r for r in csv.DictReader((s/"baselines.tsv").open(),delimiter="\t")]
        assert len(partition)==48*3*2
        state=(m["data_version"],(s/"initial.jsonl").read_bytes(),(s/"baselines.tsv").read_bytes())
        if original is not None:
            assert state==original
        original=state
    with (s/"seeds.jsonl").open("a") as f:
        f.write("tampered")
    with pytest.raises(ValueError,match="integrity"):
        verify_snapshot(s)


def test_missing_r_preserves_independent_preparations(tmp_path):
    cfg={"prepare":{"views":["V_HVG","V_LINEAR"],"graph_backend":"seurat","rscript":"nonexistent-rscript-123"}}
    s=prepare(EXAMPLE,tmp_path/"s",cfg)
    m=verify_snapshot(s)
    assert m["status"]=="partial" and m["stages"]["V_HVG"]["status"]=="blocked"
    assert m["stages"]["V_LINEAR"]["status"]=="ok" and m["stages"]["V_LINEAR_graph"]["status"]=="blocked"
    assert (s/"initial.jsonl").stat().st_size>0
    assert not (s/"views"/"V_HVG.npy").exists()


@pytest.mark.skipif(not shutil.which(os.environ.get("RARECELL_RSCRIPT","Rscript")),
                    reason="Seurat 5.3.0 integration requires Rscript")
def test_all_four_real_seurat_partitions_markers(tmp_path):
    cfg={"prepare":{"views":list(NAMES),"rscript":os.environ.get("RARECELL_RSCRIPT","Rscript"),
                   "n_features":25,"n_pcs":8,"k":7,"linear_dim":32,"tensor_dim":64}}
    s=prepare(EXAMPLE,tmp_path/"production",cfg)
    m=verify_snapshot(s)
    assert m["status"]=="complete",m["stages"]
    marker=list(csv.DictReader((s/"seurat"/"markers.tsv").open(),delimiter="\t"))
    groups=list(csv.DictReader((s/"baselines.tsv").open(),delimiter="\t"))
    hvg={(r["resolution"],r["group_id"]) for r in groups if r["view"]=="V_HVG"}
    assert len(groups)==160*4*2 and len(marker)>10
    assert any(float(r["avg_log2FC"])>0 and float(r["p_val_adj"]) < .05 for r in marker)
    assert {(r["resolution"],r["group_id"]) for r in marker}<=hvg
    assert all(0<=float(r["pct.1"])<=1 for r in marker)
    assert min(sparse.load_npz(s/"views"/(v+"_snn.npz")).nnz for v in NAMES)>160
    print("Seurat production fixture:",len(groups),"assignments;",len(marker),"positive marker rows")


def test_uniform_null_no_forced_outlier_or_confirmation(tmp_path):
    counts = np.tile([1,4,2,3,0], (30,1))
    path = fixture(tmp_path/"same",counts)
    s = prepare(path,tmp_path/"snapshot",{"prepare":{"views":["V_LINEAR"],
        "graph_backend":"igraph","linear_dim":5,"k":5,"outlier_mad":3.5}})
    outliers=list(csv.DictReader((s/"views"/"V_LINEAR_outliers.tsv").open(),delimiter="\t"))
    assert len(outliers)==30 and all(row["outlier"]=="0" for row in outliers)
    rows=[json.loads(line) for line in (s/"initial.jsonl").read_text().splitlines()]
    assert all(r.get("status")!="confirmed" and not r["markers"] for r in rows)
    assert all(source["kind"]!="outlier" for r in rows for source in r["sources"])


def test_structural_incidence_uses_exact_graph_rows_and_ignores_self_edges(tmp_path):
    from rarecell.prepare import structural_statistics
    pool=Pool(["A","B","C","D","absent"])
    select=pool.add([0,1,3,4],{"kind":"example"})
    iso=pool.add([3],{"kind":"example"})
    omitted=pool.add([4],{"kind":"example"})
    row_ids=["B","D","A","C"]
    write_tsv(tmp_path/"V_LINEAR_graph_rows.tsv",["cell_id"],[dict(cell_id=x) for x in row_ids])
    weights=np.array([[99,0,2,1],[0,99,0,0],[2,0,99,0],[1,0,0,99]],dtype=float)
    sparse.save_npz(tmp_path/"V_LINEAR_snn.npz",sparse.csr_matrix(weights))
    structural_statistics(pool,tmp_path,["V_LINEAR","V_HVG"],{"V_LINEAR_graph":{"status":"ok"}})
    group=pool.candidates[select]["structure"]
    assert group["V_LINEAR"]=={"status":"ok","n_graph_cells":3,"graph_coverage":.75,
                              "n_zero_degree":1,"within_weight":4,"incident_weight":5,"within_fraction":.8}
    assert group["V_HVG"]=={"status":"graph_unavailable"}
    assert pool.candidates[iso]["structure"]["V_LINEAR"]["within_fraction"] is None
    assert pool.candidates[omitted]["structure"]["V_LINEAR"]["graph_coverage"]==0


def test_duplicate_integer_entries_preserve_counts_and_reject_overflow(tmp_path):
    path = fixture(tmp_path / "d", np.array([[0, 1]]))
    def write(values, dtype):
        mat = sparse.csr_matrix((np.array(values,dtype=dtype), [0,0,1], [0,3]), shape=(1,2))
        sparse.save_npz(path.parent/"counts.npz", mat)
    write([230, 220, 5], np.uint8)
    data = read_dataset(path)
    np.testing.assert_equal(data.raw.toarray(), [[450,5]])
    np.testing.assert_allclose(data.x.toarray(), np.log1p(np.array([[450,5]]) / 455 * 10000))
    write([2**63, 2**63, 1], np.uint64)
    with pytest.raises(ValueError, match="storage range"):
        read_dataset(path)
