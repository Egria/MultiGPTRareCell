"""Export only counts and stable IDs from the SeuratObject reference example."""
import argparse
from pathlib import Path
import subprocess
import tempfile

from scipy import io, sparse
import yaml

SOURCE = """
suppressPackageStartupMessages(library(SeuratObject))
data("pbmc_small", package="SeuratObject")
mat <- LayerData(pbmc_small, assay="RNA", layer="counts")
args <- commandArgs(trailingOnly=TRUE)
write.table(data.frame(cell_id=colnames(mat)), file.path(args[1], "cells.tsv"),
            sep="\\t", quote=FALSE, row.names=FALSE)
genes <- rownames(mat)
write.table(data.frame(gene_id=genes, symbol=genes), file.path(args[1], "genes.tsv"),
            sep="\\t", quote=FALSE, row.names=FALSE)
Matrix::writeMM(mat, file.path(args[1], "counts.mtx"))
"""


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--rscript", default="Rscript")
    args = p.parse_args()
    dest = args.out.resolve()
    if dest.exists():
        p.error("use a new output directory")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dest.parent) as work:
        stage = Path(work) / "input"
        stage.mkdir()
        subprocess.run([args.rscript, "--vanilla", "-e", SOURCE, str(stage)], check=True)
        raw = stage / "counts.mtx"
        sparse.save_npz(stage / "counts.npz", io.mmread(raw).T.tocsr())
        raw.unlink()
        manifest = dict(cells="cells.tsv", counts="counts.npz", genes="genes.tsv", orientation="cells_by_genes",
                        layer="counts", species="human", tissue="peripheral blood", allowed_metadata=[])
        (stage / "dataset.yaml").write_text(yaml.safe_dump(manifest), encoding="utf-8")
        stage.rename(dest)
    print(dest / "dataset.yaml")


if __name__ == "__main__":
    main()
