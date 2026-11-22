args <- commandArgs(trailingOnly=TRUE)
mode <- args[[1]]
work <- args[[2]]
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
if (as.character(packageVersion("Seurat")) != "5.3.0") stop("Install pinned Seurat 5.3.0")
future::plan(future::sequential)
path <- function(x) file.path(work, x)
emit <- function(x, filename, cols=FALSE) {
  write.table(x, file=path(filename), row.names=FALSE, col.names=cols, quote=FALSE, sep="\t")
}
if (mode == "pca") {
  expr <- as(readMM(path("expression.mtx")), "CsparseMatrix")
  rownames(expr) <- paste0("G", seq_len(nrow(expr)))
  colnames(expr) <- paste0("C", seq_len(ncol(expr)))
  if (args[[8]] == "linear_expression") {
    obj <- CreateSeuratObject(CreateAssay5Object(data=expr), assay="RNA")
    selection <- "dispersion"
  } else {
    obj <- CreateSeuratObject(expr, min.cells=0, min.features=0)
    obj <- NormalizeData(obj, scale.factor=as.numeric(args[6]), verbose=FALSE)
    selection <- "vst"
  }
  obj <- FindVariableFeatures(obj, selection.method=selection,
    nfeatures=min(as.integer(args[3]), nrow(obj)), verbose=FALSE)
  genes <- VariableFeatures(obj)
  # Counts of zero/constant genes can yield uninformative features. RunPCA
  # independently checks for zero variance after scaling.
  if (length(genes)*as.double(ncol(obj))*8 > as.double(args[7])) stop("ScaleData layer exceeds scale_max_bytes")
  pcs <- min(as.integer(args[4]), length(genes)-1L, ncol(obj)-1L)
  if (pcs < 1L) stop("insufficient HVGs or informative cells for PCA")
  obj <- ScaleData(obj, features=genes, verbose=FALSE)
  var_ok <- apply(LayerData(obj, layer="scale.data"), 1L, function(x) all(is.finite(x)) && sd(x)>0)
  pcs <- min(pcs, sum(var_ok)-1L)
  if (pcs < 1L) stop("insufficient variable HVGs after scaling")
  # ScaleData already prepared the matrix. irlba 2.4.1 on R 4.3.1
  # errors when checking NULL for either of these optional transformations.
  obj <- RunPCA(obj, features=names(var_ok)[var_ok], npcs=pcs, seed.use=as.integer(args[5]),
                scale=FALSE, shift=FALSE, verbose=FALSE)
  emit(Embeddings(obj, "pca")[colnames(expr), , drop=FALSE], "pca.tsv")
  emit(sub("^G", "", genes), "hvg_indices.tsv")
  loadings <- Loadings(obj[["pca"]])
  emit(data.frame(gene_index=as.integer(sub("^G","",rownames(loadings))),loadings),
       "pca_loadings.tsv", TRUE)
  packages <- c("Seurat","SeuratObject","Matrix","irlba","leidenbase","RANN")
  emit(data.frame(package=c("R",packages),version=c(as.character(getRversion()),
       vapply(packages,function(n) as.character(packageVersion(n)),""))),"versions.tsv", TRUE)
  saveRDS(obj, path("seurat.rds"), compress=FALSE)
} else if (mode == "graph") {
  z <- as.matrix(read.table(path("embedding.tsv"), header=FALSE, sep="\t"))
  rownames(z) <- paste0("C", seq_len(nrow(z)))
  # FindNeighbors includes the query cell. Request k additional neighbors.
  graph <- FindNeighbors(z, k.param=min(nrow(z), as.integer(args[3])+1L),
    nn.method="rann", nn.eps=0, compute.SNN=TRUE, prune.SNN=1/15, verbose=FALSE)$snn
  writeMM(graph, path("snn.mtx"))
  res <- as.numeric(strsplit(args[4], ",", fixed=TRUE)[[1]])
  groups <- FindClusters(graph, resolution=res, algorithm=4L, random.seed=as.integer(args[5]),
    n.iter=10L, group.singletons=FALSE, verbose=FALSE)
  emit(groups[rownames(z), , drop=FALSE], "groups.tsv")
} else if (mode == "markers") {
  obj <- readRDS(path("seurat.rds"))
  parts <- read.delim(path("marker_partitions.tsv"), stringsAsFactors=FALSE, colClasses="character")
  outfile <- path("markers_out.tsv")
  cat("resolution\tgroup_id\tgene_index\tavg_log2FC\tpct.1\tpct.2\tp_val\tp_val_adj\n", file=outfile)
  pooled <- colnames(obj)
  if (as.character(packageVersion("SeuratObject")) != "5.2.0") {
    stop("Seurat 5.3.0 Wilcoxon marker bridge requires pinned SeuratObject 5.2.0")
  }
  for (resolution in unique(parts$resolution)) {
    r <- parts[parts$resolution == resolution, , drop=FALSE]
    for (gid in unique(r$group_id)) {
      members <- paste0("C", r$index[r$group_id == gid])
      background <- setdiff(pooled, members)
      if (length(members) < 3 || length(background) < 3) next
      res <- FindMarkers(obj[["RNA"]], cells.1=members, cells.2=background, test.use="wilcox",
        min.pct=0.01, logfc.threshold=0, only.pos=TRUE, verbose=FALSE)
      if (!nrow(res)) next
      res <- res[order(res$p_val, -res$avg_log2FC, rownames(res)),,drop=FALSE]
      tab <- data.frame(resolution=resolution, group_id=gid,
        gene_index=as.integer(sub("^G", "", rownames(res))), res[,c("avg_log2FC","pct.1","pct.2","p_val","p_val_adj")])
      write.table(tab, file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    }
  }
} else stop("unknown Seurat operation")
