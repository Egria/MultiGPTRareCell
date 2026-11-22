# 第三方实现与环境

## 使用范围

| 来源 | 固定版本／授权 | 本版调用及边界 |
| --- | --- | --- |
| [OntoAnno](https://github.com/yrsong001/OntoAnno/tree/ceda74c7d74af67a02093b07dd87623d846a59fe) | 2.0.0，commit `ceda74c7d74af67a02093b07dd87623d846a59fe`，MIT | 已实际检出源码并核对 DESCRIPTION、NAMESPACE、R/preprocess.R、R/cluster_anno.R、R/llm.R、R/score.R。调用作者 `gptcelltype` 和 `map_celltypes_to_cl`。R 进程只处理已有 marker。 |
| [GPTAnno 历史地址](https://github.com/yrsong001/GPTAnno) | 检查时 HEAD 与上行 commit 相同 | 同一作者代码的一条旧入口。仅安装一个 R 包。方法背景参见[作者预印本](https://pmc.ncbi.nlm.nih.gov/articles/PMC12687786/)。 |
| [Seurat](https://github.com/satijalab/seurat/tree/v5.3.0) | 5.3.0，commit `9938d61b484bfc7203cda1bde0e0a8d3540d038b`，MIT | 从官方 R 实现获取 HVG、PCA、FindNeighbors 的 SNN 和 FindClusters 的 Leiden、FindMarkers 的组际 marker。代码与选择参数见 `rarecell/r/seurat.R`。 |
| [GiniClust3](https://github.com/rdong08/GiniClust3/tree/a58d0c886cde690398c2a46e4320f1247f72cd4f) | 1.1.0 源码，commit `a58d0c886cde690398c2a46e4320f1247f72cd4f`，MIT | 已核对 `giniclust3/gini.py`。作者实现用密集表达、两次 LOWESS 校正 log(max) 趋势并按阈值选择。这里采用独立稀疏 Gini 及两次分箱趋势校正实现。具体差异见 DOC_VIEWS.md。 |
| [scikit-learn](https://github.com/scikit-learn/scikit-learn/tree/1.5.2) | 1.5.2，BSD-3-Clause | 调用 [SparseRandomProjection](https://scikit-learn.org/stable/modules/generated/sklearn.random_projection.SparseRandomProjection.html) 和 [PolynomialCountSketch](https://scikit-learn.org/stable/modules/generated/sklearn.kernel_approximation.PolynomialCountSketch.html)。实际执行的版本由 wheel 和 snapshot 中的 packages 记录。 |
| [igraph](https://github.com/igraph/python-igraph) | Python 0.11.8，GPL-2.0-or-later | core 使用实际 `community_leiden`。R 和 Python 后端分别记录图方法。 |

[Seurat 官方函数说明](https://satijalab.org/seurat/reference/findclusters) 和 v5.3.0 的 `R/clustering.R` 已核对。后者有串行和并行单例分组路径，使用固定版本与串行运行。关闭 `group.singletons` 仍可能产生共享字符串 `singleton`，桥接层依据原始 cell_id 拆开。Seurat 自身没有提供本版定义的全部 outlier 规则。

在 R 4.3.1 和 irlba 2.4.1 上，Seurat 的默认 PCA 调用可能在 `oknum(scale)` 处报 `LENGTH or similar applied to NULL object`。桥接脚本对经过 `ScaleData` 的矩阵调用 `RunPCA(..., scale=FALSE, shift=FALSE)`，显式关闭 irlba 的额外缩放和平移。成功的快照在 `seurat/versions.tsv` 记录 irlba 版本。更新前已在 PCA 阶段失败的任务需要重新准备快照，GUI 可用原数据和配置新建任务。

## OntoAnno 对接细节

作者的 `run_multi_resolution_clustering` 写 `markers_res_<resolution>.rds`，并把每个 resolution 的标签放入 `cluster_res.<resolution>`。作者高层 `ontoanno` 按同一约定读取，并执行 resolution 选择。本适配器直接调用已导出的注释子函数，按 prepare 产生的精确 partition_id、group_id 和成员集合匹配 `seurat/markers.tsv`。HVG 名单及原始划分留在原 snapshot 中。实际注释对象为各个 HVG 的 Seurat partition；其他视角、singleton 和无阳性 marker 的种子继续参加探索。

同组按调整后 p 值升序、log2FC 降序和稳定 gene_id 排序，选效应正的前 `top_genes` 个不同 symbol（缺失时使用 gene_id）。作者的命名 list 分支在源码中不会执行 topgenenumber 截断，因此在适配层先截断。按 resolution、稳定候选 ID 排序。默认每批只含一组，减少串行配对和跨群提示依赖。可设 `annotation.batch_size` 为 1–50 以降低 API 开销，整个批次内容及顺序计入缓存键。支持 stop 序列的模型可显式设置 `annotation.single_line_stop: true`（要求 `batch_size: 1`），通过 ellmer 告知生成器在换行时结束一条标签；默认关闭。传 `add_cl_prompt=TRUE`，允许父层所有细胞类型及 mixed/unknown。使用 `restrict_to=NULL, parent_celltype=NULL, min_log2fc=NULL`。作者的包内 `gptcelltype` 只在此唯一注释路径调用。支持作者内置的本地 CL 同义词键对齐，未映射时记录空 ontology_id。

`rarecell/r/annotation.R` 将作者注释函数的一个本地副本绑定到私有环境。这个环境把一次真实 `OntoAnno::call_llm` 的输入和原始输出来回记录，然后返回原始文本给原函数的解析器。R 包自身无需补丁。对响应行数及名称顺序再校验，才能接受标签关联。原始响应、归一化行、选中 genes、marker 和 cell set hash、完整模型/作者配置及 contract 版本进入外部注释目录。作者传输层对 429 有最多 3 次额外重试，本包只报告可观察到的已完成批次，未知 token/价格保持空。

同一组注释可由 A1 和 A3 只读复用。已落盘的 response.json 先校验并重用；只有 dispatch.json 的批次标记远端状态未知。操作者核对服务记录后可设置 `annotation.retry_unknown: true` 进行最多一次额外尝试。多个注释进程对同一目录使用文件锁。更换设置或 snapshot 时使用新的输出目录。修改可变模型别名指向的底层权重时请显式创建新的模型标识并保留 checkpoint。

## 用户独立 R 环境

Linux HPC 可使用 micromamba 创建私有环境。先进入新环境，再设置 R_LIBS_USER。下面锁定影响桥接的作者版本；R 的编译依赖、工具链及系统 BLAS 需根据节点选择，实际测试环境使用 R 4.5.3、Seurat 5.3.0、SeuratObject 5.2.0、ellmer 0.5.2、ontologyIndex 2.12、openai 0.4.1、leidenbase 0.1.36、RANN 2.6.2 和 remotes 2.5.0。SeuratObject 5.4.0 与该 Seurat 版本的 Wilcoxon marker 路径不兼容。

```bash
micromamba create -y -n rarecell-r -c conda-forge r-base=4.5.3 r-seurat=5.3.0 r-remotes=2.5.0 r-leidenbase=0.1.36 r-rann=2.6.2 make compilers
micromamba activate rarecell-r
export R_LIBS_USER="$CONDA_PREFIX/local-R-library"
mkdir -p "$R_LIBS_USER"
Rscript --vanilla -e 'options(repos=c(CRAN="https://cloud.r-project.org")); remotes::install_version("SeuratObject",version="5.2.0",upgrade="never"); remotes::install_version("ellmer",version="0.5.2",upgrade="never"); remotes::install_version("ontologyIndex",version="2.12",upgrade="never"); remotes::install_version("openai",version="0.4.1",upgrade="never"); remotes::install_github("yrsong001/OntoAnno@ceda74c7d74af67a02093b07dd87623d846a59fe",upgrade="never",build_vignettes=FALSE)'
Rscript --vanilla -e 'library(OntoAnno); print(packageDescription("OntoAnno")$RemoteSha)'
```

仅安装 CRAN 的同名未来版本无法满足固定 commit 校验。重启作业时保留私有 `R_LIBS_USER` 及 Rscript 的 PATH；也可配置 `prepare.rscript`、`annotation.rscript` 为私有解释器路径。依赖缺失可用 doctor 定位。

## 模型选项

发现接口 `provider.kind: openai_compatible` 支持非流式 Chat Completions 和 Responses。配置 `base_url`、`endpoint`、`api_style`（chat/responses）、`model`、`timeout` 和 `api_key_env`。凭据只通过环境变量传入，本地无鉴权服务可设为 null。`options` 可显式传 temperature、seed、top_p、max_tokens、max_completion_tokens、max_output_tokens、reasoning_effort、reasoning。只填写具体部署实际支持的字段。URL 使用 TLS 或本机回环。未完整结束的模型回复不会送去执行。网络中断的远端执行可能已经计费。

注释接口由作者的 ellmer 传输提供。`annotation.provider` 可以是 openai、anthropic、gemini、ollama 或 vllm。设置独立的 `model`、`api_key_env`、`timeout` 和 `top_genes`。前三个配置相应远端密钥的环境变量名。Ollama 使用 `api_key_env: null`。vLLM 的 ellmer 0.5.2 客户端要求 `VLLM_API_KEY` 存在，包括不检查鉴权的本地部署。可将 `api_key_env` 设为其他凭据环境变量名，R 子进程会传给作者客户端。两个本地接口都可设 `base_url`（例如兼容 vLLM 的 `http://127.0.0.1:8000/v1`）。固定作者版本会去掉 Ollama 的冒号 tag，因此适配器要求使用不含 tag 的唯一模型别名。锁定底层模型权重及部署版本有助于重放。

## 重建公开 PBMC 接入样例

实际 smoke 数据来自 **SeuratObject 5.2.0** 自带的 `pbmc_small`。仅抽取 RNA counts、基因名和 barcode，维度为 80 细胞 × 230 基因，4,456 个非零值。对象已经过小样本和特征筛选，无法代表完整 PBMC-TCR 队列。没有使用对象的聚类、细胞标签、图、PCA 或预处理 marker。此样例用来检查真实输入、作者接口和模型传输，未定义 target 或计算 benchmark F1。

在已安装 Python 包和私有 R 环境的 shell 中运行以下命令。它们输出符合 prepare 契约的目录。

```bash
mkdir -p pbmc_small_input
Rscript --vanilla - <<'RS'
library(SeuratObject)
data("pbmc_small", package="SeuratObject")
x <- LayerData(pbmc_small, assay="RNA", layer="counts")
dst <- "pbmc_small_input"
Matrix::writeMM(t(x), file.path(dst,"counts.mtx"))
write.table(data.frame(cell_id=colnames(x)),file.path(dst,"cells.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(gene_id=rownames(x),symbol=rownames(x)),file.path(dst,"genes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
RS
python - <<'PY'
from pathlib import Path
import yaml
from scipy import io, sparse
p = Path("pbmc_small_input")
sparse.save_npz(p/"counts.npz", io.mmread(p/"counts.mtx").tocsr())
manifest = dict(counts="counts.npz", cells="cells.tsv", genes="genes.tsv",
                orientation="cells_by_genes", layer="counts", species="human",
                tissue="peripheral blood", allowed_metadata=[])
(p/"dataset.yaml").write_text(yaml.safe_dump(manifest))
PY
```

验证采用四个表示、seed 1009、resolution 0.4、HVG/Gini 特征上限 100、8 PCs、k=8、linear 64 / TensorSketch 128。科学结果保留给独立验证。两个实测 CPU 模型通过 llama.cpp b12609 提供兼容服务，测试别名用于记录部署。权重和推理服务器均未随代码包分发。

| 实际模型文件 | SHA-256 |
| --- | --- |
| [Qwen2.5-0.5B-Instruct-Q4_K_M.gguf](https://huggingface.co/bartowski/Qwen2.5-0.5B-Instruct-GGUF) | `6eb923e7d26e9cea28811e1a8e852009b21242fb157b26149d3b188f3a8c8653` |
| [Qwen2.5-3B-Instruct-Q4_K_M.gguf](https://huggingface.co/bartowski/Qwen2.5-3B-Instruct-GGUF) | `9c9f56a391a3abbd5b89d0245bf6106081bcc3173119d4229235dd9d23253f94` |

公开小样例的 5 个可注释 partition 群各提交一次真实请求。3B 模型自由生成的首次运行只有 2 个群得到可用响应。启用单行 stop 和单群请求后，5/5 个响应都能核对群身份及输出格式，4 个是混合名称，1 个是 Unknown。部分模型命名含可疑缩写，不能当作 GT。显式 partial 读取仅接受持久记录中有完整哈希和对应成功响应的群；缺失的群保留未注释状态。
