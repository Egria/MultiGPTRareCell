# 在 DCC 使用 single-agent 的基准数据

## 导入以前的运行数据

从 MultiGPTRareCell 根目录执行。默认读取旁边的 LLMRareCell 工程，检查旧 release.lock.json 中每个被使用文件的 SHA256，核对矩阵行、源视图 cell_id 顺序及参考覆盖，再为每个 view 写一对 GUI 路径。已有完整导出会校验后复用。各数据集独立处理，某个失败时会继续处理其余数据集并返回非零退出码。

```bash
cd /hpc/group/jilab/xn24/MultiGPTRareCell
source .venv/bin/activate
python -m pip install -e .
export SINGLE_AGENT_PROJECT=/hpc/group/jilab/xn24/LLMRareCell
bash scripts/assemble_single_agent.sh
```

只准备一个数据集可以执行 `bash scripts/assemble_single_agent.sh cellsius`。标准输出每行依次给出数据集清单和评估 GT。GUI 新建任务时分别填写这两项。默认输出是 `work/single-agent-datasets/<dataset>/<view>/input/dataset.yaml` 与同一 view 下的 `evaluation_fixture/gt.tsv`。源码中 examples/datasets 的 YAML 是模板；真实导出的目录还包括非负稀疏矩阵、行列标识、独立评估文件和逐文件哈希。下载版 ZIP 含相同目录结构，解压到 `work/single-agent-datasets/` 即可直接运行。

| 用户名称 | 导出目录 | 默认 view | single-agent 来源和参考 |
|---|---|---|---|
| pbmc_tcr | pbmc_tcr | view0 至 view4 | 10x sc5p_v2_hs_PBMC_10k，10,548 个 GEX 细胞起步。原始 D4_VJ_1 顺序保留 144、72、36、18、9 个 MAIT 阳性；全部其他细胞保留。GT 含 MAIT 和 iNKT 的四态参考。 |
| cellsius | cellsius | V000 至 V005 | Zenodo 3238275 原版 replicate 0。前四个单目标条件为 200、202、205、210 个细胞，V004 是 1,000 个细胞的多目标混合，V005 是 1,000 个细胞的技术零对照。 |
| cellsius_spheroid | cellsius_spheroid | V000 | PRJNA545246 人皮层球状体。复用已冻结 D10 或原 standalone 导出；portable profile 按 G.sub_1 的 Choroid_plexus 目标评分。 |
| villani | villani | V000 | GSE94820/SCP43 的 1,078 个发现集细胞，原 TPM 和作者细标签。 |
| baron | baron | V000 | GSE84133 的 4 份人胰腺样品，8,569 个细胞、14 类作者标签。 |
| gse156490 | gse156490 | V000 | 7,393 个胰腺细胞，RNA 来自未整合表达矩阵，loom 提供作者标签。 |
| pgc | human_pgc | V000 | E-MTAB-9388，1,195 个细胞，包括 7 个作者标注的 PGC。 |
| airway | airway | V000 | GSE102580 小鼠 Fig. 2 气道上皮，14,163 个细胞，使用原始 counts。 |
| pbmc68k | pbmc68k | V000 | 10x fresh_68k_pbmc_donor_a，68,579 个细胞，使用已发表的参考相关性标签。 |

`SINGLE_AGENT_PROJECT` 指完整的旧工程根目录，包含 data/ 和 benchmark/。也可直接执行下列命令，读取准备好的工程，无需安装原 preparer。

```bash
python -m rarecell.single_agent_data --project "$SINGLE_AGENT_PROJECT" \
  --out work/single-agent-datasets pbmc_tcr cellsius cellsius_spheroid villani baron gse156490 pgc airway pbmc68k
```

PBMC-TCR 的五个 view 直接读取原版 row_indices.npy，保留旧运行用的细胞 ID。只有 V000 的旧 PBMC-TCR release 需要在新目录以固定版本重建。CellSIUS 原版以每个条件各自的矩阵为输入。准备过额外重复并需要全部导出时，设置 `INCLUDE_EXTRA_VIEWS=1` 和新的 `DATASET_OUT`；CellSIUS 的完整额外方案共有 60 个条件。已发布的 main-only release 如需新增视图，使用全新的 SINGLE_AGENT_PROJECT。

## 在线组装缺少的 release

在计算节点上运行以下命令。已有旧工程 .venv 可直接作为 PREP_PYTHON，需具备原 preparer 的依赖。全新环境可使用以下独立环境，随后用相同的脚本按数据集组装。脚本从固定的源代码版本构建、按源数据锁校验实际下载并复用缓存。旧仓库需要你已有的 GitHub 读取权限。

```bash
cd /hpc/group/jilab/xn24/MultiGPTRareCell
export SINGLE_AGENT_CODE="$PWD/work/single-agent-preparer"
export SINGLE_AGENT_PROJECT="$PWD/work/legacy-releases"
mkdir -p work
git clone https://github.com/Egria/LLMRareCell.git "$SINGLE_AGENT_CODE"
git -C "$SINGLE_AGENT_CODE" checkout --detach eed0c7a2d07f6772fbc85c0cffea2e56244c30e2
python3 -m venv work/legacy-prepare-venv
export PREP_PYTHON="$PWD/work/legacy-prepare-venv/bin/python"
"$PREP_PYTHON" -m pip install -r "$SINGLE_AGENT_CODE/module/preparation/requirements.txt"
bash scripts/assemble_single_agent.sh pbmc_tcr cellsius villani baron gse156490 pgc airway pbmc68k
```

若 D10 已在旧 LLMRareCell 工程发布，直接按第一节导入。如果只保留了 standalone 的 data/spheroid 目录，其中应有 manifest.json、counts.npz、cells.tsv、genes.tsv 和 labels.tsv。指定真实目录即可导入相同的计算参考。

```bash
export SPHEROID_EXPORT=/hpc/group/jilab/xn24/LLMRareCell/data/spheroid
bash scripts/assemble_single_agent.sh cellsius_spheroid
```

若没有旧导出，可使用原仓库的 PRJNA545246 重建流程。先用以下命令列出参考、比对及 legacy R 依赖，按原文档运行重建，完成后再导入。重建的计算标签可能与旧任务的标签有差异，要复现旧得分应选同一份冻结 release。[原版重建说明](https://github.com/Egria/LLMRareCell/blob/eed0c7a2d07f6772fbc85c0cffea2e56244c30e2/module/preparation/SPHEROID_RECONSTRUCTION.md)

```bash
"$PREP_PYTHON" "$SINGLE_AGENT_CODE/module/preparation/prepare.py" prepare \
  --dataset cellsius_spheroid --projectdir "$SINGLE_AGENT_PROJECT" --spheroid-plan
```

已配置 Cell Ranger、SRA Toolkit、MCL 和原版兼容的 R 环境后，可在计算节点启动重建。下面两个环境变量分别填本机真实的 GRCh37 Cell Ranger 参考目录和 legacy Rscript 路径。已有同一数据的 filtered 10x 矩阵时，可将参考参数换为 `--spheroid-counts-dir "$FILTERED_10X"`，跳过原始测序下载与比对。

```bash
"$PREP_PYTHON" "$SINGLE_AGENT_CODE/module/preparation/prepare.py" prepare \
  --dataset cellsius_spheroid --projectdir "$SINGLE_AGENT_PROJECT" --spheroid-reconstruct \
  --spheroid-reference "${CELLRANGER_GRCH37:?set the reference path}" \
  --spheroid-rscript "${LEGACY_RSCRIPT:?set the legacy Rscript path}" \
  --spheroid-threads 8 --spheroid-memory-gb 64
bash scripts/assemble_single_agent.sh cellsius_spheroid
```

## 运行与评估

GUI 中选择 A3_adaptive、live 和可用的模型。快速试跑可从 CellSIUS V001 开始，V000 可用于检查无目标背景。命令行使用已有的 real.yaml 和独立输出目录。

```bash
set -e
view="$PWD/work/single-agent-datasets/cellsius/V001"
task=$(mktemp -d "$PWD/work/cellsius-V001.XXXXXX")
rarecell prepare --dataset "$view/input/dataset.yaml" --out "$task/snapshot" --config real.yaml
rarecell annotate --snapshot "$task/snapshot" --out "$task/annotations" --config real.yaml
rarecell discover --snapshot "$task/snapshot" --annotations "$task/annotations" --run "$task/run" \
  --mode A3_adaptive --policy live --config real.yaml || { rc=$?; [ "$rc" -eq 2 ] && [ -f "$task/run/frozen.json" ] || exit "$rc"; }
rarecell evaluate --run "$task/run" --gt "$view/evaluation_fixture/gt.tsv" --out "$task/evaluation"
```

目标参考文件由 evaluate 自动识别，GUI 使用相同入口。完整细标签走原本的一对一分区评估；PBMC-TCR 和 portable spheroid 按每个目标在可判定细胞上选择 F1 最高的最终群，同时输出各冻结基线的同口径结果。unknown/disputed 从该目标的 TP、FP、FN 分母中排除，零阳性目标的 F1 留空。这些目标参考无法提供全体细胞类型的 ARI/NMI。

Villani TPM、PGC 的分数型 read quantification 和 GSE156490 未整合的非整数表达按原精度导出为 linear_expression，逐行归一化后 log1p，再使用 Seurat dispersion 筛特征。原始计数走 VST。在 DCC 为 PBMC68k 单独设置 `prepare.scale_max_bytes: 1073741824` 或更高并申请相应内存。默认 1,000 个特征的稠密缩放层约 523 MiB，已经超过 512 MiB 默认上限，R 的完整峰值更高。大数据优先复用旧 release，以节约公共源下载、重复验证和临时存储。
