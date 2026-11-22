# 在 DCC 使用 single-agent 的基准数据

## 已经拷好了 data 目录

`rarecell-import` 直接读取旧版的 `data/<dataset>/input.json`、`views.tsv` 和选定 view 的矩阵、细胞及基因表。可以先在有数据的 DCC 计算节点只导入一个较小的阳性条件。已有环境只需更新可编辑安装，转换全程离线。

```bash
cd /hpc/group/jilab/xn24/MultiGPTRareCell
.venv/bin/python -m pip install --no-deps -e .
DATA=/hpc/group/jilab/xn24/data
.venv/bin/rarecell-import --data-root "$DATA" --list cellsius pbmc_tcr
.venv/bin/rarecell-import --data-root "$DATA" --view V001 cellsius
```

`--list` 的四列依次是数据集名、view ID、角色、细胞数。最后一条命令打印可粘贴到 GUI 的清单绝对路径。V001 有 202 个细胞。重新打开「新建任务」，在「已导入的数据集与 view」里选择 `cellsius / V001 · 仅表达`，配置文件填自己的 `real.yaml`。只验证输入时选择「准备快照」。跑真实发现时选择「完整发现流程」、`A3_adaptive` 和已经配置的全局模型。只有表达数据时，真值字段留空。

| 拷贝中的布局 | 转换方法 |
|---|---|
| `airway`、`baron`、`pbmc68k`、`cellsius_spheroid` | 使用根目录 counts 和各 view 的行索引。Airway 使用主 raw track。 |
| `pbmc_tcr/views/view0…view4` | 按旧 view 的行索引选择矩阵行，并逐行比对原 cell_id。 |
| `cellsius/views/V000…V005` | 各自读取 view 中的 counts 和 cells，共用根目录的 genes。 |
| `villani`、`human_pgc`、`gse156490` | 分别读取 TPM 或作者提供的分数表达，保留精度并标作 `linear_expression`。 |
| `hao`、`mmochi`、`larry` | 支持表达输入。前两者为人血液 counts，LARRY 为小鼠培养细胞的归一化表达。MMoCHi 可用 `--all-views` 或指定单个 view。 |
| `spheroid/` | 含 labels、10x 和 h5ad 的独立重建产物。同级已有 `cellsius_spheroid/` 时导入后者。命令别名 `spheroid` 也指向后者。仅有独立产物时见下方 `SPHEROID_EXPORT`。 |

批量准备最初的九个基准可以使用下列命令。大数据集在有足够内存和磁盘的计算节点上转换。每个共享稀疏矩阵在单次导入该数据集时只加载一次，选定 view 分别占用输出空间。转换阶段不调用模型。

```bash
.venv/bin/rarecell-import --data-root "$DATA" \
    pbmc_tcr cellsius cellsius_spheroid villani baron gse156490 human_pgc airway pbmc68k
```

默认输出在 `work/copied-inputs/<dataset>/<view>/input/`，包括 `dataset.yaml`、`expression.npz`、`cells.tsv` 和 `genes.tsv`。同级 `import.json` 记录源文件和导出文件的 SHA256。可在同一输出目录增补 view。重新导入时检查旧导出，发现内容不一致会给出错误，改用新的 `--out` 即可保留两个版本。GUI 自动列出两个默认导出目录里的 view，自定义输出路径可手动粘贴。

### 有旧 benchmark 时启用评估

如果旧工程仍在 DCC，将 `--benchmark-root` 指向包含多个数据集子目录的 **benchmark 目录本身**。可先运行 `ls "$BENCH/cellsius/release.lock.json"` 确认位置。导入会按旧锁里的逻辑路径核对 RNA 与参考，然后为该数据集的各主 view 生成带 GT 的独立导出。

```bash
BENCH=/hpc/group/jilab/xn24/LLMRareCell/benchmark
.venv/bin/rarecell-import --data-root "$DATA" --benchmark-root "$BENCH" cellsius
```

从 GUI 的「已导入的数据集与 view」选择同名的「含评估 GT」项，数据清单和同一 view 的 `evaluation_fixture/gt.tsv` 一起填入。已完成的发现也可以在「评估结果」中指定原 run 和匹配的 GT。冻结参考导出支持本页列出的全部 12 个基准。独立 `spheroid/labels.tsv` 采用源数据的标识，其导入映射走下方 standalone 流程。

### 一次转换 DCC 的全部表达和评估文件

在 MultiGPTRareCell 根目录运行下列脚本。它默认从相邻的 `../data/` 读取 RNA，从 `../LLMRareCell/benchmark/` 读取标签和冻结校验值，遍历 12 个基准在 `views.tsv` 中发布的所有 view。MMoCHi 包含 40 个额外 view，PBMC-TCR 保留五档丰度，CellSIUS 保留六个条件。给定的两份 DCC 目录清单共对应 61 个 view。`data/spheroid` 是独立源包，规范化入口为 `cellsius_spheroid`。

```bash
cd /hpc/group/jilab/xn24/MultiGPTRareCell
bash scripts/convert_benchmarks.sh
```

在登录节点可将转换直接提交到计算节点。下方是一个初始资源申请，耗时和实际峰值内存随节点和文件系统而变。Slurm 标准输出每行是一对可用的清单及默认 GT 绝对路径，标准错误包含逐数据集进度及失败原因。

```bash
cd /hpc/group/jilab/xn24/MultiGPTRareCell
mkdir -p work
sbatch --cpus-per-task=2 --mem=64G --time=06:00:00 \
  --output="$PWD/work/import-%j.tsv" --error="$PWD/work/import-%j.log" scripts/convert_benchmarks.sh
```

也可指定子集、其他目录或 Python 解释器，例如 `DATA_ROOT=/path/to/data BENCHMARK_ROOT=/path/to/benchmark bash scripts/convert_benchmarks.sh hao larry`。默认使用项目的 `.venv/bin/python`，可用 `RARECELL_PYTHON` 覆盖。命令行等价入口是 `rarecell-import --data-root ../data --benchmark-root ../LLMRareCell/benchmark --all-views all`。全程离线，各个数据集独立发布，中途遇到文件缺失或哈希冲突会继续处理其余数据集并返回非零状态。原冻结版本和视图集合完全相同的导出会经过校验后复用。已有不同源版本或不同视图集合时可用 `DATASET_OUT` 指定新目录，GUI 中手动粘贴打印的路径。

| 基准 | GUI 默认 GT | 同一 view 的其他可评分文件 |
|---|---|---|
| Hao | 作者 L2 分区，和原稀有目标的来源层级一致 | `labels_l1.tsv`、`labels_l3.tsv`，以及 `targets.tsv` 中分开记录的 source membership 和 ADT protein surrogate |
| LARRY | pDC、Ccr7_DC、Eos 的逐细胞 source membership，按每个目标报告最佳最终群 F1 | 原始终末标签只覆盖部分细胞。此导出评价目标恢复；轨迹、克隆命运概率及终末子集 ARI 属于旧项目的其他指标。 |
| MMoCHi | 作者细标签分区 | `targets.tsv` 按目标以及 source membership、FACS sort 分别评分。两类记录来自同一分选设计。 |
| PBMC-TCR、portable spheroid | 每目标的四态 membership | unknown/disputed 在各目标评分时排除。 |
| 其余基准 | 原作者完整 fine 分区 | 一对一匹配给出 ARI、NMI 及各类 F1。 |

在 GUI「已导入的数据集与 view」中选择「含评估 GT」可填入默认文件，Hao、MMoCHi 还列出额外评分选项。命令行也可为同一个 run 更换 `rarecell evaluate --gt` 的文件，每种口径使用独立的 `--out`。这些标签、蛋白证据和目标名称均放在 `evaluation_fixture/` 中；`input/` 只包含表达矩阵、cell_id、基因表和清单。

### 在本地浏览器打开 DCC 的 GUI

已有同一计算节点的 8787 转发和 GUI 时，保留当前连接，安装新版后重新打开「新建任务」即可。首次从没有热更新功能的版本升级时，在原节点重新启动 GUI，保留 SSH 转发。在新的连接中，可按 [Duke 的配置](https://oit-rc.pages.oit.duke.edu/rcsupportdocs/dcc/vscode/) 将 `dcc-cpu` 加入本机的 SSH 配置，然后在 Windows PowerShell 执行：

```powershell
ssh -o ExitOnForwardFailure=yes -L 127.0.0.1:8787:127.0.0.1:8787 dcc-cpu
```

在这条连接打开的 DCC shell 中运行：

```bash
cd /hpc/group/jilab/xn24/MultiGPTRareCell
.venv/bin/rarecell ui --workspace "$PWD/work/tasks" --port 8787
```

保持 SSH 窗口开启，把终端打印的完整 `http://127.0.0.1:8787/#…` 链接复制到本机浏览器。SSH 目标和 GUI 进程必须在同一台机器上。已用 VS Code Remote SSH 连接该计算节点时，也可在 Ports 面板转发 8787。如果代码和数据都在自己的 Linux 或 WSL2 环境，按 README 安装依赖，将 `--data-root` 改成本机目录，在同一 Linux 环境启动 GUI 并打开其本机链接即可。

## 导入完整的冻结工程

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
