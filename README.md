# LLM Rare Cell V1

在完整稀疏 scRNA-seq 矩阵上发现稀有群体的成员。程序先冻结多视角聚类、差异表达、singleton 和 outlier 证据，再由一个 Discovery Agent 查询、扩展和修订候选。每个输入 cell_id 在最终 partition 中占一行。科学效果与工程验证的实际结果见 [VALIDATION_REPORT.md](VALIDATION_REPORT.md)。

## 安装与第一次运行

适用于 Linux CPU 和 Slurm 计算节点。已验证 Python 3.9.25 和 3.12.14。项目在独立目录安装，不接管原有 LLMRareCell 项目。解压后执行以下命令。虚拟环境名称可自行更改。

```bash
cd LLMRareCell_V1
python3.9 -m venv .venv
source .venv/bin/activate
python -m pip install '.[test]'
rarecell demo --mode core --out work/demo_core --seed 1009
cat work/demo_core/run/usage.json
cat work/demo_core/evaluation/targets.tsv
```

也可直接从 ZIP 中的 wheel 安装，运行 `python -m pip install dist/llm_rarecell-1.0.0-1-py3-none-any.whl`。该方式同样需要锁定的核心依赖。

core demo 使用真实的稀疏归一化、全基因线性与二阶 sketch、Python igraph Leiden、表达检索和隔离评估。脚本策略先列出候选，再按固定优先级对至多两个初始候选请求全局扩展，并提交通过阈值的提案。命令主动暂停一次，再从 SQLite 恢复。输入含 160 个细胞、96 个基因。Seurat 和 OntoAnno 在该模式明确关闭。每次 demo 使用新的空目录。

```bash
python -m pytest -q
rarecell doctor --config configs/real_provider.yaml
```

doctor 区分本地导入/版本检查和网络调用。上面的 Python 安装足以运行 core demo。真实 A3 还需要 R、固定作者包和已部署模型。R 的私有环境安装命令、精确 commit、R/Seurat 配对和 API 方言见 [THIRD_PARTY.md](THIRD_PARTY.md)。

## 网页任务管理

从安装了当前仓库的 Python 环境启动控制台。之前安装过初版 wheel 的环境可用 `python -m pip install -e .` 更新。首次尝试可以点「运行离线示例」，这会运行真实的 core demo 并生成可下载的 partition 和评估文件。

```bash
# 在仓库根目录运行。安装过初版的环境先执行 .venv/bin/python -m pip install -e .
.venv/bin/rarecell ui --workspace "$PWD/work/tasks" --port 8787
```

在首页「全局模型」填入部署模型 ID（例如 `gpt-6-luna`）并保存，即可同时设置发现和 OntoAnno 的模型。该值保存在当前项目的 `.rarecell-model`，GUI 和同一 checkout 提交的 Slurm 作业均可读取。清空并保存即可恢复 YAML 中的模型。优先级从高到低为任务表单的分项模型或命令行 `--model`、项目全局模型、YAML 中的模型。GUI 任务在入队时保存解析后的配置，恢复任务使用原运行的配置。

打开终端打印的完整私人链接（含 `#` 后的随机令牌）。页面只绑定本机回环地址。每次启动服务会生成新链接，请只在自己的浏览器使用。选择「新建任务 → 完整发现流程」，填写数据集清单或选择已有快照，再选择 A0–A3。真实任务可在表单中填模型名，或指定含模型名的 YAML 配置；相对路径以启动目录为准。A1/A3 可指定已有注释目录，留空则执行 OntoAnno。提供真值 TSV 时会在发现结束并产生结果后单独评估。凭据从启动服务的环境继承。

任务详情提供当前阶段、日志末尾、最近发现事件、细胞/候选计数及已生成的结果文件。可以取消排队、停止运行、用已有 `checkpoint.json` 恢复发现，也可以单独准备、注释或评估。结果未知的 API 调用只有在恢复表单中主动勾选后才会额外重试。停止正在运行的任务会终止其本地进程组；已发出的模型请求可能仍然计费。

每个工作目录有一个串行队列。队列、阶段日志和各任务输出保存在 `--workspace` 下；请选择私有、持久并支持 POSIX 文件锁的目录。关掉浏览器后，执行器继续运行。节点或作业被回收时，重新在相同目录启动服务可恢复排队项；失去执行进程的在途任务标为「已中断」，可手动从发现检查点续跑。恢复沿用原运行目录，结果链接展示该目录当前发布的文件。独立的准备和注释阶段可在详情所示的任务目录中找到产物。

**在 Duke DCC 使用。** 可以在本机按 Duke 的 Remote SSH 指南配置 `dcc-cpu` 别名，再从 Windows PowerShell 打开下列连接。进入该 SSH 会话后，在分配到的计算节点启动上述控制台，并在本机浏览器打开它打印的链接。若使用已经连接到该计算节点的 VS Code，也可在 Ports 面板转发 8787。队列内的任务共享当前 Slurm allocation 的资源和时限。[Duke Remote SSH 配置](https://oit-rc.pages.oit.duke.edu/rcsupportdocs/dcc/vscode/)

```powershell
ssh -L 8787:127.0.0.1:8787 dcc-cpu
```

### DCC 的 R 环境、密钥与批处理

在 DCC 登录或计算节点从 checkout 运行 `rarecell` 时，启动器自动执行 `module load R/4.3.1-rhel8`，并将 `R_LIBS_USER` 指向仓库同级的 `r-libs/MultiGPTRareCell/4.3.1-rhel8`。其他安装位置可用 `RARECELL_R_MODULE` 和 `RARECELL_R_LIBS_USER` 指定；将前者设为空字符串可使用自己已加载的 R 环境。在其他主机默认沿用当前 R 环境。

将已有的 `.api_keys` 放在仓库根目录，每行使用 `NAME=value` 或 `export NAME="value"`。支持引号、空行和注释，值按字面量读取。官方 OpenAI 接口默认使用 `OPENAI_API_KEY`。文件中的凭据覆盖同名凭据环境变量。启动器保留模块设置的 PATH、动态库路径和 R 库路径。DCC 上的 doctor、GUI 工作进程与批处理使用模块解析出的同一个 Rscript；显式填写的 prepare.rscript 或 annotation.rscript 优先。用 `rarecell doctor --config real.yaml` 查看两个组件各自的 `rscript` 字段。密钥由启动进程传给子进程。更新密钥后重新启动 GUI，新提交的批处理作业会在计算节点重新读取文件。可以用 `chmod 600 .api_keys` 限制文件访问。

在仓库根目录执行以下命令。网页与命令行均可设置同一份全局模型。

```bash
.venv/bin/rarecell model --set gpt-6-luna
.venv/bin/rarecell doctor --config real.yaml
sbatch --time=00:10:00 scripts/rarecell.sbatch doctor --config real.yaml
sbatch --time=02:00:00 --cpus-per-task=4 --mem=16G scripts/rarecell.sbatch \
  discover --snapshot work/prepared --annotations work/annotations --run work/A3 --mode A3 --config real.yaml
```

资源、账号和 partition 参数放在脚本名之前。批处理入口使用提交目录中的 `.venv` 和 `.api_keys`，按作业启动时的全局模型运行。要在提交时指定模型，可在子命令后加 `--model gpt-6-luna`。脚本默认使用 Slurm 的 `--export=NIL`，自定义模块和库目录可在自己的作业脚本里设置。GUI 队列在它所在的 allocation 内串行运行。

## 输入

| 文件 | 契约 |
|---|---|
| counts.npz | SciPy CSR，cells × genes。每一项有限且非负。 |
| cells.tsv | 必需 cell_id，唯一且稳定。可含原始 barcode、sample_id、batch 等允许字段。 |
| genes.tsv | 必需 gene_id，可选 symbol、species。用稳定 ID 解决重复 symbol。 |
| dataset.yaml | 四个输入字段组如下。路径相对此文件。 |

```yaml
counts: counts.npz
cells: cells.tsv
genes: genes.tsv
orientation: cells_by_genes
layer: counts
species: human
tissue: peripheral blood
allowed_metadata: [sample_id, batch]
```

示例要求 cells.tsv 确实含有 sample_id 和 batch；没有这些列时用空列表。每一行和每一列分别对齐原始细胞表和基因表。没有额外过滤步骤。counts 在本地按库大小缩放后取 log1p。已有的同类表达层应声明为 log1p_normalized，仍通过 counts 字段指定 NPZ；此时没有可还原的原始 counts，Seurat VST 会记录其输入限制。全零行保留。GT 存在发现目录之外，示例格式为含 cell_id、label 的两列 TSV。评估器要求完整覆盖，缺失的 VDJ 判读需在实验设计中先处理。

examples/minimal_sparse/input 是可读的最小文件实例。examples/datasets 提供 PBMC-TCR、Villani、PGC、Baron、CellSIUS/spheroid 和 GSE156490 入口。MAIT 144/72/36/18/9 的净化输入生成方式见 [DOC_EVALUATION.md](DOC_EVALUATION.md)。

## 真实准备与模型配置

复制配置文件，在全局模型或两个模型字段中填写实际部署的模型名。provider 控制主 Agent，annotation 控制作者的注释模块。密钥通过对应环境变量提供。自定义兼容接口须按服务配置 base_url、endpoint 和 api_style。脚本不假定模型支持 temperature、seed 或固定 token 上限。应在集群作业环境继承凭据和私有 R 库。

```bash
cp configs/real_provider.yaml real.yaml
# 若两个组件使用不同模型，可清空全局模型后分别设置 provider.model 与 annotation.model。
rarecell doctor --config real.yaml
rarecell prepare --dataset data/pbmc_tcr/dataset.yaml --out work/prepared --config real.yaml
rarecell annotate --snapshot work/prepared --out work/annotations --config real.yaml
rarecell discover --snapshot work/prepared --annotations work/annotations \
  --run work/A3 --cache work/query_cache --mode A3 --config real.yaml
rarecell resume --run work/A3
rarecell evaluate --run work/A3 --gt fixtures/pbmc_tcr/gt.tsv --out work/evaluation/A3
```

实验性二阶分支可在 real.yaml 的 prepare.views 中增加 V_TENSOR2，或在 prepare 命令单独使用 configs/tensor2.yaml。不同表示各建一张图。正式比较需要各臂指向同一个 prepare snapshot；发现参数不会反向修改该 snapshot。输入为 counts 时，快照保存 raw.npz 和 data.npz。每种方法的已用参数、依赖、阶段状态和数据哈希记在 manifest.json。线性、二阶与 Gini 的数值定义及边界见 [ALGORITHMS.md](ALGORITHMS.md)。

在依赖已配置的环境测试整个 A3 接口，可运行真实 integration demo。低能力模型可能给出无法配对的注释或无法解析的动作；此时会导出 partial 及其原因。合成输入的基因使用匿名标识，注释可为 unknown。

```bash
rarecell demo --mode integration --out work/demo_live --config real.yaml
```

## 消融与固定策略对照

| 实验臂 | OntoAnno 建议 | 新增数值查询 | 决策策略 |
|---|---|---|---|
| A0 | 关闭 | 禁止 | 真实主 Agent |
| A1 | 开启 | 禁止 | 真实主 Agent |
| A2 | 关闭 | 允许 | 真实主 Agent |
| A3 | 开启 | 允许 | 真实主 Agent |
| fixed | 与 A3 同源 | 与 A3 相同的上限 | 预定优先级 |

```bash
rarecell ablate --snapshot work/prepared --annotations work/annotations \
  --out work/comparison --arms A0 A1 A2 A3 fixed --config real.yaml
rarecell evaluate --run work/comparison/fixed --gt fixtures/pbmc_tcr/gt.tsv \
  --out work/evaluation/fixed
```

每臂创建独立数据库、队列、证据和决定。可预先运行 annotate 一次；省略 ablate 的 --annotations 时会在比较目录中生成共用注释。A0/A2 不加载标签。A0/A1 的 finalizer 仅用预计算成员和确定性回退。每臂对同一合法科学请求首次提出时都扣逻辑预算，真实计算次数和跨臂缓存命中另计。各个已成功视角的两个预定分辨率都有完整的确定性 partition 和单独的评估记录。

固定对照按 outlier 来源优先、其余单例其次、其他候选最后的顺序处理；同层优先小候选，然后用稳定成员 ID、候选 ID 排序。只对初始候选执行一次同规则扩展，按配置分批。提交所有有新增成员的扩展；多个候选争用同一细胞时沿用下述未决规则。此策略和 A3 共用新增查询上限，报告各自实际查询量，可进一步按相同实耗比较。

## 调查、预算和恢复

每次模型调用使用重新组装的完整契约、当前状态摘要、近期工具回执和剩余预算。模型只看到工具名和已登记的实体。科学计算读取本地稀疏矩阵，较长的细胞列表和数值只能分页获取。首次科学请求默认最多 25 个，最多 12 个模型回合，每轮至多 12 个动作，每条谱系默认仅一轮扩展。连续三个无进展回合结束。context_bytes 控制序列化上下文的字节预算。可设 max_reported_tokens 为正整数，收到的累计 token 用量达此值时停止下一次调用；单次回复有可能跨过该值。启用此预算且服务未报告总用量时停止。数值 0 关闭累计 token 检查。精确的扩展阈值和空背景处理见 [ALGORITHMS.md](ALGORITHMS.md)。

```bash
rarecell discover --snapshot work/prepared --mode A2 --config real.yaml \
  --run work/paused --cache work/query_cache --pause-after 2
rarecell resume --run work/paused
```

暂停导出的 frozen.json 及其文件在恢复时保存于 history，原先已评分的部分分配仍可复核。正常提交后重复 resume 只核验并读取冻结结果。中断的本地计算会检查已发布产物并恢复；同一组的成员编辑事务要么全部提交，要么全部回滚。网络断开时远端可能已执行。检查服务记录后，可用 resume --run work/paused --retry-unknown 明确发起一次额外尝试。改变发现模型、配额、方法或冻结输入时使用新的运行目录。若 A1/A3 因注释阻塞，可在原先指定的注释目录完成 annotate，然后 resume。收到格式错误时查 events.jsonl 和 SQLite 中的 Receipt，以确定哪个请求需要澄清。

## 产物与解释

| 路径 | 内容 |
|---|---|
| config_resolved.yaml、project_contract.yaml、manifest.json | 冻结的配置、17 条项目契约和 prepare 清单 |
| state.sqlite、events.jsonl、checkpoint.json | 持久科学状态、有序事件、恢复所需引用 |
| artifacts/ | 注册的数值结果、基线与可验证的内容哈希 |
| candidates.jsonl、candidate_members.tsv | 初始候选与修订、父谱系和真实成员 |
| seeds.tsv、candidate_tasks.tsv | 每个独立 seed、原始 barcode、调查状态、最终去向及全部初始候选的待查清单 |
| annotations.tsv、evidence.jsonl、decisions.jsonl | 注释及其来源、证据卡和依据 |
| partition.tsv、groups.tsv、conflicts.jsonl | 一行一个 cell_id、实际群组的成员哈希、冲突处置 |
| usage.json、run_summary.json、frozen.json | 成本与预算、实际完成状态、冻结文件 SHA-256 |
| 单独的 evaluation 目录 | 匹配的目标 P/R/F1、ARI/NMI、最佳单候选诊断、新成员与碎片统计 |

一个细胞被恰好一个提交候选认领时沿用其标签及状态。零认领或多认领时，标记 unresolved 并保留同一预计算粗群的剩余成员；存在多个粗类配方时依次优先 V_HVG、V_LINEAR、V_GINI、V_TENSOR2，并选最小分辨率。没有基线时，每个此类细胞保留独立未决 ID。partition 中的 claim 组只含无争议成员，groups.tsv 按实际分配重新计算版本，候选的原始成员版本仍可审计。此回退可提高全局 ARI，而目标稀有群 F1 仍可能为零。

评价首先对提交群体和全部 GT 类别做全局一对一最大总 F1 匹配。最佳单候选 F1 是冻结后的另一项离线诊断。不能将二者混用。真实数据部署、重复实验与完整标签要求见 [DOC_EVALUATION.md](DOC_EVALUATION.md)。

## 常见故障

- **Seurat marker 阶段失败。** 使用 Seurat 5.3.0 与 SeuratObject 5.2.0 的配对；固定环境和命令见 THIRD_PARTY.md。查看 manifest 的 stages，修正依赖后生成一个新 snapshot。
- **A3 无法启动。** 查看 annotation_status.json 中的错误、分组数和上游凭据环境变量。明确选择 A2 会新建无标签实验臂；原 A3 保留实际状态。
- **大输入或 Slurm 作业超时。** 按 cells × sketch_dim、HVG 缩放、KNN 和每次全细胞扫描分配内存与时间。scale_max_bytes 约束 Seurat 选定特征层。大样本高维精确近邻可能耗时，尚未发布生产吞吐保证。所有节点能访问同一 snapshot、运行目录和缓存时，可在下一个作业中恢复。底层文件系统需支持 POSIX 文件锁、原子 rename 和 SQLite。
- **退出码 2。** 查看标准错误或 JSON 中的 status/detail。处于 partial 的 discovery 仍可能已有覆盖全部细胞的 partition。评估结果应与 run_summary 一起读取。
- **重复运行输出目录已存在。** prepare 和 ablate 要求新的目标目录。恢复用 resume。重新评分时可复用内容相同的目录；不同冻结版本或 GT 用不同评价目录。

软协议和精确身份边界见 [PROTOCOL.md](PROTOCOL.md)，状态与证据依赖见 [DOC_STATE.md](DOC_STATE.md)，实现和测试索引见 [REQUIREMENTS_TRACEABILITY.md](REQUIREMENTS_TRACEABILITY.md)。

## 复现公开 PBMC 接口测试

SeuratObject 自带的 pbmc_small 可以单独导出原始 RNA counts 和 ID。该示例不导出对象元数据中的分群结果。使用下列命令可以复现输入；具体的真实调用配置与测得状态见 VALIDATION_REPORT.md。

```bash
python examples/export_pbmc_small.py --out work/pbmc_small_input
rarecell prepare --dataset work/pbmc_small_input/dataset.yaml --out work/pbmc_small_snapshot --config real.yaml
```

标记为 partial 的注释批次也可保留经校验的成功行。此时 A1/A3 只使用真实返回并匹配成功的建议，其结果保留 partial。没有一行可用时，将冻结完整未决分配并报告 annotation_blocked。新注释版本若已改变此前看到的证据，请新建实验臂。可以通过 discover 的 --replication-id 或 ablate 的 --replication-prefix 显式绑定预登记重复；未指明的独立运行仅有 run_id。
