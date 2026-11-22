# 冻结后的成员评价

## 从分配表计算指标

发现流程提交后，评估器先核对 `frozen.json` 中各文件的 SHA-256，再打开单独的 GT TSV。冻结清单至少覆盖 partition、候选索引与成员、各分支基线、用量和 run summary。GT 至少含 `cell_id` 和 `label` 两列，ID 应唯一，且与提交的 partition 完全一致；每行记录一个适用于本评估的真值类别。缺测 V(D)J 等信息时，需要在评估设计中区分未知与阴性。该接口面向单标签全细胞基准。完整分配中 `unresolved` 单元仍然参与评分，其最终 group_id 也参与竞争。

每个 GT 类别都是一个 target，每个提交 group 是一个预测集合。程序按 TP 交集大小计算两者的 F1，求所有类别的一对一最优匹配，使各类别 F1 之和最大。实现用稀疏二分图加每类别一个 0 分虚拟节点。零交集不能构成真匹配。矩阵按 GT 名称和 group_id 稳定排序，在相同软件版本下，等价解的选择与输入文件行序无关。缺少可分配群体时，该目标的 precision、recall、F1 记 0，group_id 留空。未出现在 GT 的类别无可定义的 recall，不生成该类一行。该匹配作用于所有类别，包括作为背景记录的类别。命名和 ontology 不参与匹配。

- `targets.tsv` 逐类别给出 TP/FP/FN、precision、recall、F1、目标分散到的组数 `fragments`、匹配组包含的 GT 类别数及未匹配状态。
- `metrics.json` 保存整体 ARI、arithmetic NMI、macro-F1、混杂组和碎片化类别数、额外组 ID、按状态计数、unresolved 比例、候选数量及初始候选数量。ARI/NMI 使用完整原始标签向量，完全独立于匹配过程。baselines 段按 `view,resolution` 逐个评分；每条基线都必须覆盖全体细胞。usage 和 run_summary 原样保留，便于结合查询预算、执行/缓存、provider 消耗、待完成状态及恢复记录一起阅读。
- `candidate_diagnostics.tsv` 对每个目标从冻结候选中选**单个** F1 最高的候选，并分别给出仅限初始库的同一指标。并列时选择最小 candidate_id。该诊断使用了 GT。候选可以重叠，一个候选可以成为多个目标的诊断匹配。主结果使用提交 partition 的一对一匹配。
- `new_members.tsv` 比较目标在最终匹配 group 里的 TP 与初始库最佳单候选的 TP，记录前者多出的真成员。候选诊断另报丢失的 TP 数。这是回顾性的成员差集，归因需要再查看活动查询和决定。
- `revision_changes.tsv` 逐个非初始候选，与声明的直接父版本成员并集比较，报告新增或去除的 barcode、GT 及该细胞的最终判定是否引用此候选。无父候选的新增成员会依定义记录为全部新增。并行候选可重复涉及同个 barcode，因此按 revision 分别解释。
- `freeze_sha256` 识别完整冻结清单，`partition_sha256` 识别实际提交分配，`snapshot_hash` 识别共同预处理，`gt_sha256` 识别这次评分的 GT 文件。

独立评估目录放在 discovery run 外。可对相同冻结结果及 GT 重新输出。不同版本使用不同评估目录。有 partial 状态的运行仍可评分；来源 run 状态保存在指标文件内。

## 重复运行及缺失结果

开跑前复制 `examples/replications.yaml`，预先填写全部预定重复的唯一 replication_id、实验臂和预期 metrics 文件，按数据集、丰度视图、策略区分 arm。各文件路径相对计划文件；计划描述重复设计，实验程序实际使用的 seed、数据和配置由各运行记录。每个 arm 只代表一个明确研究条件。

```bash
python -m rarecell.evaluation aggregate --plan examples/replications.yaml --out work/evaluation/replications.json
```

报告保存计划哈希、所有已规划条目及其完成状态、每次重复的原始分数、每个 arm 的均值/中位数/极值/样本标准差。每项指标同时记录实际贡献的重复数；单次重复的样本标准差留空。没有结果的重复列为 missing。已评分 partial run 单独标明上游状态并保留在分布中，确保读者看到实际运行结果。同一个冻结运行和同一 GT 不能作为两个独立重复供统计。

## 开发数据模板与合成输入

`examples/datasets/{pbmc_tcr,villani,pgc,baron,cellsius_spheroid,gse156490}/dataset.yaml` 为入口模板。把正确的 CSR 与 cell/gene 表放在对应目录，检查组织/物种、矩阵层及方向。需要 GSE156490 等具体样本信息时请填写已核实的样本元数据。GT、目标清单和 barcode 的判读证据另存为 evaluation fixture。模板不包含真实患者数据或作者标签。

`examples/minimal_sparse/input/` 是可读的 CSR、cell/gene TSV 和 YAML。相邻的 `evaluation_fixture` 只供离线评分。合成样例有 160 个细胞、96 个基因；在两个常见背景各植入 4 个具有共同表达信号的细胞。其数值结果只代表这项合成条件。生成器只返回 discovery manifest 路径。执行 `generate(root, seed)` 可在新的空输出位置重建独立稀疏输入及单独 GT。

## PBMC-TCR 的 MAIT 丰度压力系列

此命令属于 benchmark 数据构建。它是独立读取 GT 的程序，输出可交给 prepare 的净化输入及单独评估材料。

```bash
python -m rarecell.benchmark_views --dataset data/pbmc_tcr/dataset.yaml \
  --gt fixtures/pbmc_tcr/gt.tsv --out work/mait_series --target MAIT \
  --sizes 144 72 36 18 9 --seed 1009
```

按照种子与稳定 ID 的 SHA-256 排序生成嵌套子集，各视图保留相同的非目标细胞。原矩阵最小目标数须满足请求的最大丰度；目标大小不足会明确失败。净化结果进入 `discovery/view0…view4`，只保留通用元数据和完整 cell_id 表。稀疏矩阵、基因顺序和表达层随对应 manifest 保持一致。标明哪些 barcode 属于目标的表、GT、抽样种子和选择记录只在 `evaluation_fixture/view*/` 下。每个发现运行单独读取对应 view 的 manifest，禁止把多个 view 送入同一 Agent 进行成员反推。嵌套丰度共享来源，不能作为五个独立外部验证队列。

## 已运行的三次离线重复

先登记了数据生成种子 1009、1013、1019。登记时分布报告为 3 次计划、0 次已评分。随后在 Python 3.9.25、项目锁定的 CPU 依赖、单线程 BLAS 下分别执行 `python -m rarecell demo --mode core --out <新目录> --seed <种子>` 并调用 aggregate。每次各含 160 个细胞、96 个基因、8 个植入的稀有细胞。两种启用的表示为 128 维线性投影与 256 维 TensorSketch。该小输入的线性投影增加维数，scikit-learn 据实提示。常见群体有两种，稀有细胞均分到这两个背景。每次都完成两次逻辑科学查询，暂停并恢复一次，保留 160/160 个最终记录。

| 种子 | 提交稀有群 P / R / F1 | 初始库最佳单候选 F1 | 最终库最佳单候选 F1 | 未决细胞 | 全体 ARI | 全体 NMI |
|---|---|---:|---:|---:|---:|---:|
| 1009 | 1.000 / 0.375 / 0.545 | 0.222 | 0.545 | 157/160 | 0.937 | 0.882 |
| 1013 | 0 / 0 / 0 | 0.400 | 0.400 | 160/160 | 0.902 | 0.850 |
| 1019 | 0 / 0 / 0 | 0.091 | 0.091 | 160/160 | 0.902 | 0.850 |

三次均进入分布。稀有群 F1 的均值为 0.182、中位数为 0、样本标准差为 0.315。整体 ARI 均值为 0.913，NMI 均值为 0.860。种子 1009 时，提交集合从初始单例扩展到 3 个真成员。其余两次提交的全体分配沿用两个常见群的未决回退，稀有群未得到一对一匹配。全部查询由显式脚本策略提出；这些结果描述合成表达、数值工具和冻结评价的行为。真实 LLM 自主发现的表现仍需独立运行。

独立目标参考可用 `cell_id, target_id, reference_kind, state` 四列 TSV，每个目标与参考来源必须覆盖冻结分区中的每个细胞一次。对每个目标单独在 positive/negative 细胞上选 F1 最大的最终群和各基线群，unknown/disputed 排除，零阳性时 F1 留空。目标可重叠，每个结果保留其分母。metrics.json 以 `reference_profile=target_membership` 标识该评估，ARI/NMI 留空。
