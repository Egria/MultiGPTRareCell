# V1 验证报告 · 2026-11-21

**结论。** 稀疏输入、多视角准备、真实局部数值检索、隔离消融、恢复、全细胞分配和独立评估已执行。真实 Seurat 和固定 OntoAnno 能对接真实 CPU 模型。在本机小模型测试中，接口可用和自主稀有成员发现的表现有明显区别。

## 环境与安装

| 环境／检查 | 实测结果 |
|---|---|
| Linux CPU，Python 3.9.25，启用 R/Seurat/OntoAnno | **109 passed**。11 条 sklearn 提示均来自小型 fixture 的线性投影维数大于基因数。 |
| Linux CPU，Python 3.12.14，仅 Python 核心环境 | **106 passed，3 skipped**。两个真实 Seurat 测试及固定 OntoAnno 接口探测需要外部 R。11 条同类投影提示。 |
| 核心库 | numpy 1.26.4、scipy 1.13.1、scikit-learn 1.5.2、igraph 0.11.8、PyYAML 6.0.2，见 pyproject.toml。 |
| 真实 R 库 | R 4.5.3、Seurat 5.3.0、SeuratObject 5.2.0、Matrix 1.7-6、leidenbase 0.1.36、RANN 2.6.2、OntoAnno 2.0.0、ellmer 0.5.2。 |
| 安装验证 | 将交付源码解压至新目录，使用独立 Python 3.9.25 的 pip 安装带 test extra 的包并执行 CLI；另外构建 py3-none-any wheel，装入新的 Python 3.9 环境。依赖一致性检查通过，已安装的 wheel 在项目目录之外成功执行 core demo、暂停、恢复和评价。 |

实际测试使用项目内 pytest 套件、私有 Rscript 和 R_LIBS_USER。可用下面的命令复现。R 不在 PATH 时，RARECELL_RSCRIPT 可指向测试专用的已安装解释器；正式运行使用 prepare.rscript、annotation.rscript 或 PATH。wheel 与源码一同交付。

```bash
python -m pip install '.[test]'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q
rarecell demo --mode core --out work/check1009 --seed 1009
rarecell demo --mode core --out work/check1013 --seed 1013
rarecell demo --mode core --out work/check1019 --seed 1019
rarecell resume --run work/check1009/run
```

## 数值结果与对照

- 单个 seed 从真实稀疏计数、规范化表达与 SRP 中检回全部 **4/4** 个同信号细胞。两个普通粗群均包含目标成员，检索覆盖了未进入局部背景的同信号细胞。不同大小的植入群输出 **1、2、5、8** 个成员，孤立 seed 保持单例。
- 对同组成、不同库大小的细胞验证了归一化。全零细胞可保留，无有效背景或签名时有具体状态。固定单例在五个共享 Gamma/Poisson 表达分布的噪声对照中均未新增成员。
- 真实二阶 PolynomialCountSketch 与精确核在小矩阵上做了比较。64 维时，三个预设种子的内积 RMSE 为 0.21236、0.18241、0.12498；2048 维为 0.01970、0.04570、0.13483。维数增加改善了平均误差，第三个种子的误差略升。测试另核验非负单位输入下 s 与 s² 及其非负加权混合的精确排序一致。实际图可因有限维噪声与边权变化而变，已分别导出 baselines。没有证据证明该分支提升临床数据集的稀有群 F1。
- 160 × 96、4,708 个非零项的合成矩阵在真实四分支 Seurat 配置下生成了两个分辨率共 1,280 条基线分配及 312 条正向 Wilcoxon marker 记录。一轮代表性准备实测 26.4 秒墙钟时间。该次使用 k=7，HVG/Gini 各 8 个 PC，线性 32 维，二阶 64 维。各处理阶段状态、输入和参数均进入 snapshot。
- 实际构建的 SNN 图用于按候选计算内部权重占比、图行覆盖与零度数。测试以图行重排、缺失个体和零分母检查数值，并核验 A0/A1 可以分页读出真实 HVG 和基于固定分区的 DE 统计。OntoAnno 的解释节点关联对应 marker 观察。

### 三个预登记的 core demo

| 生成/投影 seed | 稀有群精确率 | 召回率 | 提交 F1 | 未决细胞 | 真实数值查询 |
|---|---:|---:|---:|---:|---:|
| 1009 | 1.000 | 3/8 | 0.545 | 157/160 | 2 |
| 1013 | 0 | 0/8 | 0 | 160/160 | 2 |
| 1019 | 0 | 0/8 | 0 | 160/160 | 2 |

每次保留 160 个细胞。报告按提交 partition 的全局一对一匹配计算，包含未检出重复。三次目标 F1 均值 0.182、中位数 0。运行使用显式脚本策略；表达对比、成员提案、冲突处理、恢复和评分走实际代码。生成器只返回发现输入地址，GT 位于独立目录。完整分布及各基线口径见 DOC_EVALUATION.md。

## 真实集成

从已安装的 SeuratObject 公开示例 pbmc_small 导出 80 × 230 的 RNA counts、4,456 个非零项，以及基因/细胞 ID。未导出原对象的聚类或注释列。examples/export_pbmc_small.py 可复现此输入。四个表示均成功构图，Seurat marker 阶段通过。用于比较的预处理 snapshot 为 `e89348d1d3103eee5e548a0ab52f72a94b3f358627da12db4032a518808a974c`。这是软件接口的真实表达样例，规模和特征 panel 均受原包限制。

真实作者注释使用 OntoAnno commit `ceda74c7d74af67a02093b07dd87623d846a59fe` 及 Qwen2.5-3B-Instruct Q4_K_M（本机 llama.cpp CPU 服务）。单群单批、显式换行停止条件的运行获得 **5/5** 份对齐成功的原始响应。冻结结果中有四个混合建议和一个 Unknown，ontology_id 均未映射。未提供可信 GT，未计算 PBMC 准确率或真实稀有群 F1。自由格式的早期实跑只得到 2/5 个通过校验的批次，部分状态和有效返回的区分经过复核。关于作者源码、原始响应和缓存验证见 THIRD_PARTY.md。

```bash
python examples/export_pbmc_small.py --out work/public_example
rarecell prepare --dataset work/public_example/dataset.yaml --out work/pbmc_seed1009 --config configs/pbmc_small_smoke.yaml
rarecell annotate --snapshot work/pbmc_seed1009 --out work/pbmc_labels --config configs/pbmc_small_smoke.yaml
rarecell ablate --snapshot work/pbmc_seed1009 --annotations work/pbmc_labels --out work/pbmc_arms --config configs/pbmc_small_smoke.yaml
```

该配置要求启动模型服务并提供对应模型别名，ellmer 的本地 vLLM 兼容客户端还读取 VLLM_API_KEY。代码包不含模型权重。独立的匿名合成 integration demo 在真实 Seurat 准备成功后，Qwen2.5-0.5B 两次把输入 marker 列表回显为标签。适配器拒绝两项，标记 annotation_blocked，160/160 个细胞都有未决回退分配。这个阴性实验保留在验证结论中。


### 主 Agent 和固定策略的实跑状态

开发过程用相同公开 snapshot 运行一次四臂和固定策略。主模型为本地 Qwen2.5-3B；工程预算为 4 回合、3 个科学查询、连续 2 个无进展回合停止。这次使用开发期含单个请求示例的提示。所有臂的状态如下。

| 实验臂 | 终止原因 | 模型回合／策略回合 | 新科学查询 | 实际数值执行 | 模型是否提交 | 未决细胞 |
|---|---|---:|---:|---:|---|---:|
| A0 | stalled | 3 | 0 | 0 | 否 | 80/80 |
| A1 | stalled | 3 | 0 | 0 | 否 | 80/80 |
| A2 | stalled | 4 | 0 | 0 | 否 | 80/80 |
| A3 | submitted | 4 | 0 | 0 | 是 | 80/80 |
| fixed | submitted_after_budget_exhaustion | 2 | 3 | 3 | 是 | 80/80 |

逐项核对的五个运行具有相同的初始成员集合摘要，A0/A1/A2/A3/fixed 的注释记录数分别为 0/5/0/5/5，每个最终表均有 80 个唯一细胞 ID。部分运行中模型反复模仿提示里的 list_candidates 示例。发货版移除固定动作示例，去重重复证据页，保留完整契约和最后一轮提交提示。用同一个已验证注释库做了附加 A3 smoke，并保留了协议修正前后的运行。

| 发现模型 | 终止原因 | 回合 | 新科学查询 | 已提交 | 总 token |
|---|---|---:|---:|---|---:|
| Qwen2.5-0.5B | stalled | 2 | 0 | 否 | 6349 |
| Qwen2.5-3B，嵌套动作修正前 | max_agent_turns | 3 | 0 | 否 | 8021 |
| Qwen2.5-3B，嵌套动作修正后 | max_agent_turns | 3 | 0 | 否 | 8289 |

0.5B 两次到达单次 256 个输出 token 的上限；均未执行截断动作。3B 在第一次新提示测试中主动请求了一个 2 细胞候选的扩展，使用单层 action.type 包裹动作。旧解析器将其阻塞。加入精确的容器解析后，单独回放该份真实回复产生了一次真实数值查询：24 个背景细胞、4 个表示、0 个达标新成员；这只属于回放验证。修正后又通过真实 HTTP 推理从新运行目录独立复跑，表中最后一行记录其结果。开发期和发货版结果分别保留，没有按是否成功挑选运行。所有公开 PBMC 数值只用于接口和决策行为检查。完整机器可读摘要在 validation/observed_runs.json。

## 持久化、权限及故障路径

| 检查 | 实际执行的判断 |
|---|---|
| 中断恢复 | 终止真实子进程制造“缓存产物已写，SQLite 尚未提交”，恢复取得同一数值结果。另一测试在工具回执后暂停，再恢复至最终提交；计算执行次数不增加。部分导出可按其原始哈希独立重评。 |
| 并发与重试 | 两个进程请求同一数值查询，只有一个实际执行。第二个实验臂同样扣一次逻辑预算。故障注入验证暂时失败的有限重试、永久失败局部阻断和损坏产物校验。 |
| 证据依赖 | 拒绝悬空父节点、环、同版本内容替换和发现端升级独立验证。相同证据/注释注入 1、3、10 次后，证据和传入的有效内容相同。按 marker 选中成员再重查该 marker 时保留选择依赖。 |
| 软协议与批处理 | 测试中英文别名、围栏、尾逗号、数值字符串、单元素列表、真实 ID 歧义、scope 冲突、AND/OR、截断和混合批次。单项失败后同批独立动作仍可成功。事务冲突回滚后状态可继续使用。 |
| 实验臂隔离 | A0/A2 读取不到注释哨兵；预先加热缓存不会允许 A0/A1 发起数值计算。每臂独立的候选修订和任务状态由集成测试验证。 |
| HTTP 测试性质 | 本地协议 fixture 仅验证 chat/responses 载荷、完整性、异常和计数；真实模型测试使用运行中的 llama.cpp。未知远端执行状态在显式授权重试前不会重复发送，最多一次额外尝试。 |

## 环境阻塞与尚未运行

1. PBMC-TCR/MAIT、Villani、PGC、Baron、CellSIUS/spheroid、GSE156490 的项目用真实矩阵、分析清单和冻结 GT 未随任务提供。这六个模板没有实测 F1/ARI/NMI；公开 pbmc_small 和合成输入不代替其中任何队列。
2. 未配置商业 API 凭据。Chat Completions 的真实网络调用在本地兼容服务验证；Responses 方言使用 HTTP 协议测试。OpenAI、Anthropic、Gemini 等远端部署、实际账单及模型差异尚未实测。
3. 已执行 CPU Linux 本地测试，未向 Duke DCC 的 Slurm 队列提交作业；大队列的时间、峰值内存、跨节点共享文件系统行为未压测。当前结果不能证明 LLM 相对固定搜索有增益。
