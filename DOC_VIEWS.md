# 表示、图和初始证据（reference v1）

输入是非负有限的 CSR。以稳定 cell_id 和 gene_id 排序后处理，原始行列与 manifest 在排序时同步。counts 的行按原总量缩放到 10000 再取 log1p，全零行保留。已归一化输入必须声明 log1p_normalized；直接使用该层。Seurat VST 分支要求实际 counts。

公开基准的非整数线性表达（Villani TPM、Human PGC、GSE156490）可声明为 `linear_expression`。按行总量归一化并取 log1p，Seurat 在 data 层执行 dispersion 筛选，prepare manifest 和发现证据记录选择方法。导出的输入保留原值，此分支不生成 raw.npz。

* V_LINEAR 使用整个 log1p 层做行 L2 归一化，scikit-learn SparseRandomProjection(density="auto")，默认 128 维，并按行 L2 归一化输出。
* V_TENSOR2 对同一全基因单位向量运行 PolynomialCountSketch(degree=2,gamma=1,coef0=0)，默认 256 维，然后按行 L2 归一化。变换分块写入 cells × dimension 的磁盘数组。对非负单位输入，s、s²、αs+(1−α)s²（0≤α≤1）的精确排序相同；有限维近似、归一化、边权和随机性均可影响图。表示间一致仅表明同一数据的分析结果相近。
* V_GINI 是本项目的稀疏选择器。n 个观测中将每个基因的非零值排序，赋以考虑前置零的秩 r；G=Σ(2r−n−1)x/[nΣx]。拟合变量为 log2(max(x)+0.1)，参考 GiniClust3 的两轮最大表达校正。把相同/邻近协变量归至最多 64 个分位箱（每箱目标至少 8 个基因；相同边界合并），对箱中位数插值；排除第一轮残差高于 75% 分位的基因后重新拟合，端点恒定外推。取可变基因中残差最高的前 n_features 个，平票按 gene_id。保留每个基因的 G、检出比例、拟合值和残差。本实现以稀疏闭式 Gini、分箱趋势、效应排序替代作者的 dense + ε、LOWESS 和残差正态尾概率。没有运行 GiniClust3 consensus。用稀疏隐式中心化、每个所选基因除以样本标准差的 ARPACK SVD 得到 PCA 得分，保留载荷和均值/尺度。
* V_HVG 在 R 中调用 Seurat 5.3.0 的 FindVariableFeatures、ScaleData 和 RunPCA。counts 使用 vst，linear_expression 使用 dispersion。选定特征的稠密缩放层有显式内存上限。行数或变异不足时记录阻塞。HVG 与 partition marker 分文件。
* 生产后端为 Seurat。所有启用表示在各自的行单位化嵌入上独立运行 FindNeighbors(RANN, euclidean, prune.SNN=1/15)，FindClusters(algorithm=4,group.singletons=FALSE)。固定 future sequential 和随机种子。作者源码在并行路径调用 GroupSingletons 时未传 group.singletons；串行路径显式传入。串行模式仍把不同 singleton 统一写成 singleton 字面标签，桥接层按 barcode 拆为独立成员。保留 R 原标签和实际集合 ID。
* 显式 Python 核心后端在相同表示上做精确单位向量 Euclidean kNN（self 排除）。取有向 kNN 边的并集，以两个邻居集合（各自包含自己）的 Jaccard 加权，保留权重 >1/15，用 igraph 的 C 实现运行 Leiden 模块度优化。两个后端各自有可复用配置；Seurat 的 SNN 边集合可与 Python 后端不同。缺少生产依赖会被记录。
* 零表达、零嵌入者不会参与该视角的图或 outlier 距离拟合，在 partition 中各有独立 singleton。其他 outlier 规则为第 k 邻居距离大于中位数 + 3.5 × 1.4826 × MAD + 1e−8。固定余量抑制完全重合的数值噪声；不会强取固定比例。每个 singleton/outlier 对应一个独立 seed，跨视角的同一 barcode 可有多个来源。互为第一邻居且 cosine ≥0.90 的非孤立对是局部候选，默认取得分最高的 64 对；平票按 barcode。达到该阈值表示可提交探索。
* 完全相同集合按排序后的全体成员内容去重，并联合来源；保留近似重叠候选。普通群 marker 在规范化层上用组均值减全集剩余成员均值，保留组内外检出比例。singleton 优先同视角局部邻居对照，无邻居时使用全部其他细胞；背景为空或查询集合表达全零时无 marker。初始卡展示绝对效应最大的 12 个基因（等值按 gene_id）。Seurat HVG partition 有至少 3 个成员和至少 3 个有效背景细胞时，另行在 RNA assay 的可用基因上运行 FindMarkers(wilcox)。用于注释的正 marker 独立保存，HVG 选择仅提供另一份特征证据。标签不参与特征选择、建图、outlier 或成员整理。

## 输入、冻结与运行记录

输入允许 cells.tsv 的 cell_id、可选原始 barcode 和显式列出的 sample_id、batch、donor_id、library_id、chemistry、timepoint、assay。未获准的细胞列不复制。存在原始 barcode 时，排序键采用 (barcode, cell_id)，允许不同样本重复使用 10x barcode。基因表保留 gene_id 和可选 symbol、species。供发现端读取的 manifest 只有净化的数据、物种、组织及分析记录。已有输入目录不能被 prepare 原地改写。

snapshot 以整目录原子发布，保存原始 CSR（输入为 counts 时）、规范化 CSR、每个表示的 SNN、图行标识、所有 partition、候选与种子。稠密文件只涉及低维表示或选定的特征。Seurat 可调整 scale_max_bytes，默认单个选定特征缩放层上限为 512 MiB，进程还有 R 对象和临时矩阵开销。prepare 配置中 rscript 和 r_timeout 控制外部 R 进程。R marker 在非零输入细胞内比较指定集合和剩余者。Seurat 5.3.0 应搭配 SeuratObject 5.2.0；5.4.0 删除了 Wilcoxon 路径使用的 PackageCheck。默认用户参数无需在已有 R 全局库中修改。

manifest.json 中 initial_evidence_hash 是全部已冻结科学文件的内容索引摘要，snapshot_hash 还包含代码、依赖和该次执行耗时。通过 verify_snapshot 校验清单和各文件的 SHA-256。不同实验组使用同一 snapshot 时这两个哈希均相同。生产 Seurat、Python igraph 的成本和结果分别记录；对大样本精确 kNN 的运行时间依赖细胞数和表示维度。

## 主 Agent 对预计算结果的读取

每个实验组都登记相同的两种 Seurat 初始观察。全局 `global_feature_selection` 给出该输入对应的 VST 或 dispersion 特征列表的 evidence_id；候选卡的 `seurat_wilcoxon` 列出该细胞集合在相应分辨率的实际正向 Wilcoxon marker 表及 evidence_id，和描述性均值差字段 `markers` 分开。小于 3 个成员或不足 3 个参考细胞时没有群体 DE 表。各组都能在 `get_evidence` 中传入单个 evidence_id、offset、limit 分页读取已冻结的基因、效应、检出率和 P 值；只读操作不产生科学查询。背景是所有非零输入细胞中该群之外的成员。每张表保留来源文件哈希、集合哈希和上游选择依赖。跨分辨率的完全相同对照及数值结果共用一项观察。OntoAnno 解释按 partition_id 和其输入 marker_hash 绑定到具体 DE 观察。

每个初始候选在每张可用 SNN 独立取得结构摘要，保存在 initial.jsonl 的 structure、候选证据和各实验组的候选卡中。以该视角 graph_rows 中实际存在的候选成员集合 S 计算：分子为从 S 出发、也终止于 S 的权重和；分母为从 S 出发的全部权重和。两项均排除自环。无向图内部边在分子和分母各算两次，向外的边在分母算一次。分母为零时 within_fraction 为 null。graph_coverage 是候选成员在该图有节点的比例，n_zero_degree 只计已表示且排除自环后度为零的节点。图未生成时使用 graph_unavailable。同一输入的多个视角可供比较，这些数值本身不增加独立验证次数。计算只发生在 prepare。
