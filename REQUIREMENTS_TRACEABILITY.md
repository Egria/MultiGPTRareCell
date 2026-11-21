# 需求到实现及验证的索引

契约原文的执行版为 PROJECT_CONTRACT.yaml，安装包在 rarecell/resources 中携带相同内容。每次启动均验证完整的 C01–C17、版本及字段，再写入运行快照。以下路径相对代码根目录。测试可用 `python -m pytest tests/文件名.py -k 关键词` 单独运行。

| ID | 实现与强制边界 | 对应测试或可观察记录 |
|---|---|---|
| C01 成员优先 | discovery.make_revision 与 finalize.assignments 接受 unknown、provisional，成员哈希独立于标签。 | test_agent 的 total_conflicts_and_seed_can_be_reassigned；test_lineage 的 annotation_consistency |
| C02 全细胞保留 | data.read_dataset、prepare.record_partitions 与 finalize.assignments 分别验证输入、保留无表达行并检查完整互斥。 | test_prepare 的 depth_control、singleton_label_recovery；test_agent 的 resume_after_real_queries |
| C03 单例身份 | prepare.Pool.seed_rows 保存逐细胞 seed_id、barcode、候选版本；seeds.tsv 另报最终归属。 | test_prepare 的 singleton_label_recovery_and_traceable_seeds |
| C04 开放搜索 | Science 对所有可测细胞打分。团内外、粗类、预注释不限制检索和成员修订。 | test_science 的 singleton_recovers_across_coarse_classes；test_scientific_permissions 的 annotation_visibility |
| C05 最终互斥 | 批量修订使用事务，finalize 对所有提交候选建立逐细胞认领表；冲突回退并写理由。 | test_agent 的 total_conflicts；test_state 的 immutable_records_and_atomic_revisions |
| C06 唯一注释 | annotation 与 r/annotation.R 只调用固定作者包的注释/ontology 子函数。 | test_integrations 的 real_ontoanno_interface_if_installed；已执行的 PBMC 注释在验证报告中 |
| C07 可追溯 | prepare snapshot、FrozenSeurat、State 与 Controller 记录输入、真实集合、方法、父证据、结构与 marker 的内容引用。 | test_base_evidence；test_lineage 的 scientific_evidence_pages；test_state 的 scope_method_and_replication_identity |
| C08 防重复强化 | 同键查询复用数值，Evidence 不可变，QueryUse 保存新用途；注释重述关联父观察，未定义独立投票数。 | test_scientific_permissions 的 repeated_annotation_input 与 equivalent_requests；test_state 的 evidence_deduplication_and_dag |
| C09 探索/验证 | post-selection 查询保留真实选中基因、背景、成员、参数和 selection_history_ref，禁止发现端发放独立验证。 | test_scientific_permissions 的 marker_recheck；test_state 的 selection_dependent_claim_cannot_promote |
| C10 GT 隔离 | data 白名单、注册工具、只读准备快照。evaluate 的 GT 在冻结校验后单独加载。 | test_evaluation 的 freeze_before_reading_gt 和 frozen_escape；test_scientific_permissions 的 unlicensed_queries |
| C11 持久运行 | SQLite、有序事件、不可变产物和已派发完整消息；resume 重建上下文与部分结果。 | test_state 的 process_crash_between_artifact_write_and_db_commit、parallel_runs_single_physical_execution；test_agent 的 interrupted_export |
| C12 软协议 | protocol.parse_message 确定性解析，受限字段和动作别名，歧义局部阻塞；精确解析控制器实体。 | test_protocol 全部；test_scientific_permissions 的 equivalent_requests_share_numeric_result |
| C13 有限往复 | 有界 turn、逻辑查询和通信，重复科学请求不增加新进展；停滞与耗尽导出未决任务。 | test_agent 的 budget_or_stall_terminates 和 warm_cache；test_state 的 run_budget_applies_before_cache |
| C14 直接沟通 | 每轮 GOAL 与契约要求简短依据和具体反证，局部失败返回原因，解释字段保留在消息/Receipt。 | test_protocol 的 ambiguity_and_conditions；test_agent 的 mixed_batch。模型的实际语言质量需另看真实回执。 |
| C15 允许零扩展 | 标签未知、零向量、无对照和低分不会触发强制取满；收紧上限有显式 incomplete。 | test_science 的 distinct_population_sizes、iid_count_noise、zero_vector_small_background |
| C16 排序与身份 | 按 barcode/cell_id 对齐，统一截断顺序，数值排序含稳定 tie-break。 | test_prepare 的 row_reordering_graph_candidates_markers；test_science 的 original_barcode_controls_ties |
| C17 真实状态 | prepare/annotate/run 分阶段记录实际实现及 blocked/partial，外部未知响应保留审计，未测费用留空。 | test_prepare 的 missing_r；test_integrations 的 unknown_dispatch_and_receipt_recovery；test_agent 的 ambiguous_api |

## 主要阶段

| 任务书范围 | 实际入口与交付 | 检查 |
|---|---|---|
| §3–4 安装和输入 | pyproject.toml、rarecell.data、examples/minimal_sparse、六个开发入口模板 | test_prepare 验证输入类型、维度、数值与归一化；独立 wheel 安装和 CLI demo 见验证报告 |
| §5–8 真实表示与注释 | r/seurat.R、views、prepare、base_evidence、r/annotation.R、annotation；作者版本和接口见 THIRD_PARTY.md | test_prepare 的四视角真实 Seurat、Gini 公式、TensorSketch 内积；test_integrations 的真实作者 API 探测和缓存；test_base_evidence 检查实际 DE 和 HVG 的独立来源 |
| §9 一轮检索 | Science.execute、Controller.scientific、lineage_roots、真实候选修订 | test_science 的跨粗类、孤立、丰度、技术缩放及噪声对照；test_lineage 及 test_scientific_permissions 检查通过别名、继承和拆分后的一轮深度限制 |
| §10–12 Agent/协议/缓存 | Controller.context、run.process_message、providers.complete、State.query；受约束分页，单主 Agent | test_protocol、test_state、test_lineage、test_scientific_permissions、test_agent 的中断恢复与预算 |
| §13 partition | finalize.assignments 与 export；群组版本、冲突表、独立种子状态一并冻结 | test_agent 检查覆盖、冲突、释放种子、恢复时哈希和已提交结果 |
| §14 消融/评估 | cli.ablate、固定策略、全部固定聚类基线、evaluation.evaluate/aggregate、benchmark_views | test_evaluation 稀疏一对一匹配与密集参照、缺失/重复 run、未检出；跨臂缓存计费和标签可见性见 test_scientific_permissions |
| §15–18 CLI 与交付 | 8 个 CLI 子命令、独立评估工具、示例 YAML、中文指南、验证报告、ZIP | VALIDATION_REPORT.md 按实际运行区分科学计算、线协议 fixture、脚本策略、真实模型和未测试数据集 |
