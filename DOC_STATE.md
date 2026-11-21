# 持久状态与恢复

每次发现运行拥有自己的 `state.sqlite`、`events.jsonl` 和 `artifacts/`。不同消融可以将本地科学工具的确定性结果写入同一计算缓存。调用前由 controller 判定模式、作用域和种子 lineage 权限。缓存只按已获准的请求返回数值，候选、注释、标签、任务、决策和对话均在 run 内管理。

## 状态对象和标识

| 表或标识 | 作用 |
|---|---|
| `record(kind,id,body)` | 带版本的对象，同一 `kind,id` 上再次写入必须字节等价。候选成员修订使用新 ID。 |
| `meta` | 可更新的恢复游标、当前任务等；与多对象编辑放在同一事务提交。 |
| `query.key` | `canonical_query_key`。对版本域和规范化请求的 UTF-8 紧凑 JSON 做 SHA-256。键顺序不敏感，数组顺序保留。controller 负责在此之前解析集合、基因、方向、背景、参数、方法和数据版本。 |
| `attempt.id` | `execution_attempt_id`。每次准备实际调用本地计算分配新的随机唯一 ID。表保存开始时间、状态、CPU/耗时和失败原因。缓存命中不会创建计算尝试。 |
| `Evidence.evidence_id` | 证据内容、祖先、选择历史和版本的标识。默认用内容哈希；从可信上游导入时也可提供指定 ID，绑定后不可换内容。 |
| `message_id` | controller 对通信或协议修复的标识。多个消息可映射到同一个逻辑查询。 |
| `replication_id` | 上游用于标识一次预定种子/样本/扰动重复的实验元数据。影响数值的改变也须写入规范查询。技术重试不创建生物学重复。 |
| `event.seq` | run 内已提交事件的顺序号。JSONL 从此表增量追加，恢复时与 SQLite 对齐。 |
| `artifacts/<sha256>.json` | 本地工具输出按内容寻址，读取时验证哈希，只支持已登记的精确引用。 |

同一证据登记 1、3、10 次都只保留一条。父节点存在性和不可变 ID 共同保证新边只会指向已存在的祖先，自指和用改写制造回路均拒绝。视角、统计量和解释保留各自的数据来源、方法、父证据和选择历史。V1 的 discovery 端没有注册独立验证协议，拒绝 `kind=validation` 或 `purpose=independent_validation`；GT 评估从冻结输出另行执行。

## 本地查询的提交顺序

1. 用 SQLite 写事务在本 run 扣一次逻辑预算。相同键重读不会再次扣费；另一个实验组必须自行提出请求并占用其预算。
2. 获取 run 和共享缓存的每键文件锁，先查本地成功记录，再校验共享产物的请求键和结果哈希。
3. 实际计算前持久登记 `running` attempt。成功后先将带结果哈希、时间及 attempt ID 的 JSON 写入临时文件、`fsync`，然后以硬链接原子发布缓存。
4. 将结果注册为 run 内内容寻址的 artifact，在同一写事务内提交 attempt、查询结果引用和事件。嵌套控制操作使用 SQLite savepoint；本地计算期间不持有控制器写事务。
5. 恢复时检查待定请求，若有已发布结果则完成登记并复用。若进程在计算中断且没有产物，同键可有限重试。打开另一个活跃进程的 run 时跳过尚被锁住的请求。损坏或内容冲突的成功产物报告完整性错误。

初始默认每个查询最多 2 次本地实际尝试。超时和 OS 错误允许在限额内重试。非有限数值、类型错误等永久问题局部 blocked。用尽临时重试额度后状态保留为 retryable，恢复后重复消息不会再自行增加次数。条件修复后，受信任的调用方可以显式提高 `max_attempts`，或使用随后生成的匹配缓存。不存在成功证据的尝试仍计成本。远程模型调用在独立的 provider/controller 路径记录，不能通过本地计算缓存保证供应商侧只执行一次。

`usage()` 返回 `new_queries`、`compute_executions`、`global_cache_hits`、`run_reuses`、`technical_retries` 和 `recovered_results`。`cpu_seconds` 与 `wall_seconds` 汇总可测量的本地工具执行时间，包括有计时的失败。`unmeasured_attempts` 表示仍在运行或死于提交前且没有计时产物的尝试数。其他进程、独立 R 阶段和 provider 的时间/费用按各自的运行记录解释。

## 存储要求

使用支持 SQLite POSIX 锁、`flock`、同目录硬链接和 `fsync` 的文件系统。每个进程单独打开连接，使用 SQLite rollback journal 和 FULL 同步；不跨主机宣称缺少锁语义的挂载能安全并发。Slurm 上可选经管理员验证的项目卷，或同节点任务用的持久目录。备份时先让写入进程退出，再复制运行目录、准备快照及配置中引用的文件。冻结结果可单独评价。

## 查询使用轨迹和成员证据

每个数值查询只登记一个 `QueryObservation`。通过另一个修订或显式成员表重读时，`QueryUse` 保存这次选择关联的证据 ID。Agent 的修订和带证据的提交会生成解释节点并链接到已有事实。相同成员来自多个 Seurat 分辨率时，每份 OntoAnno 注释保留自己的分辨率和证据 ID；有冲突时卡片列出备选，提交的单一标签默认为 unknown。

`get_evidence` 传入一个科学证据 ID 和 `offset,limit` 时，返回 `sections`。每个段落含 `total` 与此页的 `rows`。共表达可读 `matched_members`、`gene_counts`，对比可读 `statistics`，扩展可读 `scores`、`signature`、`background`。只在当前运行内校验观察证据、查询和已登记产物三者的绑定，整个操作不执行新的数值计算。多个证据 ID 对应证据卡列表分页。关联到 seed 的扩展深度沿父修订与相同成员集合传播。重新命名、登记同一成员或拆分不能重置已消耗的深度。

恢复曾冻结的 partial 运行时，history 目录保留原冻结清单及只读版本。candidate_tasks.tsv 另存每个初始候选的任务状态。usage 只在所有已派发请求均返回相应用量时给出该项总数，可用 reported_ 前缀字段查看已报告部分。
