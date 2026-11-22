# Discovery 通信协议（soft-protocol-1.4）

模型的一条完整回复包含一个动作，或者一个包含多个独立动作的数组。也接受 `{"actions":[...]}`。控制器提供当前 run、task、数据和契约版本，对回复中显式给出的绑定逐项比对。控制器按当前消融组许可、manifest 与证据登记表执行。解析器只输出结构化请求。

下列两个请求都要求在 H017 的实际成员内，按 log1p 层的默认阈值 0 检查两基因共同检出。H017 必须能在当前 run 中精确解析。

```json
{"action":"coexpress","genes":["CD3D","TRAC"],"candidate_id":"H017"}
```

```json
{"operation":"check co-expression","gene_list":"CD3D; TRAC","candidate":"H017"}
```

## 可调用动作

| 动作 | 请求参数 | 语义和边界 |
|---|---|---|
| `list_candidates` | 可选 `offset,limit` | 候选按稳定顺序分页。offset 从 0 起计，limit 为正整数。 |
| `get_evidence` | 可选 `candidate_id` 或 `evidence_ids`，以及 `offset,limit` | 按候选读成员页与 marker，或按证据 ID 读证据卡。省略选择器时浏览已登记证据。 |
| `compare_groups` | `candidate_id` 或 `members`。可选 `background_id` 或 `background`、`genes` | background 是明确的 cell_id 集合，background_id 引用候选。省略背景时取 foreground 的补集。空对照由科学工具标记。 |
| `coexpress` | `candidate_id`、`members`、`scope="all"` 三选一，`genes` 必填。可选 `logic`、`threshold` | 逻辑为 AND（默认）或 OR。每个指定基因的 log1p 表达须严格大于阈值（默认 0）。 |
| `expand_candidate` / `retrieve_cells` | `candidate_id`。可选 `background_id` 或 `background` | 跨全数据检索，背景只用于对照。一次扩展沿原始候选 lineage 记账。预算、门槛由运行配置决定。 |
| `propose_revision` | `parents,groups`，以及动作或组内的 `reason,evidence_ids` | 每组可给 `candidate_ids`、`members` 或二者的并集，可附 `proposed_label` 和 `status`。每组都需要理由和可验证的证据。来源明确时按下文规则恢复父候选，各组在同一事务中更新。 |
| `submit` | 可选 `candidate_ids,reason,evidence_ids,status` | 省略 candidate_ids 时提交当前选中的修订。显式 `[]` 表示没有选中的候选；剩余细胞交给确定性的回退分配。 |

组的 status 可取 `provisional`、`unresolved`、`confirmed`，生物学名称可以为 `unknown`。submit 的 status 可取 `provisional` 或 `final`。提交意图与运行是否耗尽预算、是否存在外部阻塞分别记录。科学查询可以显式声明固定表达层 `layer="log1p"`。原始 counts 需要通过不同的计算规则处理，因此该层请求会得到局部错误。

修订中的顶层理由记录在 Decision 上，组内理由记录在对应成员版本的解释证据上。顶层引用适用于所有组，组内引用只适用于该组。控制器按精确 ID 检查并去重。当前工作选择为空、提案只有一个组、组中显式给出来源 `candidate_ids` 时，允许省略 `parents` 或填 `[]`。控制器把来源登记在 `Decision.parents_from_candidate_ids`。已有工作选择时，模型须显式指定要替换的父候选。成功的修订更新工作选择，由后续 `submit` 正式提交。

主动查询模式会在每轮上下文的 `recent_expansions` 中优先保留最近查看的扩展候选 ID。通常每项还包括种子 ID、证据 ID 和成员数。控制器从本次运行已保存的调查重建这些引用，读取证据或从检查点恢复后仍可直接提交对应 ID。其余引用在字节预算内提供，最多五项。

每种动作都可显式带 `run_id,task_id,data_version,contract_version,message_id`。这些字段必须与控制器绑定值完全一致。控制器根据候选实际成员、实际基因 ID、背景、逻辑、门槛及影响结果的配置构造科学查询键。消息 ID、字段顺序和说明文字不参与科学查询键。

`expand_candidate`、`coexpress` 和 `compare_groups` 接受可选的 `evidence_ids`，用于注明发起查询的依据。控制器先核对这些 ID 属于当前运行，再把引用记入 `QueryUse.cited_evidence_ids`。复用相同计算时只产生一条数值观测，引用已有证据不会增加独立支持。

读取候选列表或证据表时，控制器按返回的页面内容判断进展。读取同一次科学查询结果中的新一页可以推进调查，再读相同内容不会重置停机计数。模型可以在 `remaining.consecutive_no_progress_turns_before_stop` 中看到还允许连续多少轮无进展。接近停机时，提示会引导模型提交有依据的成员选择，或执行能改变成员判断的新调查。主动发现模式还有查询预算时，停滞反馈也会列出待调查种子的数量和最多三个候选 ID。没有得到模型的有效提交时，运行仍记录为 `partial`。A3_adaptive 的暂定提交见下文。

## 可以本地归一化的输入

- 动作和注册字段支持中英文别名、大小写及空格、下划线和连字符差异。例如 `检查共表达`、`gene_list`、`候选`。实体名称逐字符保留。
- 允许 JSON 围栏、前后说明文字及可唯一修复的尾逗号。解析器只提取一个平衡的 JSON 信封。多个独立动作应写在同一数组内。
- 支持将 `args` 或 `parameters` 对象展开；同名别名在语义归一化后须一致。单个 ID 可以写成一个字符串或单元素数组。基因字符串允许以中英文分号分隔。在 JSON 数组中每个字符串均为完整 ID。AND/OR 的基因及无序成员集按字面值去重排序；组对象的提交顺序保留。
- 明确的数值字符串可转换为有限非负数值。分页参数须为可表示的整数。空值、未知值、布尔数值和无穷大产生错误。省略可选字段仅启用文档中约定的默认值。
- `note,comment,explanation,uncertainty,counterevidence,summary` 接受解释文字并保存在 metadata。用 `args` 或顶层字段明确表达科学条件。未注册字段（包括 metadata 中的 scope/filter/weights）要求澄清；解析器不自行选择条件。
- 在明确的动作名称独占首行时，也接受逐行 `字段：值`。全消息只进行一次本地语义归一化。该受限语法使用与 JSON 完全相同的合法字段和实体检查。

例如同时比较两个群体、再查询一条已有证据可写成一个独立动作数组。

```json
[
  {"action":"compare_groups","candidate":"H017","background_id":"H022","genes":["CD3D"]},
  {"action":"get_evidence","evidence_id":"E001"}
]
```

## 冲突、追溯与恢复

`{"candidate":"H017","scope":"all"}` 同时给出两个范围，会返回 needs_clarification。相似 barcode、证据 ID、候选 ID、gene_id、数据版本不做拼写修复。基因可以使用当前 genes.tsv 的精确 symbol；映射到多个稳定 ID 时要求指定 gene_id。表达“这些细胞”时须给实际候选或明确成员。带符号权重或方向等未注册计算不能转写成无符号共表达。

单条消息允许最多 512,000 个字符、32 层 JSON 和 128 个动作。只有来源为 Discovery Agent 且传输已完成的消息可以进入动作解析。包裹对象解析失败或截断时整条消息不执行。数组内可以解析的各个动作分别检查，单个动作出错只影响该动作。工具返回、引用文章、证据内容和说明字段内的 JSON 不作为后续动作重新解析。

解析器使用本地确定性修复，不调用额外模型。尾逗号和逐行语法修复在解析结果的 metadata.format_repair 中标记。无法消除的歧义返回 `needs_clarification` 和具体字段。重新解释同一条证据与修复通信格式均不增加 observation 数量；执行与重试状态由控制器和状态层记录。逐项冲突测试见 `tests/test_protocol.py`。

相同原始消息及完成标志的本地解析结果存在 ParsedMessage 中。再次出现时复用解析结果，执行前仍检查当前精确 ID、实体版本、模式、谱系和预算。明确截断的远端答复不会产生可执行动作，服务报告的相应用量保留。

真实模型也会将操作放入单层 action 对象，如 `{"action":{"type":"expand_candidate","candidate_id":"H017"}}`。仅在 action 这个确定的容器内，type/name 映射到动作选择符。内外参数同样经过实体、冲突与范围检查；工具产物里的结构不会被提取为命令。

## A3_adaptive

- `expand_candidate` 总是在全体细胞中搜索。省略参照或使用 `background_mode: "local"` 会构造多视图局部参照。`background_mode: "complement"` 明确使用全体非成员作参照。兼容模型以标量 `background: "all"` 表达全局搜索的写法；如果数据中确有名为 `all` 的细胞，它仍是精确的细胞 ID。显式数组始终逐项解析 ID，未知项整项拒绝。回执记录 anchor_id、search_scope 和实际 background_rule。
- 自适应臂返回常规扩展提案，同时保存全体可评分细胞（含种子）的完整排序、各分量以及分数间隙。使用 `get_evidence(evidence_ids=[E], metric="positive_expression", offset=0, limit=10)` 可查看另一个排序；默认按综合分。`score`、`expression`、`positive_expression` 对应综合分、正负 marker 合并表达分和仅正向 marker 表达分。无有效结构或签名时回执说明原因。正向分不可用的细胞不进入该分量的排序。
- `apply_boundary(candidate_id=H, evidence_ids=[E], metric=..., threshold=..., reason=...)` 使用原始锚点 H 对应的扩展观测 E，保留锚点原成员，并加入所选分量达到阈值的所有可评分细胞。阈值在 [0,1] 内，由模型决定；其他分量没有隐含下限。超出配置的 max_cells 新增数时返回 over_limit 和完整命中数，模型可以修订规则。成功的提案为 provisional，需另行选择。复用分数不会新建独立观察或再耗扩展深度。
- `coexpress(scope="all", genes=[...], min_detected=k, threshold=t)` 检索在至少 k 个基因上 log1p 表达严格大于 t 的细胞。k 的范围为去重并解析后的 1 到基因数；k 与 AND/OR 二选一。阈值默认 0。可在 `propose_revision.groups` 中以 `matched_evidence_id` 引用某次共检出的全部成员；可通过 `exclude_candidate_ids` 从成员并集中减掉已存在的候选，再执行原子修订。
- `preview_submission(candidate_ids=[...], candidate_id=H, reference_id=R)` 预览拟提交选择。candidate_id 可填一个未被选中的待查假设，以判断其细胞会进入哪个类别（独立归属、被其他候选认领、冲突或回退）；reference_id 指定参考候选。两个字段均可省略。局部占比以候选与参考集合交集大小除以参考集合大小。实际 submit 按相同成员认领规则记录结果。每个修订群可选填 biological_rarity（yes/no/uncertain）和 rarity_reason，作为模型解释关联到该成员版本。


### 边界比较与停滞提交（discovery/adaptive-2）

首次对某个种子成功执行 `apply_boundary` 时，控制器保存模型指定的成员集。如果它与当前已选群互不重叠，即刻纳入工作草案。控制器对同一指标的非种子细胞排序，在允许的新增人数内寻找更宽边界下方最大的正分数间隙，同分细胞一起纳入。没有可用间隙时只保留模型方案。宽方案作为可选群保存，模型后续指定的其他阈值也可暂存。各方案沿用同一观察 ID。

`boundary_reviews` 比较核心、新增成员、宽方案之外的全部可评分细胞，以及其中指标最高的一小组细胞。后者人数最多等于新增人数，未扩宽时显示一个。表内给出人数及综合分、正向表达分、表达分和结构分的均值。两种成员集各自与其他工作草案组合，分别报告独占成员数、全局冲突数和未认领数。上下文优先保留近期比较，以及最近读过的证据页。可用 `get_evidence(candidate_id=...)` 回看对应的比较。

`select_boundary(candidate_id=..., reason=...)` 替换同一种子的工作边界，保留其他互不重叠的草案。跨群重叠可用 `propose_revision` 分配成员。模型可以用 `submit` 提交当前选择，用显式候选列表指定最终选择，或用空列表放弃全部候选。同一种子的互斥边界只能选一个，除非先修订为互不重叠的群。

### Live 决策轮

所有 live 模式从同一动作许可表生成提示并执行请求。总轮数至少为 3 时保留最后两轮决策，设为 1 或 2 时保留最后一轮。已有充分证据可提前提交。提前停滞时，预算允许即进入一次决策轮。此时可用已有证据修订和提交；自适应臂也可选择已保存的边界，在轮数收尾时还可对已有排名应用阈值。翻页与新科学查询会留下拒绝回执。同批回复中合法的决策仍可执行。提交已选草案可以省略 candidate_ids，明确放弃全部候选请给出空数组。

在停滞或轮数、已报告 token 预算停止时，有非空工作草案则暂定提交，来源标记 controller_stall 或 controller_budget。若模型已进入决策阶段且始终没有提交或留下草案，记录 decision_missing 和 closeout_trigger。服务派发失败继续使用既有的恢复协议。

### 单细胞种子初筛和补查

对 live 主动臂，初筛按冻结的前三个 marker 比较基因 ID 和效应方向，优先较长的共享组合及较多提名视图。每个种子保存一个代表性同伴；可从证据 ID 引用原观察。该顺序只用于安排调查。同一 RNA 的多个视图仍共用证据来源。

探索阶段整轮只请求 list_candidates、至少一次成功且本轮还有动作额度时，控制器可补查一个待查种子。补查使用标准 expand_candidate，受同一查询数、深度和缓存约束，不调用语言模型。budget.max_seed_followups 控制全程上限。开始计算前固定锚点，重启后复用同一计划。逐轮输出中的 controller_seed_review 标出来源。seeds.tsv 保留初筛、实际调查状态以及最终归属。
