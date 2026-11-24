# 强制逐种子明文调查

在 GUI 的「完整流程」或「从快照发现」中，选择 A2、A3 或 A3_adaptive，再把「单细胞种子调查」设为「强制逐种子调查并写明判断」。默认沿用 YAML，项目默认关闭。也可以在配置顶层设置 `exhaustive_seed_review: true`。对已有任务复用快照和注释时，只需新建一次发现运行。以下命令中的 `old` 指向原任务目录。

```bash
.venv/bin/rarecell discover --snapshot "$old/snapshot" --annotations "$old/annotations" \
  --config "$old/settings.yaml" --mode A3_adaptive --policy live --exhaustive-seed-review \
  --run "$PWD/work/exhaustive-seed-run"
```

启用后，控制器先逐一扩展快照 `seeds.jsonl` 中的单细胞种子，再将测量结果分批给模型。每个种子须有明文判断、处置建议和对应观测。可以建议继续调查、保留、暂缓或排除。模型随后可以综合多个种子的证据，修改成员边界并提交。`seed_reviews.tsv` 保存全部种子，包括中途停止时尚未写完的行；GUI「诊断」显示明文覆盖数，逐轮记录和任务文件提供原文。

默认每批最多 2 个，可用 `budget.seed_review_batch_size` 调整。对 19 个种子，批次装得下且模型一次写完整时，需要 19 次本地逻辑查询和 10 次模型回复，然后才进入常规调查。专项查询额度随种子数增加，批次回复单独计轮，并额外开放一层扩展深度。常规查询、调查轮数和最终决策轮仍然可用。配置中的 `max_reported_tokens` 对整个运行生效。普通模式使用配置中的常规预算，命令行可用 `--no-exhaustive-seed-review` 覆盖已启用的 YAML。

需要在评估后核对某一真值类型的上游注释时，可运行 `.venv/bin/python scripts/audit_seed_annotations.py --task "$old" --gt /path/to/gt.tsv --label "$target_label"`。输出按细胞和 V_HVG 分辨率列出所属群、该群的真实 OntoAnno 标签及当时送入模型的 marker；空标签表示缺少对应注释，明确返回的 unknown 按原文显示。

### 从种子评语交接到候选池调查

在 A3 和 A3_adaptive 的强制逐种子模式中，同一全局 CL 标签形成的种子池，只要有一个成员被标为 `follow_up`，就会在常规调查前得到一次池级调查机会。单细胞扩展加入 0 个细胞也保留该机会。池按种子数从小到大处理，同大小按 CL ID 排序。

每批给模型显示池的 marker 出现次数和书面待查理由。模型可以从观察到的 marker 中选择 `seed_genes`，另外提名正向验证基因 `support_genes` 和竞争谱系基因 `counter_genes`，每组最多 6 个。一次本地查询返回各基因的池/背景差异，以及三组各自从 1 到 m 个基因共检出的种子数和池外细胞数。有原始计数时还返回逐基因 NB 拟合，并按种子的 log2 文库大小分层报告深度匹配的背景检出率与覆盖数。种子提名基因的富集是描述性证据。单独指定验证基因可检查其他信号。

旧式请求只给 `genes` 时，程序补算种子 marker 程序，并列出所提名基因中 effect 为正的一组。若旧式请求显式给出 `genes + min_detected`，沿用全局共检出查询。漏答或基因无法解析时，回退查询保留来源记录。正式做成员修订时可引用 `matched_evidence_id + program + min_detected`，直接取得选定程序在该 k 下的全体成员，也可以用 `matched_evidence_id + nb_high_gene` 取得指定基因的拟合高表达组。

每个全局 CL 池最多测量一次。成员明确支持该名字的池和单独命名引导的种子可使用各自的专项查询，其他池占用常规逻辑查询。交接批次单独计模型轮数，共享全程 token 预算。查询额度用完时，未处理的池保留为 `pending`。每个池的请求来源、待查种子、状态、观测和证据 ID 导出在 `pool_followups.jsonl` 中。`run_summary.json` 的 `pool_followup_coverage` 单列未形成共享 CL 池的待查种子数。

当同一个种子的全局与邻居注释给出不同 CL、已扩展到其他细胞、且种子评语标为 `follow_up`，会单列一次调查。模型收到两个 CL、各自的证据 ID、原评语、扩展候选 ID、成员和邻居基因。候选超过 50 个成员时，卡片显示前 50 个，测量使用完整成员。每个符合条件的扩展最多新增一次专用本地查询，按现有批大小投递给发现模型，共享全程 token 上限。可以在 `pool_followups.jsonl` 按 `kind=conflicting_seed_expansion` 核对。此类记录的 `n_foreground` 以及各程序的 `seed_hits` 对应整个扩展组，其中包含原种子。

普通发现轮优先发送模型主动请求的证据。上下文中的 `requested_evidence` 关联原始回执 ID。较大的结果可压缩显示，在返回可用结果前会保留待交付记录，已提交的判断和完整成员可以追溯。查询候选后，其池级测量优先进入后续调查目录。`get_evidence(candidate_id)` 同样适用于上述独立扩展。

已有 Baron 任务可复用同一任务配套的快照、注释和数值查询缓存。更新后端后，在仓库根目录运行以下命令。新发现运行会重新调用模型做逐种子判断。新运行被暂停或可恢复地中断时，可在相同代码版本下执行 `.venv/bin/rarecell resume --run "$run"`。

```bash
old="$PWD/work/tasks/30fabf9eaefb4b3f96af186a72c72741"
field() {
  .venv/bin/python -c 'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])' "$old/run/checkpoint.json" "$1"
}
run=$(mktemp -d "$old/discovery-handoff.XXXXXX")
echo "$run"
.venv/bin/rarecell discover --config "$old/settings.yaml" --mode A3_adaptive --policy live \
  --exhaustive-seed-review --snapshot "$(field snapshot)" --annotations "$(field annotations)" \
  --cache "$(field cache)" --run "$run"
```

有暂定名字的种子在书面评估中记录 `guide_label` 和 `label_reliability`。agent 用组织、物种和已有的标记知识提出解释，并以该种子的测量说明名字是否可信。首次漏写命名判断会再询问一次，仍缺失时记录为 unassessed，并在覆盖统计中报告。可用的名字引导对应全局池或单独种子的标记查询，允许提议未出现在 OntoAnno 面板里的支持基因和竞争基因。

A3_adaptive 对已测的外扩方案安排单独复核。最短池摘要保留各面板的基因及全部 k 阈值的池内、池外命中数。有新增成员并在核心中富集时，每个正向面板按“命中核心数 /（核心大小 + 方案大小）”展示一个候选，另展示至多两个正向 NB 高表达候选。同一成员集只登记一次。agent 比较核心与新增细胞的支持及竞争程序，选择一个暂定边界、仅保留核心或暂缓，随后仍可修订。种子判断在 `seed_reviews.tsv`，候选与实际复核理由在 `pool_followups.jsonl` 的 `boundaries` 和 `expansion_review`。每批默认两个候选，专项复核单独计轮，仍计入总 token 预算。已保存的完整配置沿用其中指定的预算。使用旧配置新开 discovery 时，可传 `--context-bytes 72000 --seed-review-batch-size 2 --max-agent-turns 18`。已有的可信核心若有成员落在拟合边界外，agent 可显式选择 `include_core: true` 将核心并入。
