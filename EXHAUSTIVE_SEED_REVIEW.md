# 强制逐种子明文调查

在 GUI 的「完整流程」或「从快照发现」中，选择 A2、A3 或 A3_adaptive，再把「单细胞种子调查」设为「强制逐种子调查并写明判断」。默认沿用 YAML，项目默认关闭。也可以在配置顶层设置 `exhaustive_seed_review: true`。对已有任务复用快照和注释时，只需新建一次发现运行。以下命令中的 `old` 指向原任务目录。

```bash
.venv/bin/rarecell discover --snapshot "$old/snapshot" --annotations "$old/annotations" \
  --config "$old/settings.yaml" --mode A3_adaptive --policy live --exhaustive-seed-review \
  --run "$PWD/work/exhaustive-seed-run"
```

启用后，控制器先逐一扩展快照 `seeds.jsonl` 中的单细胞种子，再将测量结果分批给模型。每个种子须有明文判断、处置建议和对应观测。可以建议继续调查、保留、暂缓或排除。模型随后可以综合多个种子的证据，修改成员边界并提交。`seed_reviews.tsv` 保存全部种子，包括中途停止时尚未写完的行；GUI「诊断」显示明文覆盖数，逐轮记录和任务文件提供原文。

默认每批最多 4 个，可用 `budget.seed_review_batch_size` 调整。对 19 个种子，批次装得下且模型一次写完整时，需要 19 次本地逻辑查询和 5 次模型回复，然后才进入常规调查。专项查询额度随种子数增加，批次回复单独计轮，并额外开放一层扩展深度。常规查询、调查轮数和最终决策轮仍然可用。配置中的 `max_reported_tokens` 对整个运行生效。普通模式使用原有预算，命令行可用 `--no-exhaustive-seed-review` 覆盖已启用的 YAML。

需要在评估后核对某一真值类型的上游注释时，可运行 `.venv/bin/python scripts/audit_seed_annotations.py --task "$old" --gt /path/to/gt.tsv --label "$target_label"`。输出按细胞和 V_HVG 分辨率列出所属群、该群的真实 OntoAnno 标签及当时送入模型的 marker；空标签表示缺少对应注释，明确返回的 unknown 按原文显示。
