# Tentative OntoAnno labels for singleton seeds

In the dashboard choose **注释快照**, select an existing snapshot and set **OntoAnno 单细胞注释** to **注释所有 singleton 和 outlier 种子**. For a full pipeline choose a mode that uses annotations and leave the final annotations input empty to run annotation. If the same snapshot already has a successful annotation, enter its directory under **复用同一快照的成功注释批次**. The task creates its own output directory.

From the project root, with the same R environment, model settings and credentials as the original task:

```bash
rarecell annotate --snapshot "$TASK/snapshot" --config "$TASK/settings.yaml" \
  --include-singletons --reuse-annotations "$TASK/annotations" \
  --out "$TASK/annotations-with-seeds"
```

Set `TASK` to the task directory. Omit `--reuse-annotations` for a first annotation. Alternatively set `annotation.include_singletons: true` in YAML. The default is false. A matching receipt reuses exactly the same snapshot, group composition, panel, prompt, model and batch settings. A completed output with different settings requires a new directory. Supply the combined output as the annotations input to a new discovery run.

Each frozen one-cell seed is considered once, including outliers shared by multiple views. The panel ranks positive differences in mean log-normalized expression against the seed’s frozen neighbourhood, or its complement if that was the stored background. The names of distinct genes sent to OntoAnno are limited by `annotation.top_genes` (default 10). The original prepared marker card also retains negative effects and ranks by magnitude. Single-cell panels have no within-group replication. The upstream model can answer `unknown`. A seed with no positive contrast, no background or no expression is recorded as skipped. Each seed label is an interpretation derived from that seed's RNA observation. Enabling this option neither selects nor merges cells.

Inspect `annotation_status.json` for coverage, skips and errors. `annotations.jsonl` contains the candidate and seed IDs, positive effects, detection fractions, background identity and receipt reference for each annotated singleton. Missing Seurat partition markers are reported separately when singleton annotations can still be produced. The dashboard diagnostics and transcript show the singleton count and the actual requests and responses.

After discovery, compare one gene against reference labels using the original frozen snapshot. This offline check reads the recorded Seurat feature table and each target seed’s exact background and reports its positive and absolute effect ranks.

```bash
python scripts/audit_seed_gene.py --snapshot "$TASK/snapshot" \
  --gt "$GT" --gene "$GENE_ID" --label "$TARGET_LABEL"
```
