# OntoAnno labels for each Seurat group and seed

Run annotation on an existing frozen snapshot, then supply its new output directory to discovery:

```bash
.venv/bin/python -m rarecell annotate --snapshot "$TASK/snapshot" --out "$TASK/annotations-dual" --config "$TASK/settings.yaml" --include-singletons
```

The default workflow covers every V_HVG partition group at every resolution and every proposed singleton seed, across all nomination views. An identical member set shares one pair of labels, retaining its seed ID and all partition/resolution links. Existing snapshots work without another prepare. Set `annotation.include_singletons: true` for seed coverage when reusing an older configuration that explicitly turned it off.

* **Global** compares mean log-normalized RNA expression inside the target with the mean across all other input cells.
* **Neighbor** uses the seed’s frozen local neighbours when available. For a group it uses unique outside cells adjacent to its members in the frozen HVG SNN graph. For an isolated group with a valid HVG embedding it uses up to prepare.k external cells nearest the group's centroid by cosine similarity of the unit PCA vectors. A missing neighborhood is recorded explicitly.
* Each contrast ranks positive mean differences across measured genes and keeps up to `annotation.top_genes` distinct symbols (default 10). Ranked effects and detection fractions are retained. Gene order is deterministically shuffled for the model. Neighbor and global requests run separately, even when their input sets coincide. A no-expression, no-background, or no-positive-marker result produces a local unknown with the reason. An unsuccessful request is distinguishable from an upstream model answer of unknown.
* A1, A3 and A3_adaptive expose both labels and gene lists on candidate cards. Active seed triage and exhaustive review also receive them. These are interpretations of the same RNA. No label changes seed membership.
* Inspect `annotations.jsonl` and the run's `annotations.tsv` for both rows per target. Per-panel background rules, sample sizes, marker effects, cell-set hashes, and original model receipts are retained in the annotation directory. `annotation_status.json` distinguishes completed labels, skipped panels and failed batches. Run the same command on a partial directory to reuse successful receipts under the existing retry policy.

For historical experiments set `annotation.dual_backgrounds: false` to use the previous Seurat Wilcoxon input and optional single-seed contrast. Historical annotation directories remain readable by discovery. Use a fresh output directory when changing the panel protocol. Historical single-panel receipts belong to a different request.
