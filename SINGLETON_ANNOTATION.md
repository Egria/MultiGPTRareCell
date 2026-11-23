# OntoAnno marker panels and Cell Ontology seed pools

Each distinct Seurat HVG group and each proposed singleton seed gets one global marker panel by default. The contrast compares mean log-normalized expression inside the candidate with the mean across all other cells. Up to ten distinct symbols with the largest positive differences enter the panel (`annotation.top_genes`). OntoAnno receives them once, in a reproducible shuffled order with equal weight. The annotation directory records the effects, detection fractions, original response, CL mapping and any failure or reason for an empty panel.

To use a frozen snapshot and reuse the matching global batches of a completed dual-background annotation run:

```bash
.venv/bin/python -m rarecell annotate --snapshot "$TASK/snapshot" \
  --config "$TASK/settings.yaml" --panel-mode global --include-singletons \
  --reuse-annotations "$TASK/annotations" --out "$TASK/annotations-global"
```

Omit `--reuse-annotations` for a fresh annotation. Reuse requires matching snapshot, model, prompt, batch composition and settings. A new output directory keeps the original experiment reproducible. YAML accepts `annotation.panel_mode: global` or `dual`. The explicit panel mode takes precedence over older configurations that set `dual_backgrounds`. Without a panel mode, an explicit historical `dual_backgrounds: true` requests separate global and neighbor labels, and `false` selects legacy Seurat markers. The dashboard also offers the global and dual modes.

In A1, A3 and A3_adaptive, each exact CL ID assigned to at least two singleton seeds by the global panel creates one provisional candidate. Names that map to the same CL ID contribute to the same pool. Unmapped seeds retain their individual candidates. A Seurat group can receive a label but does not vote in the seed counts. For previously saved dual annotations discovery uses their global rows for pooling.

Ordinary discovery turns show the pool IDs, sizes and most frequent genes. `list_ontology_pools` pages the full index. `get_evidence(candidate_id, offset, limit)` pages the seed IDs, cell IDs and the union of marker genes with counts. Each count is the number of distinct member seeds whose global panel included that gene. It describes panel selection, not the expression frequency across all cells. Active discovery can use the pool as a foreground for comparison and query global marker expression to find further members. The agent can select, revise or leave out the pool.
