"""Paired positive RNA marker panels for frozen Seurat groups and singleton seeds."""
import csv
import re
from collections import Counter, defaultdict

import numpy as np
from scipy import sparse

from .config import fingerprint, read_lines

SOURCE = "positive_mean_log_expression_difference"
PROMPT = "ontoanno-dual-background-v1"
GUIDANCE = ("Marker panels and their OntoAnno labels derive from the same RNA. "
            "Consider competing identities when their evidence differs.")


def _table(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class GroupNeighbours:
    """Outside neighbours in the frozen HVG graph, with an embedding fallback for isolated groups."""

    def __init__(self, snapshot, manifest, cell_index):
        folder = snapshot / "views"
        self.graph = None
        self.embedding = None
        if "V_HVG" not in manifest.get("views", []):
            return
        graph_path = folder / "V_HVG_snn.npz"
        if graph_path.exists():
            self.graph = sparse.load_npz(graph_path).tocsr()
            self.graph.eliminate_zeros()
            self.to_cell = np.array([cell_index[row['cell_id']] for row in _table(folder / 'V_HVG_graph_rows.tsv')], dtype=int)
            self.to_graph = np.full(len(cell_index), -1, dtype=int)
            self.to_graph[self.to_cell] = np.arange(len(self.to_cell))
        path = folder / "V_HVG.npy"
        if path.exists():
            self.embedding = np.load(path, mmap_mode="r", allow_pickle=False)
            self.norms = np.linalg.norm(self.embedding, axis=1)
            self.k = manifest["prepare"]["k"]

    def outside(self, foreground):
        if self.graph is not None:
            rows = self.to_graph[foreground]
            adjacent = self.graph[rows[rows >= 0]].indices
            outside = np.setdiff1d(self.to_cell[adjacent], foreground)
            if len(outside):
                return outside, "HVG_SNN_boundary"
        # An internally connected but isolated cluster still has nearby cells in PCA space.
        if self.embedding is not None:
            valid = self.norms > 1e-10
            members = foreground[valid[foreground]]
            valid[foreground] = False
            others = np.flatnonzero(valid)
            if len(members) and len(others):
                center = np.mean(self.embedding[members] / self.norms[members, None], axis=0)
                if np.linalg.norm(center) > 1e-10:
                    scores = (self.embedding @ center)[others] / self.norms[others]
                    rank = np.lexsort((others, -scores))[:self.k]
                    return np.sort(others[rank]), "HVG_nearest_external_to_centroid"
        return np.array([], dtype=int), "no_local_reference"


def positive_markers(genes, indices, inside, outside, detected_in, detected_out, top):
    effects = inside - outside
    selected, seen = [], set()
    for j in np.lexsort((indices, -effects)):
        if effects[j] <= 1e-12:
            break
        gene = genes[indices[j]]
        symbol = gene.get("symbol") or gene["gene_id"]
        if symbol not in seen:
            seen.add(symbol)
            selected.append(dict(gene_id=gene["gene_id"], symbol=symbol, effect=float(effects[j]),
                                 pct_in=float(detected_in[j]), pct_out=float(detected_out[j])))
            if len(selected) == top:
                break
    return selected


def panels(snapshot, manifest, settings):
    seeds = {s["candidate_id"]: s for s in read_lines(snapshot / 'seeds.jsonl')} if settings.get('include_singletons') else {}
    cell_ids = [row['cell_id'] for row in _table(snapshot / 'cells.tsv')]
    index = {c: i for i, c in enumerate(cell_ids)}
    genes = _table(snapshot / 'genes.tsv')
    x = sparse.load_npz(snapshot / 'data.npz').tocsr().astype(np.float64, copy=False)
    totals = np.asarray(x.sum(axis=0)).ravel()
    detected = x.getnnz(axis=0)
    global_only = settings.get("panel_mode") == "global"
    neighbours = None if global_only else GroupNeighbours(snapshot, manifest, index)
    for row in read_lines(snapshot / 'initial.jsonl'):
        cid = row['id']
        seed = seeds.pop(cid, None)
        partitions = [{k: s[k] for k in ('partition_id', 'resolution')} for s in row['sources']
                      if s['kind'] == 'partition' and s.get('view') == 'V_HVG']
        if not (seed or partitions):
            continue
        members = row['members']
        if seed and members != [seed['cell_id']]:
            raise ValueError('singleton seed membership mismatch')
        fg = np.array([index[c] for c in members], dtype=int)
        inside = x[fg]
        if len(fg) == 1:
            ix, sums, counts = inside.indices, inside.data, np.ones(inside.nnz, dtype=int)
        else:
            sums = np.asarray(inside.sum(axis=0)).ravel()
            ix = np.flatnonzero(sums > 0)
            sums, counts = sums[ix], inside.getnnz(axis=0)[ix]
        common = dict(candidate_id=cid, n_members=len(fg), cell_set_hash=fingerprint(sorted(members)),
                      partitions=partitions, partition_id=partitions[0]['partition_id'] if partitions else None,
                      resolution=partitions[0]['resolution'] if partitions else None, marker_source=SOURCE)
        if seed:
            common['seed_id'] = seed['id']
        bg = row.get('marker_background', {})
        if neighbours is None:
            local, rule = None, None
        elif len(fg) == 1 and bg.get('type') == 'local_neighbours':
            local = np.array(sorted({index[c] for c in bg['members']}), dtype=int)
            if np.intersect1d(fg, local).size:
                raise ValueError('local reference overlaps foreground')
            rule = 'frozen_seed_neighbours'
        else:
            local, rule = neighbours.outside(fg)
        for scope in (('global',) if global_only else ('neighbor', 'global')):
            if scope == 'global':
                n = x.shape[0] - len(fg)
                bg_hash = fingerprint(['complement', manifest['snapshot_hash'], common['cell_set_hash']])
                bg_sum, bg_count = totals[ix] - sums, detected[ix] - counts
                bg_rule = 'all_other_cells'
            else:
                n = len(local)
                bg_hash = fingerprint(sorted(cell_ids[i] for i in local))
                bg_rule = rule
                reference = x[local][:, ix]
                bg_sum, bg_count = np.asarray(reference.sum(axis=0)).ravel(), reference.getnnz(axis=0)
            markers = positive_markers(genes, ix, sums / len(fg), bg_sum / n, counts / len(fg), bg_count / n,
                                settings['top_genes']) if n and len(ix) else []
            # A deterministic permutation keeps ranks out of the annotation prompt and cache resumable.
            symbols = sorted((m['symbol'] for m in markers), key=lambda g: fingerprint([cid, g]))
            panel = dict(common, background_type=scope, background_rule=bg_rule, background_hash=bg_hash,
                         n_background=n, marker_hash=fingerprint(markers), markers=markers, genes=symbols)
            if not symbols:
                panel['reason'] = 'no_expression' if not len(ix) else 'no_background' if not n else 'no_positive_markers'
            yield panel
    if seeds:
        raise ValueError('singleton seed missing from frozen candidates')


def jobs(snapshot, manifest, settings, head):
    batches, skipped = defaultdict(list), []
    for panel in panels(snapshot, manifest, settings):
        if not panel['genes']:
            skipped.append(dict(panel, status='not_annotatable'))
            continue
        subject = 'single cell' if panel['n_members'] == 1 else 'cell group'
        batches[(panel['background_type'], subject)].append(panel)
    result = []
    size = settings['batch_size']
    for (scope, subject), group in sorted(batches.items()):
        background = 'neighboring cells outside it' if scope == 'neighbor' else 'all other cells in the dataset'
        note = (f'Each row describes one {subject}. These are its most positively enriched genes compared with {background}. '
                'Treat each shuffled list as an unordered marker pool. Position carries no weight. '
                'Give one tentative cell type per row. Use unknown when evidence is insufficient.')
        for pos in range(0, len(group), size):
            result.append(dict(head, prompt_version=PROMPT, resolution=None, groups=group[pos:pos+size],
                               system_prompt=head['system_prompt'] + '\n' + note))
    return result, skipped


def ontology_seed_pools(seed_ids, annotations):
    """Group the frozen global annotations by exact CL term, counting a gene once per seed."""
    groups = defaultdict(dict)
    for row in annotations:
        cid = row["candidate_id"]
        term = row.get("ontology_id")
        if (cid not in seed_ids or row.get("background_type") != "global"
                or row.get("status") != "annotated" or not isinstance(term, str) or not re.fullmatch(r"CL:[0-9]{7}", term)
                or row["label"].strip().casefold() == "unknown"):
            continue
        groups[term][cid] = row
    for term, rows in sorted(groups.items()):
        if len(rows) < 2:
            continue
        counts, labels = Counter(), Counter()
        for row in rows.values():
            counts.update(set(row["selection_history"]["genes"]))
            labels[row["label"]] += 1
        label = min(labels, key=lambda name: (-labels[name], name))
        yield term, label, sorted(rows), sorted({r["evidence_id"] for r in rows.values()}), dict(counts)
