"""Local focal gene annotation, ORA and gene/pathway Jaccard analysis.

Inputs are independently annotated CNA tables; no CNA alignment is performed.
Coordinates are 1-based inclusive. Only autosomal genes (chr1-22) are used.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import sys
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

LIBRARIES = ['MSigDB_Hallmark_2020', 'Reactome_2022', 'KEGG_2021_Human']
KEY = ['sample', 'caller', 'bin_size', 'hcCNA_type']
GENES = {}
AUTOSOMES = {f'chr{i}' for i in range(1, 23)}
DEFAULT_RELEASES = {'hg19': 'v19', 'hg38': 'v46'}
ASSEMBLIES = {'hg19': 'GRCh37', 'hg38': 'GRCh38'}

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_csv(rows, path, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows, columns=columns)
    frame.to_csv(path, index=False)


def parse_gtf(path, output, genome="hg38"):
    """Read gene features from a local plain/gzipped GENCODE GTF."""
    rows = []
    headers = []
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                headers.append(line.rstrip())
                assemblies = set(re.findall(r'GRCh(?:37|38)|hg(?:19|38)', line))
                allowed = {genome, ASSEMBLIES[genome]}
                if assemblies and not assemblies <= allowed:
                    raise ValueError(f'GTF assembly header conflicts with {genome}: {line.strip()}')
                continue
            if not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9:
                raise ValueError(f'Invalid GTF row in {path}')
            if fields[2] != 'gene' or fields[0] not in AUTOSOMES:
                continue
            a = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            rows.append((fields[0], int(fields[3]), int(fields[4]), a['gene_id'], a['gene_name'], a['gene_type']))
    genes = pd.DataFrame(rows, columns=['Chrom', 'Start', 'End', 'gene_id', 'gene_name', 'gene_type'])
    if genes.empty:
        raise ValueError(f'No chr1-22 gene features in {path}')
    genes.attrs['gtf_header'] = headers
    write_csv(genes, output)
    return genes


def init_worker(gene_path):
    global GENES
    genes = pd.read_csv(gene_path)
    GENES = {chrom: sub.sort_values('Start').reset_index(drop=True) for chrom, sub in genes.groupby('Chrom')}


def overlaps(chrom, start, end):
    g = GENES.get(chrom)
    if g is None:
        return pd.DataFrame(columns=['Chrom', 'Start', 'End', 'gene_id', 'gene_name', 'gene_type'])
    right = np.searchsorted(g.Start.to_numpy(), end, side='right')
    sub = g.iloc[:right]
    return sub.loc[sub.End.to_numpy() >= start]


def process_file(job):
    source = Path(job['source'])
    out = Path(job['output'])
    sample, caller, bin_size = job['sample'], job['caller'], job['bin_size']
    mode, resume = job['mode'], job['resume']
    rel = Path(sample) / caller / bin_size / job['file_id']
    cache = out / 'file_cache' / rel.with_suffix('.json')
    signature = job['signature']
    dest = out / 'annotated' / rel.with_suffix('.annotated.csv.gz')
    event_path = out / 'events' / rel.with_suffix('.focal_events.csv.gz')
    if resume and cache.exists() and dest.exists() and event_path.exists():
        saved = json.loads(cache.read_text())
        if saved['signature'] == signature:
            return str(cache)
    dest.parent.mkdir(parents=True, exist_ok=True)
    event_path.parent.mkdir(parents=True, exist_ok=True)
    chromosomes = set()
    interval_cache = {}
    focal_loci = set()
    events = {}  # (chrom, cell, Segment, type) -> [number of input rows, set of covered intervals]
    all_cells, callable_genes, group_genes = set(), set(), defaultdict(set)
    row_counts, type_counts, focal_counts = Counter(), Counter(), Counter()
    input_rows = 0
    first = True
    with gzip.open(dest, 'wt', compresslevel=1, newline='') as output:
        for d in pd.read_csv(source, chunksize=200_000, keep_default_na=False, dtype={'cell': str, 'Chrom': str, 'Segment': str, 'size': str, 'type': str}):
            expected = ['cell', 'Chrom', 'Start', 'End', 'value', 'Segment', 'size', 'type']
            if not set(expected).issubset(d.columns):
                raise ValueError(f'Unexpected columns: {source}: {list(d.columns)}')
            if (not d.Chrom.isin(AUTOSOMES).all() or (d.Start < 1).any()
                    or (d.End < d.Start).any() or (d.Start % 1 != 0).any() or (d.End % 1 != 0).any()):
                raise ValueError(f'Invalid coordinates in {source}')
            for typ in d['type'].unique():
                _label(typ, 'hcCNA type')
            if d['type'].eq('ALL_TYPES').any():
                raise ValueError('ALL_TYPES is reserved for the pooled output')
            chromosomes.update(d.Chrom)
            d['source_row'] = np.arange(input_rows + 2, input_rows + len(d) + 2)
            input_rows += len(d)
            all_cells.update(d.cell)
            row_counts.update(d['size'])
            type_counts.update(d['type'])
            for chrom, start, end in d[['Chrom', 'Start', 'End']].drop_duplicates().itertuples(index=False, name=None):
                ik = (chrom, start, end)
                if ik not in interval_cache:
                    hits = overlaps(chrom, start, end)
                    all_symbols = set(hits.gene_name)
                    ids = set(hits.gene_id)
                    pc = set(hits.loc[hits.gene_type.eq('protein_coding'), 'gene_name'])
                    interval_cache[ik] = (all_symbols, ids, pc)
                callable_genes.update(interval_cache[ik][2])
            f = d.loc[d['size'].str.lower().eq('focal')]
            focal_counts.update(f['type'])
            for chrom, cell, segment, typ, start, end in f[['Chrom', 'cell', 'Segment', 'type', 'Start', 'End']].itertuples(index=False, name=None):
                m = re.fullmatch(r'(chr\d+):(\d+)-(\d+)', segment)
                if not m or m[1] != chrom or start < int(m[2]) or end > int(m[3]):
                    raise ValueError(f'Row outside Segment: {source}: {segment} {start}-{end}')
                ik = (chrom, start, end)
                ek = (chrom, cell, segment, typ)
                if ek not in events:
                    events[ek] = [0, set()]
                events[ek][0] += 1
                events[ek][1].add(ik)
                focal_loci.add((typ, chrom, start, end))
                group_genes[typ].update(interval_cache[ik][2])
            target = d.copy() if mode == 'full' else f.copy()
            for col in ['focal_genes', 'focal_gene_ids', 'focal_protein_coding_genes']:
                target[col] = ''
            target['n_focal_genes'] = 0
            target['n_focal_protein_coding_genes'] = 0
            # Annotate the actual row span, not the Segment envelope.
            if len(f):
                tuples = [interval_cache[x] for x in f[['Chrom', 'Start', 'End']].itertuples(index=False, name=None)]
                for col, j in [('focal_genes', 0), ('focal_gene_ids', 1), ('focal_protein_coding_genes', 2)]:
                    target.loc[f.index, col] = [';'.join(sorted(t[j])) for t in tuples]
                target.loc[f.index, 'n_focal_genes'] = [len(t[0]) for t in tuples]
                target.loc[f.index, 'n_focal_protein_coding_genes'] = [len(t[2]) for t in tuples]
            target.to_csv(output, index=False, header=first)
            first = False
    group_event_counts, group_event_gene_counts = Counter(), defaultdict(Counter)
    group_cells = defaultdict(set)
    gap_events = 0
    with gzip.open(event_path, 'wt', compresslevel=1, newline='') as handle:
        w = csv.writer(handle)
        w.writerow(['cell', 'Chrom', 'Segment', 'hcCNA_type', 'n_input_rows', 'n_unique_intervals',
                    'covered_bp', 'segment_bp', 'coverage_complete', 'focal_genes', 'focal_protein_coding_genes'])
        for (chrom, cell, segment, typ), (nrows, intervals) in sorted(events.items()):
            all_symbols, pc = set(), set()
            for ik in intervals:
                all_symbols.update(interval_cache[ik][0])
                pc.update(interval_cache[ik][2])
            merged = []
            for _, start, end in sorted(intervals):
                if merged and start <= merged[-1][1] + 1:
                    merged[-1][1] = max(end, merged[-1][1])
                else:
                    merged.append([start, end])
            covered = sum(b - a + 1 for a, b in merged)
            a, b = map(int, segment.split(':')[1].split('-'))
            complete = covered == b - a + 1
            gap_events += not complete
            w.writerow([cell, chrom, segment, typ, nrows, len(intervals), covered, b - a + 1, complete,
                        ';'.join(sorted(all_symbols)), ';'.join(sorted(pc))])
            group_event_counts[typ] += 1
            group_event_gene_counts[typ].update(pc)
            group_cells[typ].add(cell)
    summary = {
        'signature': signature, 'source': str(source), 'sample': sample, 'caller': caller, 'bin_size': bin_size,
        'chromosome': ';'.join(sorted(chromosomes)), 'chromosomes': sorted(chromosomes), 'sha256': signature['source_sha256'], 'input_rows': input_rows,
        'size_counts': dict(row_counts), 'type_counts': dict(type_counts), 'focal_row_counts': dict(focal_counts),
        'cells': sorted(all_cells), 'callable_genes': sorted(callable_genes),
        'group_genes': {t: sorted(g) for t, g in group_genes.items()},
        'group_event_counts': dict(group_event_counts),
        'group_event_gene_counts': {t: dict(c) for t, c in group_event_gene_counts.items()},
        'group_cells': {t: sorted(c) for t, c in group_cells.items()},
        'unique_focal_loci_counts': dict(Counter(t for t, chrom, a, b in focal_loci)),
        'incomplete_segment_events': int(gap_events), 'annotated_file': str(dest.relative_to(out)),
        'event_file': str(event_path.relative_to(out)),
    }
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(summary))
    return str(cache)


def bh(p):
    p = np.asarray(p, float)
    if not len(p):
        return p
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    ans = np.empty(len(p))
    ans[order] = np.minimum(q, 1)
    return ans


def read_gmt(path):
    terms = {}
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        if len(fields) >= 3:
            terms[fields[0]] = {s.split(',')[0].strip() for s in fields[2:] if s.strip()}
    if not terms:
        raise ValueError(f'Empty GMT: {path}')
    return terms


def metrics(a, b):
    union, inter = a | b, a & b
    return {'n_a': len(a), 'n_b': len(b), 'intersection': len(inter), 'union': len(union),
            'jaccard': len(inter) / len(union) if union else np.nan,
            'overlap_coefficient': len(inter) / min(len(a), len(b)) if a and b else np.nan,
            'retained_from_a': len(inter) / len(a) if a else np.nan,
            'both_empty': not union}


def analyze(records, args):
    os.environ.setdefault('MPLCONFIGDIR', str(args.output / 'mpl_cache'))
    import gseapy as gp
    types = sorted(set().union(*(set(r['type_counts']) for r in records)))
    grouped = defaultdict(list)
    for r in records:
        grouped[(r['sample'], r['caller'], r['bin_size'])].append(r)
    observed, cells = {}, {}
    for k, rr in grouped.items():
        observed[k] = set().union(*(set(r['callable_genes']) for r in rr))
        cells[k] = set().union(*(set(r['cells']) for r in rr))
    samples = list(dict.fromkeys(k[0] for k in grouped))
    settings = list(dict.fromkeys((s, b) for s, c, b in grouped))
    callers_by_setting = {(s, b): [c for ss, c, bb in grouped if (ss, bb) == (s, b)] for s, b in settings}
    resolutions_by_caller = {(s, c): [b for ss, cc, b in grouped if (ss, cc) == (s, c)]
                             for s, c, b in grouped}
    shared_observed = {s: set.intersection(*(gs for k, gs in observed.items() if k[0] == s)) for s in samples}
    gene_table = pd.read_csv(args.gene_path)
    genome_genes = set(gene_table.loc[gene_table.gene_type.eq('protein_coding'), 'gene_name'])
    backgrounds = {s: genome_genes.copy() for s in samples} if args.background == 'genome' else shared_observed
    write_csv([{'sample': s, 'gene': g} for s, gs in shared_observed.items() for g in sorted(gs)],
              args.output / 'shared_observed_coordinate_genes.csv', columns=['sample', 'gene'])
    bgrows = [{'sample': s, 'gene': g} for s, gs in backgrounds.items() for g in sorted(gs)]
    write_csv(bgrows, args.output / 'background_genes.csv', columns=['sample', 'gene'])
    library_sets = {lib: read_gmt(path) for lib, path in args.gene_sets.items()}
    group_rows, genes_long, all_ora = [], [], []
    gene_sets = {}
    for s, c, b in grouped:
        rr = grouped[(s, c, b)]
        for typ in types:
            k = (s, c, b, typ)
            raw = set().union(*(set(r['group_genes'].get(typ, [])) for r in rr))
            query = raw & backgrounds[s]
            gene_sets[k] = query
            event_counts = Counter()
            for r in rr:
                event_counts.update(r['group_event_gene_counts'].get(typ, {}))
            n_events = sum(r['group_event_counts'].get(typ, 0) for r in rr)
            row = dict(zip(KEY, k))
            row.update(n_input_rows=sum(r['input_rows'] for r in rr), n_cells=len(cells[(s, c, b)]),
                       n_focal_rows=sum(r['focal_row_counts'].get(typ, 0) for r in rr), n_focal_events=n_events,
                       n_focal_cells=len(set().union(*(set(r['group_cells'].get(typ, [])) for r in rr))),
                       n_focal_genes_before_background=len(raw), n_focal_genes=len(query),
                       n_background=len(backgrounds[s]), n_excluded_from_background=len(raw - query),
                       focal_genes=';'.join(sorted(query)))
            group_rows.append(row)
            for g in sorted(raw):
                genes_long.append({**dict(zip(KEY, k)), 'gene': g, 'in_common_background': g in query,
                                   'n_supporting_cell_events': event_counts[g]})
            # Every group/library has an explicit status, including no focal events.
            for lib, terms in library_sets.items():
                eligible = {term: sorted(gs & backgrounds[s]) for term, gs in terms.items()
                            if args.min_set_size <= len(gs & backgrounds[s]) <= args.max_set_size}
                status = 'tested' if len(query) >= args.min_query_size else ('no_focal_genes' if not query else 'too_few_genes')
                if not eligible:
                    status = 'no_eligible_terms'
                if status != 'tested':
                    all_ora.append({**dict(zip(KEY, k)), 'library': lib, 'term': '', 'status': status,
                                    'query_size': len(query), 'background_size': len(backgrounds[s]),
                                    'n_tested_terms': 0, 'n_eligible_terms': len(eligible)})
                    continue
                # Local GSEApy ORA; no gene lists are sent to any external service.
                # GSEApy omits zero-overlap terms: restore them with p=1 before BH.
                enr = gp.enrich(gene_list=sorted(query), gene_sets=eligible, background=sorted(backgrounds[s]),
                                outdir=None, no_plot=True, cutoff=1, verbose=False)
                result = enr.results
                # GSEApy 1.1.9 leaves results as [] when no term has an overlap.
                # Such terms are still valid tests with p=1, not missing data.
                if result is None or (isinstance(result, list) and not result):
                    lookup = {}
                elif isinstance(result, pd.DataFrame):
                    lookup = {} if result.empty else result.set_index('Term')['P-value'].to_dict()
                else:
                    raise TypeError(f'Unexpected GSEApy result type: {type(result)}')
                library_rows = []
                for term, gs in eligible.items():
                    hits = query & set(gs)
                    pval = float(lookup.get(term, 1.0))
                    expected_p = float(hypergeom.sf(len(hits) - 1, len(backgrounds[s]), len(gs), len(query)))
                    if not np.isclose(pval, expected_p, rtol=1e-8, atol=1e-14):
                        raise AssertionError(f'GSEApy/hypergeometric mismatch: {k}, {term}')
                    library_rows.append({**dict(zip(KEY, k)), 'library': lib, 'term': term, 'status': status,
                        'query_size': len(query), 'background_size': len(backgrounds[s]),
                        'n_tested_terms': len(eligible), 'n_eligible_terms': len(eligible),
                        'pathway_size': len(gs), 'overlap_count': len(hits), 'overlap_genes': ';'.join(sorted(hits)),
                        'gene_ratio': len(hits) / len(query), 'pathway_coverage': len(hits) / len(gs),
                        'fold_enrichment': (len(hits) / len(query)) / (len(gs) / len(backgrounds[s])), 'p_value': pval})
                qvals = bh([r['p_value'] for r in library_rows])
                for r, q in zip(library_rows, qvals):
                    r['fdr_library'] = q
                all_ora.extend(library_rows)
    ora_columns = KEY + ['library', 'term', 'status', 'query_size', 'background_size',
        'n_tested_terms', 'n_eligible_terms', 'pathway_size', 'overlap_count', 'overlap_genes',
        'gene_ratio', 'pathway_coverage', 'fold_enrichment', 'p_value', 'fdr_library']
    ora = pd.DataFrame(all_ora, columns=ora_columns)
    valid = ora.status.eq('tested')
    ora['fdr_global'] = np.nan
    ora.loc[valid, 'fdr_global'] = bh(ora.loc[valid, 'p_value'])
    ora['fdr_group'] = np.nan
    for _, sub in ora.loc[valid].groupby(KEY, sort=False):
        ora.loc[sub.index, 'fdr_group'] = bh(sub.p_value)
    ora['significant'] = (ora.fdr_group < args.fdr) & (ora.overlap_count >= args.min_overlap)
    summary = pd.DataFrame(group_rows, columns=KEY + ['n_input_rows', 'n_cells', 'n_focal_rows',
        'n_focal_events', 'n_focal_cells', 'n_focal_genes_before_background', 'n_focal_genes',
        'n_background', 'n_excluded_from_background', 'focal_genes'])
    sigcounts = ora.loc[ora.significant].groupby(KEY).size().rename('n_significant_pathways').reset_index()
    summary = summary.merge(sigcounts, on=KEY, how='left')
    summary['n_significant_pathways'] = summary.n_significant_pathways.fillna(0).astype(int)
    write_csv(summary, args.output / 'group_summary.csv')
    write_csv(genes_long, args.output / 'focal_genes_long.csv', columns=KEY + ['gene', 'in_common_background', 'n_supporting_cell_events'])
    write_csv(ora, args.output / 'pathways_all_tests.csv.gz')
    write_csv(ora.loc[ora.significant], args.output / 'focal_pathways_significant.csv')
    for k, sub in ora.groupby(KEY):
        path = args.output / 'pathways_by_group' / k[0] / k[1] / k[2] / (k[3] + '.csv')
        write_csv(sub, path)
    statuses = ora.groupby(KEY + ['library', 'status'], dropna=False).agg(
        query_size=('query_size', 'first'), n_eligible_terms=('n_eligible_terms', 'first'),
        n_tested_terms=('n_tested_terms', 'first'), n_significant=('significant', 'sum')).reset_index()
    write_csv(statuses, args.output / 'enrichment_status.csv')
    path_sets = {k: set() for k in gene_sets}
    for k, sub in ora.loc[ora.significant].groupby(KEY):
        path_sets[k] = set(zip(sub.library, sub.term))
    comparisons = []
    for kind, sets in [('gene', gene_sets), ('pathway', path_sets)]:
        for s, b, typ in ((s, b, typ) for s, b in settings for typ in types + ['ALL_TYPES']):
            for a, c in itertools.combinations(callers_by_setting[(s, b)], 2):
                aa = sets[(s, a, b, typ)] if typ != 'ALL_TYPES' else set().union(*(sets[(s, a, b, t)] for t in types))
                cc = sets[(s, c, b, typ)] if typ != 'ALL_TYPES' else set().union(*(sets[(s, c, b, t)] for t in types))
                comparisons.append({'comparison': 'between_callers', 'feature': kind, 'sample': s, 'bin_size': b,
                                    'hcCNA_type': typ, 'caller_a': a, 'caller_b': c, **metrics(aa, cc)})
        for (s, c), bins in resolutions_by_caller.items():
            for typ in types + ['ALL_TYPES']:
                for bin_a, bin_b in itertools.combinations(bins, 2):
                    aa = sets[(s, c, bin_a, typ)] if typ != 'ALL_TYPES' else set().union(*(sets[(s, c, bin_a, t)] for t in types))
                    bb = sets[(s, c, bin_b, typ)] if typ != 'ALL_TYPES' else set().union(*(sets[(s, c, bin_b, t)] for t in types))
                    comparisons.append({'comparison': 'between_bin_sizes', 'feature': kind, 'sample': s, 'caller_a': c,
                                        'caller_b': c, 'bin_a': bin_a, 'bin_b': bin_b, 'hcCNA_type': typ, **metrics(aa, bb)})
    # Preserve the reference's first-seen column order for nonempty comparisons.
    write_csv(comparisons, args.output / 'consistency_jaccard.csv', columns=None if comparisons else
              ['comparison', 'feature', 'sample', 'bin_size', 'hcCNA_type', 'caller_a', 'caller_b',
               'n_a', 'n_b', 'intersection', 'union', 'jaccard', 'overlap_coefficient', 'retained_from_a',
               'both_empty', 'bin_a', 'bin_b'])
    support = []
    for kind, sets in [('gene', gene_sets), ('pathway', path_sets)]:
        for s, b, typ in ((s, b, typ) for s, b in settings for typ in types + ['ALL_TYPES']):
            by_caller = {c: sets[(s, c, b, typ)] if typ != 'ALL_TYPES' else set().union(*(sets[(s, c, b, t)] for t in types)) for c in callers_by_setting[(s, b)]}
            for feature in sorted(set().union(*by_caller.values())):
                callers = [c for c in callers_by_setting[(s, b)] if feature in by_caller[c]]
                support.append({'feature': kind, 'sample': s, 'bin_size': b, 'hcCNA_type': typ,
                                'feature_id': feature if kind == 'gene' else feature[0] + '::' + feature[1],
                                'n_callers': len(callers), 'callers': ';'.join(callers)})
    support_frame = pd.DataFrame(support, columns=['feature', 'sample', 'bin_size', 'hcCNA_type', 'feature_id', 'n_callers', 'callers'])
    write_csv(support_frame, args.output / 'caller_support.csv.gz')
    write_csv(support_frame.loc[support_frame.feature.eq('gene') & support_frame.n_callers.ge(3)],
              args.output / 'consensus_genes_at_least_3_callers.csv')
    write_csv(support_frame.loc[support_frame.feature.eq('pathway') & support_frame.n_callers.ge(2)],
              args.output / 'shared_pathways_at_least_2_callers.csv')
    pooled = []
    for s, c, b in grouped:
        ss = summary.loc[summary['sample'].eq(s) & summary.caller.eq(c) & summary.bin_size.eq(b)]
        pooled.append({'sample': s, 'caller': c, 'bin_size': b,
                       'n_focal_events': int(ss.n_focal_events.sum()),
                       'n_genes': len(set().union(*(gene_sets[(s, c, b, t)] for t in types))),
                       'n_pathways': len(set().union(*(path_sets[(s, c, b, t)] for t in types)))})
    write_csv(pooled, args.output / 'all_types_counts.csv')
    write_csv([{**{k: r[k] for k in ['sample', 'caller', 'bin_size', 'chromosome']}, 'size_label': size, 'n_rows': n}
               for r in records for size, n in r['size_counts'].items()], args.output / 'size_label_audit.csv', columns=['sample', 'caller', 'bin_size', 'chromosome', 'size_label', 'n_rows'])
    input_audit = [{k: r[k] for k in ['source', 'sample', 'caller', 'bin_size', 'chromosome', 'sha256',
                   'input_rows', 'incomplete_segment_events', 'annotated_file', 'event_file']} for r in records]
    write_csv(input_audit, args.output / 'input_audit.csv')
    # Match complete files only by sample/caller and the same chromosome coverage.
    # Never hash reconstructed, concatenated or sorted CSV contents.
    aa = pd.DataFrame(input_audit)
    aa = aa.loc[aa.chromosome.ne('')]
    pairs = aa.pivot(index=['sample', 'caller', 'chromosome'], columns='bin_size', values='sha256')
    pairs = pairs.reindex(columns=sorted({r['bin_size'] for r in records}))
    matched = pairs.notna().sum(axis=1) >= 2
    pairs = pairs.loc[matched]
    identical = pairs.nunique(axis=1).eq(1)
    # Arbitrary resolution labels may collide with the reference's metadata fields.
    reserved = {'sample', 'caller', 'chromosome', 'byte_identical'}
    renamed = {}
    for b in pairs.columns:
        if b in reserved:
            name = 'sha256:' + b
            while name in pairs.columns or name in renamed.values():
                name = 'sha256:' + name
            renamed[b] = name
    pairs = pairs.rename(columns=renamed)
    pairs = pairs.reset_index()
    pairs['byte_identical'] = identical.to_numpy()
    write_csv(pairs, args.output / 'resolution_input_identity.csv')
    cellrows = []
    for s in samples:
        common = set.intersection(*(gs for k, gs in cells.items() if k[0] == s))
        for _, c, b in (k for k in grouped if k[0] == s):
            cellrows.append({'sample': s, 'caller': c, 'bin_size': b, 'n_cells': len(cells[(s, c, b)]),
                             'n_cells_shared_all_settings': len(common), 'cells_shared_all_settings': ';'.join(sorted(common))})
    write_csv(cellrows, args.output / 'cell_coverage_audit.csv')
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(), 'assembly': ASSEMBLIES[args.genome], 'genome': args.genome,
                'assembly_basis': 'Explicit API genome; input CNA assembly must match' ,
                'gtf': str(args.gtf), 'gtf_sha256': sha256(args.gtf), 'gtf_version': args.annotation_release, 'gtf_header': args.gtf_header,
                'coordinate_convention': '1-based inclusive', 'focal_rule': 'input size label equals focal',
                'annotation_mode': args.annotation_mode, 'ora': 'local GSEApy enrich; hypergeometric ORA, not ranked GSEA',
                'gene_scope': 'all biotypes annotated; unique protein-coding symbols tested',
                'background_mode': args.background,
                'background': ('all protein-coding symbols on chromosomes 1-22 from the selected GTF; fixed across all samples/callers/resolutions/types'
                               if args.background == 'genome' else 'sensitivity only: per-sample intersection of genes overlapping observed CNA coordinates across all supplied settings of that sample'),
                'significance': f'BH across all eligible terms from all libraries within each sample/caller/bin/type: FDR < {args.fdr}; overlap >= {args.min_overlap}',
                'min_query_size': args.min_query_size, 'min_set_size': args.min_set_size, 'max_set_size': args.max_set_size,
                'empty_set_jaccard': 'NA when both sets empty, zero when exactly one empty',
                'libraries': {lib: {'sha256': sha256(args.gene_sets[lib]), 'n_terms': len(library_sets[lib]),
                  'source': str(args.gene_sets[lib])} for lib in args.gene_sets},
                'versions': {'python': sys.version, 'pandas': pd.__version__, 'numpy': np.__version__, 'gseapy': gp.__version__},
                'inputs': args.inputs, 'analysis_signature': args.run_signature,
                'n_files': len(records), 'n_rows': sum(r['input_rows'] for r in records)}
    (args.output / 'run_manifest.json').write_text(json.dumps(manifest, indent=2))
    print(summary.groupby(['sample', 'caller', 'bin_size'])[['n_focal_events', 'n_significant_pathways']].sum().to_string(), flush=True)
    return {'output_dir': str(args.output), 'group_summary': summary, 'enrichment_status': statuses,
            'consistency_jaccard': pd.read_csv(args.output / 'consistency_jaccard.csv'),
            'manifest': str(args.output / 'run_manifest.json')}


def _label(value, field):
    """Validate labels also used as output path components."""
    if not isinstance(value, str) or not value or value in {'.', '..'} or any(c in value for c in '/\\\x00'):
        raise ValueError(f'{field} must be a nonempty string without path separators: {value!r}')
    return value


def _normalize_inputs(annotated_inputs, output):
    inputs = []
    seen_groups, seen_files = set(), set()
    for record in annotated_inputs:
        if 'comparison_group' in record:
            raise ValueError('Use independently annotated CNA, not pairwise comparison_group inputs')
        sample = _label(record['sample'], 'sample')
        caller = _label(record['tool'], 'tool')
        bin_size = _label(record['bin_size'], 'bin_size')
        key = (sample, caller, bin_size)
        if key in seen_groups:
            raise ValueError(f'Duplicate input group {key}; put its chromosome files in one files list')
        seen_groups.add(key)
        files = record['files']
        if isinstance(files, (str, Path)):
            files = [files]
        files = [Path(f).expanduser().resolve() for f in files]
        if not files:
            raise ValueError(f'No input files for {key}')
        for path in files:
            if not path.is_file():
                raise FileNotFoundError(path)
            if path in seen_files:
                raise ValueError(f'Input file supplied more than once: {path}')
            if path == output or output in path.parents:
                raise ValueError('Input files must be outside the analysis output directory')
            seen_files.add(path)
        inputs.append({'sample': sample, 'caller': caller, 'bin_size': bin_size,
                       'files': [str(p) for p in files]})
    if not inputs:
        raise ValueError('annotated_inputs must contain at least one group')
    return inputs


def run_focal_gene_pathway_analysis(
    annotated_inputs, output_dir, *, genome='hg38', resource_dir=None,
    gtf_path=None, annotation_release=None, gene_sets=None, libraries=None,
    workers=2, annotation_mode='full', resume=False, background='genome',
    min_query_size=3, min_set_size=10, max_set_size=500, min_overlap=3, fdr=0.05,
):
    """Analyze independent annotated CNA files using local GTF and GMT resources.

    Each input record has sample, tool, bin_size and files (one path or a list
    of disjoint chromosome files). No alignment or CNA classification occurs.
    Default GTF: genes/hg19/gencode.v19.annotation.gtf.gz or
    genes/hg38/gencode.v46.annotation.gtf.gz under resource_dir
    (default ~/.cache/hcbench/resources). GMTs default to bundled lib/pathways;
    an explicit resource_dir selects its pathways/<library>.gmt instead.
    Custom GTF versions must
    be named explicitly with annotation_release. Downloads are never performed.
    Existing output can only be resumed with the same inputs and parameters.
    """
    genome = {'GRCh37': 'hg19', 'GRCh38': 'hg38'}.get(genome, genome)
    if genome not in ASSEMBLIES:
        raise ValueError('genome must be hg19/GRCh37 or hg38/GRCh38')
    if annotation_mode not in {'full', 'focal'} or background not in {'genome', 'shared-observed'}:
        raise ValueError('Invalid annotation_mode or background')
    if not isinstance(workers, int) or workers < 1:
        raise ValueError('workers must be a positive integer')
    if not (0 < fdr <= 1) or not (1 <= min_set_size <= max_set_size) or min_query_size < 1 or min_overlap < 1:
        raise ValueError('Invalid ORA thresholds')
    output = Path(output_dir).expanduser().resolve()
    inputs = _normalize_inputs(annotated_inputs, output)
    resources = Path(resource_dir or '~/.cache/hcbench/resources').expanduser().resolve()
    if gtf_path is None:
        expected_release = DEFAULT_RELEASES[genome]
        if annotation_release not in (None, expected_release, 'GENCODE ' + expected_release):
            raise ValueError('Nondefault annotation_release requires an explicit gtf_path')
        annotation_release = 'GENCODE ' + expected_release
        gtf = resources / 'genes' / genome / f'gencode.{expected_release}.annotation.gtf.gz'
        if not gtf.exists() and gtf.with_suffix('').exists():
            gtf = gtf.with_suffix('')
    else:
        gtf = Path(gtf_path).expanduser().resolve()
        if not isinstance(annotation_release, str) or not annotation_release.strip():
            raise ValueError('Provide annotation_release for a custom gtf_path')
    if gene_sets is None:
        libraries = LIBRARIES if libraries is None else libraries
        if isinstance(libraries, str):
            raise ValueError('libraries must be a sequence of names')
        pathway_dir = resources / 'pathways' if resource_dir is not None else Path(__file__).resolve().parent / 'lib' / 'pathways'
        gene_sets = {_label(lib, 'library'): pathway_dir / (lib + '.gmt') for lib in libraries}
    elif libraries is not None:
        raise ValueError('Specify either gene_sets or libraries, not both')
    gene_sets = {_label(lib, 'library'): Path(path).expanduser().resolve() for lib, path in gene_sets.items()}
    if not gene_sets:
        raise ValueError('At least one local GMT library is required')
    for path in [gtf, *gene_sets.values()]:
        if not path.is_file():
            raise FileNotFoundError(f'Missing local resource: {path}; see docs/realbench/focal_genes_and_pathways.md')
        if path == output or output in path.parents:
            raise ValueError('Reference resources must be outside the analysis output directory')
    gtf_hash = sha256(gtf)
    module_hash = sha256(Path(__file__))
    source_hashes = {path: sha256(path) for item in inputs for path in item['files']}
    config = {'schema_version': 1, 'implementation_sha256': module_hash, 'genome': genome,
              'annotation_release': annotation_release, 'gtf': str(gtf), 'gtf_sha256': gtf_hash,
              'libraries': {lib: {'path': str(path), 'sha256': sha256(path)} for lib, path in gene_sets.items()},
              'inputs': inputs, 'input_sha256': source_hashes, 'annotation_mode': annotation_mode,
              'background': background, 'min_query_size': min_query_size, 'min_set_size': min_set_size,
              'max_set_size': max_set_size, 'min_overlap': min_overlap, 'fdr': fdr}
    run_signature = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    config_path = output / 'analysis_config.json'
    if output.exists() and any(output.iterdir()):
        if not resume or not config_path.is_file() or json.loads(config_path.read_text()) != config:
            raise ValueError('Output is not empty or its input/resource configuration changed; choose a new output_dir')
    # Import only when the feature is requested, before creating result files.
    import gseapy  # noqa: F401
    output.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(config, indent=2))
    gene_path = output / 'resources' / f'{genome}_{gtf_hash}_genes.csv.gz'
    # Re-parse the selected GTF; never trust a table solely because it exists.
    genes = parse_gtf(gtf, gene_path, genome)
    gtf_header = genes.attrs['gtf_header']
    declared_version = re.fullmatch(r'(?:GENCODE )?v(\d+)', annotation_release)
    header_versions = {v for line in gtf_header for v in re.findall(r'\bversion (\d+)\b', line)}
    if declared_version and header_versions and header_versions != {declared_version[1]}:
        raise ValueError('annotation_release conflicts with the GTF header version')
    for path in gene_sets.values():
        read_gmt(path)
    jobs = []
    for item in inputs:
        for path in item['files']:
            file_id = hashlib.sha256(path.encode()).hexdigest()[:20]
            jobs.append({'source': path, 'output': str(output), 'sample': item['sample'],
                         'caller': item['caller'], 'bin_size': item['bin_size'], 'file_id': file_id,
                         'mode': annotation_mode, 'resume': resume,
                         'signature': {'source_sha256': source_hashes[path], 'gtf_sha256': gtf_hash,
                                       'genome': genome, 'release': annotation_release,
                                       'run_signature': run_signature}})
    if workers == 1:
        init_worker(str(gene_path))
        caches = [process_file(job) for job in jobs]
    else:
        with ProcessPoolExecutor(workers, initializer=init_worker, initargs=(str(gene_path),)) as pool:
            caches = list(pool.map(process_file, jobs))
    records = [json.loads(Path(cache).read_text()) for cache in caches]
    seen_chromosomes = defaultdict(set)
    for record in records:
        key = tuple(record[field] for field in KEY[:3])
        overlap = seen_chromosomes[key] & set(record['chromosomes'])
        if overlap:
            raise ValueError(f'Overlapping chromosome files in {key}: {sorted(overlap)}; do not mix a combined file and its splits')
        seen_chromosomes[key].update(record['chromosomes'])
    args = SimpleNamespace(output=output, gene_path=gene_path, gtf=gtf, genome=genome,
        annotation_release=annotation_release, gtf_header=gtf_header, gene_sets=gene_sets, inputs=inputs,
        run_signature=run_signature, annotation_mode=annotation_mode, background=background,
        min_query_size=min_query_size, min_set_size=min_set_size, max_set_size=max_set_size,
        min_overlap=min_overlap, fdr=fdr)
    return analyze(records, args)
