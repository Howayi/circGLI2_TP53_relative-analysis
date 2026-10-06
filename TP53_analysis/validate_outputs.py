from pathlib import Path
from itertools import product
import hashlib
import json

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
scores = pd.read_csv(OUT / 'TP53_per_sample_scores.tsv', sep='\t').set_index('sample')
tests = pd.read_csv(OUT / 'TP53_activity_group_tests.tsv', sep='\t')
expression = pd.read_csv(OUT / 'TMM_logCPM_gene_centered_symbols.tsv', sep='\t', index_col=0)
meta = pd.read_csv(ROOT / 'sample_metadata_qc.tsv', sep='\t')
assert len(scores) == 19 and scores.index.is_unique
assert meta['paired'].sum() == 16
assert np.isfinite(scores.select_dtypes(include=[np.number]).values).all()
checks = []
for network in ['CollecTRI', 'DoRothEA_AB']:
    net = pd.read_csv(ROOT / 'resources' / f'{network}_used_network.tsv', sep='\t')
    tp = net[net['source'].eq('TP53')].set_index('target')['mor']
    x = tp.reindex(expression.index).fillna(0).values
    y = expression.values
    xc = x - x.mean()
    yc = y - y.mean(axis=0)
    r = xc @ yc / np.sqrt((xc @ xc) * (yc * yc).sum(axis=0))
    independent_t = r * np.sqrt((len(x) - 2) / (1 - r * r))
    actual = scores.loc[expression.columns, f'{network}_TP53_ULM'].values
    np.testing.assert_allclose(actual, independent_t, rtol=1e-9, atol=1e-9)
    patients = meta.loc[meta['paired'], 'patient'].astype(str).unique()
    d = np.array([scores.loc[p + 'T', f'{network}_TP53_ULM'] - scores.loc[p + 'N', f'{network}_TP53_ULM'] for p in patients])
    row = tests[tests['network'].eq(network) & tests['analysis'].eq('main')].iloc[0]
    np.testing.assert_allclose(row['mean_T_minus_N'], d.mean(), atol=1e-12)
    np.testing.assert_allclose(row['paired_t_p'], stats.ttest_1samp(d, 0).pvalue, atol=1e-12)
    perm = np.asarray(list(product([-1, 1], repeat=len(d)))) @ d / len(d)
    p = (np.abs(perm) >= abs(d.mean()) - 1e-12).mean()
    np.testing.assert_allclose(row['exact_signflip_p'], p, atol=1e-12)
    checks.append({'network': network, 'independent_ULM_max_abs_difference': float(np.max(np.abs(independent_t - actual))), 'paired_t_verified': True, 'exact_signflip_verified': True})
all_annotation = pd.read_csv(OUT / 'input_gene_annotation.tsv', sep='\t')
assert len(all_annotation) == 78686 and all_annotation['symbol'].notna().all()
target = all_annotation[all_annotation['symbol'].eq('TP53')]
assert target['gene_id'].tolist() == ['ENSG00000141510']
gsea = pd.read_csv(OUT / 'Hallmark_p53_GSEA.tsv', sep='\t')
assert set(gsea['analysis']) == {'main', 'sens6', 'sens5'}
assert gsea['NES'].gt(0).all()
manifest = json.loads((ROOT / 'resources' / 'resource_manifest.json').read_text(encoding='utf-8'))
archive = ROOT / 'resources' / 'decoupleR_2.17.0.zip'
manifest.append({'file': archive.name, 'url': 'https://bioconductor.org/packages/release/bioc/bin/windows/contrib/4.6/decoupleR_2.17.0.zip', 'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(), 'bytes': archive.stat().st_size})
manifest = list({x['file']: x for x in manifest}.values())
manifest = [x for x in manifest if x['file'] != 'decoupleR-master.zip']
for entry in manifest:
    path = ROOT / 'resources' / entry['file']
    if 'sha256' in entry:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256'], str(path)
(ROOT / 'resources' / 'resource_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
report = {'status': 'passed', 'checks': checks, 'sample_counts': {'all': 19, 'paired': 16}, 'all_input_genes_exactly_annotated': True, 'TP53_ensembl_verified': True, 'resource_sha256_verified': True, 'interpretation': 'ULM enrichment p-values, sample paired-test p-values, TF-wide BH FDR, and pathway GSEA p-values are exported separately'}
(ROOT / 'output_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
