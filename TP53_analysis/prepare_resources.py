from pathlib import Path
import gzip
import hashlib
import io
import os
import json
import re
import shutil
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
RES = ROOT / 'resources'
RES.mkdir(parents=True, exist_ok=True)
INPUT = Path(os.environ.get('TP53_INPUT_DIR', str(ROOT.parent / 'input_data' / 'counts')))


def download(url, filename):
    path = RES / filename
    if not path.exists():
        response = requests.get(url, timeout=(20, 180))
        response.raise_for_status()
        path.write_bytes(response.content)
    return path


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    counts = pd.read_csv(INPUT / 'gene_counts.matrix.tsv', sep='\t', index_col=0)
    assert counts.index.is_unique and counts.columns.is_unique
    assert counts.apply(pd.to_numeric).ge(0).all().all()
    assert np.equal(counts.values, np.floor(counts.values)).all()
    full = pd.read_csv(INPUT / 'gene_counts.featureCounts.txt', sep='\t', comment='#', usecols=lambda col: col not in ['Chr', 'Start', 'End', 'Strand'])
    full = full.set_index('Geneid')
    extracted = full.iloc[:, 1:].copy()
    extracted.columns = [col.rsplit('/', 1)[-1].replace('.sorted.bam', '') for col in extracted.columns]
    assert counts.equals(extracted), 'The simplified matrix differs from the featureCounts source'
    summary = pd.read_csv(INPUT / 'gene_counts.featureCounts.txt.summary', sep='\t', index_col=0)
    summary.columns = [col.rsplit('/', 1)[-1].replace('.sorted.bam', '') for col in summary.columns]
    assert list(summary.columns) == list(counts.columns)
    assert np.array_equal(counts.sum().values, summary.loc['Assigned'].values)
    samples = pd.DataFrame({'sample': counts.columns})
    samples['patient'] = samples['sample'].str[:-1]
    samples['condition'] = samples['sample'].str[-1]
    paired = samples.groupby('patient')['condition'].nunique().eq(2)
    samples['paired'] = samples['patient'].map(paired)
    samples['assigned_fragments'] = summary.loc['Assigned'].values
    samples['summary_total'] = summary.sum().values
    for label, row in [('assigned', 'Assigned'), ('multimapping', 'Unassigned_MultiMapping'), ('singleton', 'Unassigned_Singleton'), ('no_feature', 'Unassigned_NoFeatures')]:
        samples[label + '_fraction'] = summary.loc[row].values / summary.sum().values
    samples['detected_genes'] = counts.gt(0).sum().values
    samples.to_csv(ROOT / 'sample_metadata_qc.tsv', sep='\t', index=False)
    full[['Length']].to_csv(RES / 'featurecounts_gene_lengths.tsv', sep='\t')
    stripped = counts.index.str.replace(r'\.[0-9]+$', '', regex=True)
    assert stripped.is_unique, 'Version stripping produces duplicate Ensembl IDs'
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(), 'input_shape': list(counts.shape), 'input_mode': 'raw_paired_fragment_counts', 'matrix_matches_featurecounts': True, 'column_sums_match_Assigned': True, 'files': {str(INPUT / name): {'sha256': sha256(INPUT / name)} for name in ['gene_counts.matrix.tsv', 'gene_counts.featureCounts.txt', 'gene_counts.featureCounts.txt.summary']}, 'paired_patients': samples.loc[samples['paired'], 'patient'].drop_duplicates().tolist()}
    script = ROOT.parent / 'upstream' / 'rnaseq_upstream922.sh'
    manifest['files'][str(script)] = {'sha256': sha256(script)}
    shutil.copyfile(script, RES / script.name)
    (ROOT / 'input_validation.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(samples.to_string(index=False), flush=True)
    print('Input validation passed', flush=True)

    urls = {
        'gencode.v48.annotation.gtf.gz': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz',
        'collectri_regulons.csv': 'https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1',
        'dorothea_hs.rda': 'https://raw.githubusercontent.com/saezlab/dorothea/master/data/dorothea_hs.rda',
        'decoupleR_2.17.0.zip': 'https://bioconductor.org/packages/3.23/bioc/bin/windows/contrib/4.6/decoupleR_2.17.0.zip',
        'hallmark_p53.gmt': 'https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=HALLMARK_P53_PATHWAY&fileType=gmt',
    }
    statuses = []
    for name, url in urls.items():
        try:
            path = download(url, name)
            statuses.append({'file': name, 'url': url, 'sha256': sha256(path), 'bytes': path.stat().st_size})
            print('Downloaded', name, path.stat().st_size, flush=True)
        except Exception as exc:
            statuses.append({'file': name, 'url': url, 'error': str(exc)})
            print('DOWNLOAD FAILED', name, str(exc), flush=True)
    (RES / 'resource_manifest.json').write_text(json.dumps(statuses, indent=2), encoding='utf-8')
    path = RES / 'gencode.v48.annotation.gtf.gz'
    if path.exists():
        records = []
        with gzip.open(path, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if fields[2] != 'gene':
                    continue
                attr = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
                records.append({'gene_id_version': attr['gene_id'], 'gene_id': attr['gene_id'].split('.')[0], 'symbol': attr.get('gene_name', ''), 'gene_type': attr.get('gene_type', ''), 'chromosome': fields[0]})
        pd.DataFrame(records).to_csv(RES / 'gencode_v48_gene_mapping.tsv', sep='\t', index=False)
        print('Exact GENCODE v48 gene mapping:', len(records), flush=True)


if __name__ == '__main__':
    main()
