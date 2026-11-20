"""Hash every build source and deliverable, including code and local feature inputs."""
import datetime
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/causal_modelling/cosmos_mouse_pkn'


def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            sha.update(chunk)
    return dict(path=str(path.relative_to(ROOT)), bytes=path.stat().st_size, sha256=sha.hexdigest())


def main():
    raw = sorted((ROOT / 'support/cosmos_mouse/raw').glob('*'))
    raw += [ROOT / 'results/causal_modelling/cosmos_pkn_audit/raw' / name for name in
            ['om-fresh-grn.tsv.gz', 'om-fresh-ppi.tsv.gz', 'om-mouse-pkn.json.gz',
             'swisslipids-enzymes.tsv.gz']]
    inputs = [ROOT / p for p in ['results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/gonade_type_ttop.csv',
              'results/current/03_gonadal_sex_biased_regulatory_circuitry/mofa/networks/gene_lipid_correlations_factor1.csv']]
    code = sorted((ROOT / 'scripts/pkn').glob('*.py')) + sorted((ROOT / 'scripts/pkn').glob('*.R'))
    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name != 'build_manifest.json')
    result = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  task='Mouse COSMOS prior construction; no biological input mapping or scoring',
                  raw_sources=[digest(p) for p in raw if p.is_file()],
                  local_feature_inputs=[digest(p) for p in inputs],
                  build_code=[digest(p) for p in code], outputs=[digest(p) for p in outputs],
                  documentation=[digest(ROOT / p) for p in ['docs/cosmos_mouse_pkn_build.md',
                                  'support/cosmos_mouse/README.md']],
                  upstream_version_manifests=[str(p.relative_to(ROOT)) for p in outputs
                                             if 'manifest' in p.name or p.name=='mapping_metadata.json'],
                  cached_protein_export_urls_manifest='results/causal_modelling/cosmos_pkn_audit/manifest.json')
    (OUT / 'build_manifest.json').write_text(json.dumps(result, indent=2) + '\n')
    print(f'Hashed {len(result["raw_sources"])} source snapshots, {len(code)} scripts and {len(outputs)} outputs.')


if __name__ == '__main__':
    main()
