"""Rebuild the complete mouse PKN from the frozen source snapshots.

Requires pygoslin 2.2.5 plus its pinned dependencies, and R cosmosR/jsonlite.
Example: uv run --with pygoslin==2.2.5 --with scipy==1.13.1
         --with numpy==2.0.2 python scripts/pkn/run_build.py
This intentionally does not refresh remote resources.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    scripts = ['build_mouse_gem.py', 'build_rhea_lipid_reactions.py',
               'build_lipid_regulation.py', 'map_lipids.py', 'build_lipid_rules.py',
               'assemble_mouse_pkn.py', 'audit_mouse_pkn.py']
    for script in scripts:
        subprocess.run([sys.executable, str(ROOT / 'scripts/pkn' / script)], cwd=ROOT, check=True)
    subprocess.run(['Rscript', 'scripts/pkn/validate_export.R'], cwd=ROOT, check=True)
    subprocess.run([sys.executable, 'scripts/pkn/write_build_manifest.py'], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
