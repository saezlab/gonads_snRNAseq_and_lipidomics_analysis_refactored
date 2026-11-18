"""Index public CSV pseudobulk inputs; serialized original objects stay local."""
import csv
import hashlib
from pathlib import Path
from repository_paths import ROOT, analysis_path

ROLES = ("current_fcg_zonation", "current_gdx_zonation", "current_wdsd_zonation",
         "wdsd_hepatocytes", "historical_wdsd_zonation", "wd_cell_types", "wd_zonation")
LIMIT = 50_000_000


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main():
    rows = []
    for role in ROLES:
        manifest = analysis_path(role + "_transport")
        assert manifest.is_file(), manifest
        files = sorted(p for p in manifest.parent.iterdir() if p.is_file() and p.suffix.lower() in (".csv", ".json"))
        assert any(p.name == "manifest.csv" for p in files), role
        for path in files:
            if path.suffix.lower() == ".csv":
                assert path.stat().st_size < LIMIT, path
            rows.append(dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path), bytes=path.stat().st_size,
                             transport_role=role, format=path.suffix[1:], source_manifest=manifest.relative_to(ROOT).as_posix()))
    # Existing small preprocessor sources and legacy wide CSV interfaces remain
    # public. The old ledger retains the exact local RDS identities separately.
    old = ROOT / "docs/provenance/current_input_manifest_20261116.csv"
    with old.open(encoding="utf-8-sig", newline="") as f:
        originals = list(csv.DictReader(f))
    for row in originals:
        path = ROOT / row["path"]
        if path.suffix.lower() not in (".csv", ".r"):
            continue
        assert path.is_file() and sha(path) == row["sha256"], path
        assert path.stat().st_size < LIMIT, path
        rows.append(dict(path=row["path"], sha256=row["sha256"], bytes=path.stat().st_size,
                         transport_role="preprocessing_or_compatibility", format=path.suffix[1:], source_manifest=""))
    assert len({r["path"] for r in rows}) == len(rows)
    target = ROOT / "docs/provenance/current_input_manifest.csv"
    with target.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"Public input index: {len(rows)} CSV/JSON/source files across {len(ROLES)} pseudobulk transports; no serialized RNA objects")


if __name__ == "__main__":
    main()
