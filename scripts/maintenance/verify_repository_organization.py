"""Verify organization, input copies and approved result bytes without refitting.

Full local archives are checked when present. --require-full-archive also fails
on intentionally external/ignored historical payloads absent from a clone.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from repository_paths import ROOT, analysis_path, resolve_recorded_path


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--require-full-archive", action="store_true")
    p.add_argument("--output", type=Path, default=ROOT / "miscellaneous/local/repository_cleanup")
    args = p.parse_args()
    checks, missing, mismatches = [], [], []

    def verify(label, records, field, digest="sha256", required=False, resolve=False):
        available = 0
        for row in records:
            path = resolve_recorded_path(row[field]) if resolve else ROOT / row[field]
            if not path.is_file():
                missing.append(dict(group=label, recorded_path=row[field], resolved_path=str(path)))
            else:
                available += 1
                if sha(path) != row[digest]:
                    mismatches.append(dict(group=label, path=str(path), expected=row[digest]))
        absent = len(records) - available
        passed = not any(r["group"] == label for r in mismatches) and not ((required or args.require_full_archive) and absent)
        checks.append(dict(check=label, passed=passed, total=len(records), available=available, missing=absent))

    verify("Repository moved/snapshot bytes", read(ROOT / "docs/provenance/repository_migration_20261116.csv"), "new_path")
    verify("Causal moved bytes", read(ROOT / "docs/provenance/causal_cleanup_20261116.csv"), "new_path")
    approved = read(ROOT / "docs/provenance/current_release_artifacts.csv")
    local_metadata = [r for r in approved if r.get("role") == "local_metadata"]
    verify("Approved scientific/run-record bytes", [r for r in approved if r.get("role") != "local_metadata"], "path")
    current_inputs = read(ROOT / "docs/provenance/current_input_manifest.csv")
    verify("Canonical expanded-marker inputs", current_inputs, "path", required=True)
    historical = analysis_path("current_root") / "QA/historical_sha256.csv"
    if historical.exists():
        verify("Original 104 historical references via relocation resolver", read(historical), "path", resolve=True)
    registry = json.loads((ROOT / "scripts/config/analysis_release.json").read_text())
    for key in ("current_root", "current_fcg_zonation_transport", "current_gdx_zonation_transport", "current_wdsd_zonation_transport",
                "wdsd_hepatocytes_transport", "historical_wdsd_zonation_transport", "wd_cell_types_transport", "wd_zonation_transport",
                "engine_fcg05", "engine_validation07"):
        raw = registry["paths"][key]
        ok = not Path(raw).is_absolute() and (ROOT / raw).exists()
        checks.append(dict(check="Registry: " + key, passed=ok, total=1, available=int(ok), missing=int(not ok)))
    panels = read(ROOT / "docs/provenance/current_panel_manifest.csv")
    panel_paths = [name for row in panels for key in ("files", "source_tables") for name in row[key].split("; ") if name]
    checks.append(dict(check="Portable panel/source paths", passed=all((ROOT / n).is_file() for n in panel_paths),
                       total=len(panel_paths), available=sum((ROOT / n).is_file() for n in panel_paths), missing=sum(not (ROOT / n).is_file() for n in panel_paths)))
    required_doc = ("README.md", "PIPELINE.md", "AGENTS.md", "docs/AGENT_GUIDE.md", "docs/REPOSITORY_ORGANIZATION.md", "docs/CAUSAL_ANALYSIS_MAP.md", "docs/GITHUB_COMMIT_GUIDE.md", "miscellaneous/README.md")
    checks.append(dict(check="Maintained navigation documents", passed=all((ROOT / n).is_file() for n in required_doc), total=len(required_doc), available=sum((ROOT / n).is_file() for n in required_doc), missing=sum(not (ROOT / n).is_file() for n in required_doc)))
    summary = dict(passed=all(c["passed"] for c in checks), checks=checks, byte_mismatches=mismatches,
                   missing_payloads=missing, full_archive_required=args.require_full_archive,
                   excluded_local_metadata=[r["path"] for r in local_metadata],
                   no_models_fitted=True, no_source_or_results_modified=True)
    output = args.output.resolve()
    if not output.is_relative_to(ROOT):
        p.error("Output must be inside the checkout")
    output.mkdir(parents=True, exist_ok=True)
    (output / "organization_verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (output / "organization_checks.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(checks[0])); w.writeheader(); w.writerows(checks)
    print(json.dumps({"passed": summary["passed"], "groups": len(checks), "files_hashed": sum(c["available"] for c in checks[:5]), "byte_mismatches": len(mismatches), "missing_payloads": len(missing)}, indent=2))
    return int(not summary["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
