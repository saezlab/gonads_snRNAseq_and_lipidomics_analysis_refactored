"""Resolve recorded repository paths without rewriting frozen provenance records."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def release(root=ROOT):
    return json.loads((Path(root) / "scripts/config/analysis_release.json").read_text())


def analysis_path(key, root=ROOT):
    return Path(root) / release(root)["paths"][key]


def resolve_recorded_path(recorded, root=ROOT):
    root = Path(root).resolve()
    path = Path(recorded)
    if path.is_absolute():
        try:
            relative = path.relative_to(root).as_posix()
        except ValueError:
            # Frozen records were created in this named checkout on another host.
            marker = "/gonads_snRNAseq_and_lipidomics_analysis_refactored/"
            if marker not in path.as_posix():
                return path
            relative = path.as_posix().split(marker, 1)[1]
    else:
        relative = path.as_posix().removeprefix("./")
    for move in sorted(release(root).get("relocations", []), key=lambda x: -len(x["from"])):
        before = move["from"].rstrip("/")
        if relative == before or relative.startswith(before + "/"):
            relative = move["to"].rstrip("/") + relative[len(before):]
            break
    return root / relative


def resolve_manifest_path(recorded, base, root=ROOT):
    """Support both portable bundle-relative and original production paths."""
    path = Path(recorded)
    if not path.is_absolute() and (Path(base) / path).exists():
        return Path(base) / path
    return resolve_recorded_path(recorded, root)
