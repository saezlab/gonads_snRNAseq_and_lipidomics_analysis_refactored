"""Render and check the code-audited workflow contracts without running analyses."""
import argparse
import csv
import glob
import json
from pathlib import Path
import re
import subprocess

from repository_paths import ROOT, release

SPEC = ROOT / "docs/workflow_dependencies.json"
REFERENCE = ROOT / "docs/WORKFLOW_REFERENCE.md"
CATALOG = ROOT / "docs/analysis_catalog.csv"
EDGES = ROOT / "docs/workflow_dependencies.csv"


def tracked_files():
    return set(filter(None, subprocess.check_output(
        ["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")))


def expand_braces(value):
    match = re.search(r"\{([^{}]+,[^{}]+)\}", value)
    if not match:
        return [value]
    return [result for item in match[1].split(",") for result in
            expand_braces(value[:match.start()] + item + value[match.end():])]


def declared_paths(item, registry):
    paths = list(item.get("paths", []))
    paths += [registry["paths"][key] for key in item.get("registry_keys", [])]
    return list(dict.fromkeys(paths))


def availability(item, registry, tracked):
    declarations = declared_paths(item, registry)
    if not declarations:
        return ("External package/resource" if item.get("resources") else "See contract"), []
    files, virtual, missing = set(), False, False
    for declaration in declarations:
        if any(token in declaration for token in ("<", ">", "${", "file.path(")):
            virtual = True
            continue
        for name in expand_braces(declaration):
            p = Path(name)
            if p.is_absolute():
                virtual = True
                continue
            matches = glob.glob(str(ROOT / name), recursive=True)
            if not matches:
                missing = True
            for value in matches:
                candidate = Path(value)
                if candidate.is_dir():
                    files.update(x for x in candidate.rglob("*")
                                 if x.is_file() and x.name != ".DS_Store" and "__pycache__" not in x.parts)
                elif candidate.is_file():
                    files.add(candidate)
                    # A pseudobulk manifest represents every declared payload,
                    # rather than only the small JSON descriptor itself.
                    if candidate.name == "manifest.json":
                        try:
                            manifest = json.loads(candidate.read_text())
                            if manifest.get("format") == "pseudobulk_csv":
                                files.update(candidate.parent / f["file"] for f in manifest["files"])
                        except (ValueError, KeyError):
                            pass
    local = [p for p in files if p.is_file()]
    missing = missing or len(local) != len(files)
    shipped = [p for p in local if p.relative_to(ROOT).as_posix() in tracked]
    if files and len(shipped) == len(files):
        state = "Shipped in Git"
    elif shipped:
        state = "Mixed: some files are local-only"
    elif local:
        state = "Local-only; restore or regenerate"
    elif virtual:
        state = "External / runtime-supplied"
    else:
        state = "Absent locally; restore or generate"
    if missing and (local or virtual):
        state += "; incomplete: declared paths absent"
    if virtual and local:
        state += "; also runtime/external paths"
    return state, [p.relative_to(ROOT).as_posix() for p in sorted(local)]


def cell(value):
    return str(value or "").replace("|", "\\|").replace("\n", "<br>")


def source_reference(value):
    match = re.fullmatch(r"(.+?):(\d+)(?:[-–](\d+))?", value)
    if not match:
        return value, ""
    anchor = "#L" + match[2] + ("-L" + match[3] if match[3] else "")
    return match[1], anchor


def code_paths(values):
    return "<br>".join(f"`{cell(p)}`" for p in values) or "See mode notes"


def display_paths(item, registry):
    values = list(item.get("paths", []))
    values += [f"registry:{key} → {registry['paths'][key]}" for key in item.get("registry_keys", [])]
    values += [f"resource:{name}" for name in item.get("resources", [])]
    return code_paths(values)


def csv_text(rows, fields):
    import io
    stream = io.StringIO(newline="")
    w = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    w.writeheader(); w.writerows(rows)
    return stream.getvalue()


def render(spec, registry, tracked):
    workflows = spec["workflows"]
    lines = ["# Workflow execution and data reference", "",
             f"Code audit: **{spec['audited_on']}**, working tree based on commit `{spec['source_commit']}`.", "",
             "Generated from [workflow_dependencies.json](workflow_dependencies.json). "
             "Use [PIPELINE.md](../PIPELINE.md) for the short execution map and the linked per-workflow plans for methods. "
             "This reference distinguishes reading an existing result from producing it, and orchestration from data flow.", "",
             "**Availability is a checked-out-file inventory, not a full execution guarantee.** "
             "A shipped derived result can be consumed without rerunning its producer. Local model caches and "
             "private/external inputs must be regenerated or restored as stated. Validation flags have the exact scope "
             "shown below; none is a substitute for a scientific reproduction run.", "",
             "**Numbering is not execution order.** Workflow17 orchestrates the current zonation/diet branch; "
             "compatibility wrappers and external/reference scripts are not extra sequential main-pipeline steps.", "",
             "## Workflow index", "",
             "| ID | Role | Entry point | Detailed plan |", "| --- | --- | --- | --- |"]
    catalog, edges = [], []
    for w in workflows:
        entries = "<br>".join(f"[{Path(p).name}](../{p})" for p in w["entrypoints"])
        plan = w.get("plan", "")
        plan_link = f"[{Path(plan).name}]({str(Path(plan).relative_to('docs'))})" if plan else "See notes"
        lines.append(f"| [{w['id']}](#workflow-{w['id'].lower()}) | {cell(w['status'])} | {entries} | {plan_link} |")
    for w in workflows:
        lines += ["", f"## Workflow {w['id']}", "", f"**{w['name']} — {w['status']}**", "", w["purpose"], "",
                  "**Entry:** " + "; ".join(f"`{p}`" for p in w["entrypoints"]), "",
                  "**Validation:** " + w["validation"], "", "**Writes / overwrite behavior:** " + w["overwrite"], ""]
        if w.get("calls"):
            lines += ["**Invokes other entries/phases:** " + ", ".join(f"`{p}`" for p in w["calls"]) + ". These are control-flow calls, not additional cohorts.", ""]
        lines += ["| Required input / role | File or registry key | Producer / mode | Availability and action |",
                  "| --- | --- | --- | --- |"]
        for item in w.get("inputs", []):
            state, local = availability(item, registry, tracked)
            note = item.get("availability_note", "")
            lines.append(f"| {cell(item['name'])}<br>{cell(item.get('role', ''))} | {display_paths(item, registry)} | "
                         f"{cell(item.get('producer', 'curated'))}<br>{cell(item.get('mode', 'all'))} | {cell(state)}<br>{cell(note)} |")
            edges.append(dict(consumer=w["id"], input=item["name"], producer=item.get("producer", "curated"),
                              mode=item.get("mode", "all"), role=item.get("role", ""),
                              paths="; ".join(item.get("paths", [])), registry_keys="; ".join(item.get("registry_keys", [])),
                              resources="; ".join(item.get("resources", [])), availability=state, action=note))
        lines += ["", "| Output | Location | Mode |", "| --- | --- | --- |"]
        for item in w.get("outputs", []):
            lines.append(f"| {cell(item['name'])} | {display_paths(item, registry)} | {cell(item.get('mode', 'all'))} |")
        if w.get("commands"):
            lines += ["", "**Command examples** (not executed by this documentation audit):", ""]
            for command in w["commands"]:
                lines += [command["effect"], "", "```sh", command["command"], "```", ""]
        if w.get("notes"):
            lines += ["**Important distinctions:**", ""] + ["- " + note for note in w["notes"]] + [""]
        if w.get("evidence"):
            lines += ["**Code checked:** " + ", ".join(f"[{p}](../{source_reference(p)[0]}{source_reference(p)[1]})" for p in w["evidence"]) + ".", ""]
        catalog.append(dict(id=w["id"], name=w["name"], status=w["status"], entry="; ".join(w["entrypoints"]),
                            inputs="; ".join(i["name"] for i in w.get("inputs", [])),
                            outputs="; ".join(" / ".join(declared_paths(o, registry)) for o in w.get("outputs", [])),
                            dependency_note="; ".join(w.get("notes", [])), plan=w.get("plan", "")))
    return {REFERENCE: "\n".join(lines).rstrip() + "\n",
            CATALOG: csv_text(catalog, ["id", "name", "status", "entry", "inputs", "outputs", "dependency_note", "plan"]),
            EDGES: csv_text(edges, ["consumer", "input", "producer", "mode", "role", "paths", "registry_keys", "resources", "availability", "action"])}


def check_spec(spec, registry, tracked):
    workflows = spec["workflows"]
    workflow_ids = {w["id"] for w in workflows}
    assert len(workflow_ids) == len(workflows), "Duplicate workflow ID"
    producer_ids = workflow_ids | {"curated", "external", "frozen"}
    covered = set()
    for w in workflows:
        for key in ("id", "name", "entrypoints", "purpose", "status", "validation", "overwrite"):
            assert w.get(key), (w.get("id"), "Missing contract", key)
        for call in w.get("calls", []):
            assert call in workflow_ids, (w["id"], "Unknown called workflow ID", call)
        for p in w["entrypoints"] + w.get("evidence", []) + ([w["plan"]] if w.get("plan") else []):
            assert (ROOT / source_reference(p)[0]).exists(), (w["id"], "Missing source/document", p)
        covered.update(w["entrypoints"])
        for item in w.get("inputs", []) + w.get("outputs", []):
            for key in item.get("registry_keys", []):
                assert key in registry["paths"], (w["id"], "Unknown registry key", key)
        for item in w.get("inputs", []):
            assert item.get("paths") or item.get("registry_keys") or item.get("resources"), (w["id"], "Input has no file or resource location", item)
            assert item.get("producer") and item.get("mode") and item.get("role"), (w["id"], "Incomplete input contract", item)
            assert item["producer"] in producer_ids, (w["id"], "Unknown producer ID", item["producer"])
            state, _ = availability(item, registry, tracked)
            if state != "Shipped in Git":
                assert item.get("availability_note"), (w["id"], "Missing restoration/generation note", item["name"], state)
    entries = {p.relative_to(ROOT).as_posix() for p in (ROOT / "scripts/workflows").rglob("*")
               if p.is_file() and p.suffix.lower() in (".r", ".py") and "__pycache__" not in p.parts}
    assert not (entries - covered), "Undocumented workflow entries: " + ", ".join(sorted(entries - covered))
    return {"workflows": len(workflows), "workflow_entry_files": len(entries),
            "input_groups": sum(len(w.get("inputs", [])) for w in workflows),
            "output_groups": sum(len(w.get("outputs", [])) for w in workflows),
            "source_and_registry_contracts_checked": True}


def build(check=False):
    spec = json.loads(SPEC.read_text()); registry = release(); tracked = tracked_files()
    summary = check_spec(spec, registry, tracked)
    outputs = render(spec, registry, tracked)
    if check:
        stale = [p.relative_to(ROOT).as_posix() for p, text in outputs.items()
                 if not p.exists() or p.read_text() != text]
        assert not stale, "Stale generated documentation: " + ", ".join(stale)
    else:
        for p, text in outputs.items():
            p.write_text(text)
    print(json.dumps(dict(summary, generated_documents=len(outputs), check_only=check, scientific_models_run=False), indent=2))
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check", action="store_true", help="Validate coverage/paths and detect stale generated docs without writing")
    build(p.parse_args().check)
