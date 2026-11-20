#!/usr/bin/env python3
"""Replay the manuscript workflows into a parallel results tree."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from repository_paths import release, resolve_recorded_path


REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "results"
RUNS = RESULTS / "runs"
SECTIONS = (
    "01_gonadal_sex_and_hepatic_lipid_homeostasis",
    "02_gonadal_sex_and_hepatocyte_gene_expression",
    "03_gonadal_sex_biased_regulatory_circuitry",
    "04_pericentral_and_periportal_hepatocyte_bias",
    "05_sex_polarised_zonation_gene_lipid_interactions",
    "06_western_diet_and_gene_lipid_circuitry",
)
WORKFLOWS = {
    "01": "01_data_preparation_qc.R",
    "02": "02_fcg_pseudobulk_analysis.R",
    "03": "03_fcg_lipidomic_analysis.R",
    "04": "04_mofa_integration.R",
    "06": "06_structural_analysis.R",
    "07": "07_validation_experiments.R",
    "08": "08_variance_analysis.R",
    "11": "11_cyp2c38_in_vivo_validation.R",
    "17": "17_expanded_marker_zonation_update.R",
}
RNA_FILES = ("meta_data_df.csv", *(f"scaled_RNA_{i}.csv" for i in range(1, 4)))
GDX_FILES = ("meta_data_df_gdxhr.csv", *(f"scaled_RNA_gdxhr_{i}.csv" for i in range(1, 4)))


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def link(source: Path, target: Path) -> None:
    if target.is_symlink() and target.resolve() == source.resolve():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    target.symlink_to(os.path.relpath(source, target.parent), target_is_directory=source.is_dir())


def json_write(path: Path, data: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


class Replay:
    def __init__(self, output: Path, scratch: Path):
        self.output = output
        self.scratch = scratch
        self.logs = scratch / "logs"
        self.prep = scratch / "prepare_01"
        self.cwd = scratch / "workspace"
        self.zonation = self.cwd / "results/expanded_marker_replay"
        self.paths = {
            "results/qc": output / SECTIONS[1] / "qc",
            "results/pseudobulk": output / SECTIONS[1] / "pseudobulk",
            "results/lipidomic": output / SECTIONS[0] / "lipidomic",
            "results/mofa": output / SECTIONS[2] / "mofa",
            "results/validation": output / SECTIONS[2] / "in_vitro_validation",
        }
        relocations = {item["from"]: item["to"] for item in release(REPO)["relocations"]}
        self.routes = {relocations[relative]: target for relative, target in self.paths.items()}
        self.zonal_paths = {
            "fcg05": output / SECTIONS[4] / "fcg05",
            "fcg07": output / SECTIONS[4] / "fcg07",
            "gdx_zonation": output / SECTIONS[3] / "gdx_zonation",
            "tables": output / SECTIONS[4] / "supplementary_tables",
        }
        for branch in ("cf453_models", "cf453_direct_diet", "cf453_rewiring", "cf453_layouts"):
            self.zonal_paths[branch] = output / SECTIONS[5] / branch
        self.models = resolve_recorded_path("results/mofa/RNA", REPO)
        self.env = os.environ.copy()
        libraries = (scratch / "R-library", REPO / ".local/R-library")
        configured = [str(path) for path in libraries if path.is_dir()]
        configured += self.env.get("R_LIBS_USER", "").split(os.pathsep)
        if any(configured):
            self.env["R_LIBS_USER"] = os.pathsep.join(filter(None, configured))
        self.env.setdefault("VECLIB_MAXIMUM_THREADS", "1")
        self.env.setdefault("OPENBLAS_NUM_THREADS", "1")
        self.logs.mkdir(parents=True, exist_ok=True)
        for section in SECTIONS:
            (output / section).mkdir(parents=True, exist_ok=True)

    def setup(self, workflow: str) -> Path:
        if workflow == "01":
            for name in ("fcg_pseudobulks.rds", "gdx_pseudobulks.rds", "pseudobulk_csv"):
                link(REPO / "data/RNA" / name, self.prep / "data/RNA" / name)
            (self.prep / "data/RNA/zonation/current/fcg").mkdir(parents=True, exist_ok=True)
            link(REPO / "data/Suppl_Table_1_OverviewData.xlsx",
                 self.prep / "data/Suppl_Table_1_OverviewData.xlsx")
            self.paths["results/qc"].mkdir(parents=True, exist_ok=True)
            qc_route = next(route for route, target in self.routes.items() if target == self.paths["results/qc"])
            link(self.paths["results/qc"], self.prep / qc_route)
            cwd = self.prep
        else:
            cwd = self.cwd
            for source in (REPO / "data").iterdir():
                if source.name != "RNA":
                    link(source, cwd / "data" / source.name)
            # The split CSVs approved with the release are the scientific input
            # contract. Workflow 01 independently re-exports them for comparison.
            link(REPO / "data/RNA", cwd / "data/RNA")
            for relative, destination in self.routes.items():
                destination.mkdir(parents=True, exist_ok=True)
                link(destination, cwd / relative)
            # Cache metadata stores canonical paths after resolving the output
            # links. Make those paths resolvable from this working directory too.
            link(self.output, cwd / self.output.relative_to(REPO))
            xenium = resolve_recorded_path("results/xenium", REPO)
            link(xenium, cwd / xenium.relative_to(REPO))
            for relative, destination in self.zonal_paths.items():
                destination.mkdir(parents=True, exist_ok=True)
                link(destination, self.zonation / relative)
        for name in ("scripts", "support", "docs", "miscellaneous", "tests"):
            link(REPO / name, cwd / name)
        if (REPO / ".local/R-library").is_dir():
            link(REPO / ".local/R-library", cwd / ".local/R-library")
        return cwd

    def command(self, workflow: str) -> list[str]:
        command = ["Rscript", "--vanilla", f"scripts/workflows/{WORKFLOWS[workflow]}"]
        if workflow == "04":
            command += ["--reuse-models-only", f"--model-dir={self.models}"]
        if workflow == "07":
            command += ["--global-only", f"--output={self.paths['results/validation']}"]
        if workflow == "11":
            command += [f"--output={self.output / SECTIONS[2] / 'cyp2c38_in_vivo'}"]
        if workflow == "17":
            command += [f"--output={self.zonation}", f"--python={sys.executable}"]
        return command

    def check_inputs(self, workflow: str) -> None:
        lipid = self.paths["results/lipidomic"]
        mofa = self.paths["results/mofa"]
        kd = self.paths["results/validation"]
        requirements = {
            "06": [lipid / "ttop_full_list.RData"],
            "07": [mofa / "networks/hepatocyte_gene_lipid_network_full_SIF.csv",
                   lipid / "structural_model/lipids_with_structure.csv"],
            "11": [mofa / "networks/hepatocyte_gene_lipid_network_full_SIF.csv",
                   self.paths["results/pseudobulk"] / "PCA_sample_coord_RNA_HEPA.csv"],
            "17": [lipid / "ttop_full_list.RData",
                   lipid / "enrichment/lipid_class_enrichment_scores.RData",
                   lipid / "structural_model/lipids_with_structure.csv",
                   kd / "knockdown_t_table.csv", kd / "knockdown_ttop_list.RData"],
        }
        for path in requirements.get(workflow, ()):
            if not path.is_file():
                raise FileNotFoundError(f"Workflow {workflow} requires replay output {path}")

    def check_fresh(self, workflow: str) -> None:
        lipid = self.paths["results/lipidomic"]
        destinations = {
            "01": [self.paths["results/qc"], self.prep / "data/RNA/scaled_RNA_1.csv"],
            "02": [self.paths["results/pseudobulk"]],
            "03": [lipid / "pca", lipid / "ttop_full_list.RData"],
            "04": [self.paths["results/mofa"]],
            "06": [lipid / "structural_model"],
            "07": [self.paths["results/validation"]],
            "08": [lipid / "variance_analysis"],
            "11": [self.output / SECTIONS[2] / "cyp2c38_in_vivo"],
            "17": list(self.zonal_paths.values()),
        }
        for path in destinations[workflow]:
            if path.is_file() or (path.is_dir() and any(path.iterdir())):
                raise FileExistsError(f"Preserve prior replay at {path}. Use a fresh replay destination.")

    def run(self, workflow: str) -> None:
        self.check_fresh(workflow)
        self.check_inputs(workflow)
        cwd = self.setup(workflow)
        command = self.command(workflow)
        record_name = f"workflow_{workflow}"
        attempt = 1
        while (self.logs / f"{record_name}.log").exists():
            attempt += 1
            record_name = f"workflow_{workflow}_attempt_{attempt}"
        log = self.logs / f"{record_name}.log"
        status = log.with_suffix(".json")
        target_paths = {alias: str(path) for alias, path in self.routes.items()}
        target_paths.update({f"expanded_marker_replay/{key}": str(value)
                             for key, value in self.zonal_paths.items()})
        record = {
            "workflow": workflow, "attempt": attempt, "cwd": str(cwd), "command": command,
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "output_routes": target_paths, "log": str(log),
            "entry_sha256": sha256(cwd / command[2]),
        }
        if workflow == "04":
            models = [self.models / f"mofa_res_{i}factor.hdf5" for i in range(5, 16)]
            record["saved_mofa_models_sha256"] = {str(p): sha256(p) for p in models}
        json_write(status, record)
        print(f"Workflow {workflow} → {log}", flush=True)
        with log.open("x") as stream:
            stream.write(f"cwd: {cwd}\n{shlex.join(command)}\n\n")
            stream.flush()
            result = subprocess.run(command, cwd=cwd, env=self.env, stdout=stream,
                                    stderr=subprocess.STDOUT, check=False)
        record["finished_utc"] = datetime.now(timezone.utc).isoformat()
        record["exit_status"] = result.returncode
        if workflow == "04":
            record["mofa_model_bytes_unchanged"] = all(
                sha256(Path(p)) == digest for p, digest in record["saved_mofa_models_sha256"].items()
            )
        json_write(status, record)
        result.check_returncode()
        if workflow == "04" and not record["mofa_model_bytes_unchanged"]:
            raise RuntimeError("Saved MOFA models changed during replay")
        if workflow == "01":
            self.record_csv_hashes()

    def record_csv_hashes(self) -> None:
        names = (*RNA_FILES, *GDX_FILES)
        with (self.logs / "workflow_01_csv_hashes.csv").open("x", newline="") as stream:
            csv_writer = csv.writer(stream)
            csv_writer.writerow(("file", "approved_sha256", "replay_sha256", "identical_bytes"))
            for name in names:
                approved = REPO / "data/RNA" / name
                replay = self.prep / "data/RNA" / name
                if approved.is_file() and replay.is_file():
                    before, after = sha256(approved), sha256(replay)
                    csv_writer.writerow((name, before, after, before == after))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULTS / "manuscript")
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--workflow", nargs="+", choices=WORKFLOWS, default=list(WORKFLOWS))
    args = parser.parse_args()
    output, scratch = args.output.resolve(), args.scratch.resolve()
    if RUNS not in scratch.parents:
        parser.error("--scratch must be inside results/runs")
    if output != RESULTS / "manuscript" and RUNS not in output.parents:
        parser.error("--output must be results/manuscript or inside results/runs")
    if output == scratch or output in scratch.parents or scratch in output.parents:
        parser.error("--output and --scratch must be separate directories")
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".replay.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error(f"A replay is already writing to {output}")
        replay = Replay(output, scratch)
        for workflow in dict.fromkeys(args.workflow):
            replay.run(workflow)


if __name__ == "__main__":
    main()
