#!/usr/bin/env python3
"""Inventory manuscript claims without editing the DOCX or fitting any model.

Run with the bundled Python runtime. Paragraph locators are one-based direct
children of w:body, including blank paragraphs; they are not page numbers.
Optional report/manuscript_numeric_updates.json adds sourced replacement values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from zipfile import ZipFile
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "maintenance"))
from repository_paths import analysis_path

from lxml import etree
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
    "w14": "http://schemas.microsoft.com/office/word/2010/wordml",
}
W = "{" + NS["w"] + "}"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def text_of(paragraph) -> str:
    """Current displayed text, excluding deleted text and nested paragraphs."""
    pieces = []
    for element in paragraph.iter():
        if element.tag not in {W + "t", W + "tab", W + "br"}:
            continue
        nearest = next((p for p in element.iterancestors() if p.tag == W + "p"), None)
        if nearest is not paragraph or any(p.tag == W + "del" for p in element.iterancestors()):
            continue
        pieces.append(element.text or "" if element.tag == W + "t" else "\t" if element.tag == W + "tab" else "\n")
    return "".join(pieces)


def extract_docx(path: Path, report: Path) -> tuple[list[dict], dict[int, str], list[dict]]:
    locations, images = [], []
    with ZipFile(path) as archive:
        document = etree.fromstring(archive.read("word/document.xml"))
        body = document.find("w:body", NS)
        body_paragraphs = body.findall("w:p", NS)
        body_ids = {id(p): i for i, p in enumerate(body_paragraphs, 1)}
        body_text = {i: text_of(p) for i, p in enumerate(body_paragraphs, 1)}
        stories = ["word/document.xml"] + sorted(
            n for n in archive.namelist()
            if re.fullmatch(r"word/(?:footnotes|endnotes|comments|header\d+|footer\d+)\.xml", n)
        )
        for part in stories:
            root = document if part == "word/document.xml" else etree.fromstring(archive.read(part))
            tables = {id(t): i for i, t in enumerate(root.iter(W + "tbl"), 1)}
            for ordinal, paragraph in enumerate(root.iter(W + "p"), 1):
                body_i = body_ids.get(id(paragraph)) if part == "word/document.xml" else None
                style = paragraph.find("w:pPr/w:pStyle", NS)
                location = f"body.p{body_i:04d}" if body_i is not None else f"xml.p{ordinal:04d}"
                cell = next((p for p in paragraph.iterancestors() if p.tag == W + "tc"), None)
                if cell is not None:
                    row = cell.getparent()
                    table = next(p for p in row.iterancestors() if p.tag == W + "tbl")
                    row_i = table.findall("w:tr", NS).index(row) + 1
                    cell_i = row.findall("w:tc", NS).index(cell) + 1
                    location = f"table{tables[id(table)]}.row{row_i}.cell{cell_i}.xmlp{ordinal}"
                locations.append({
                    "source_part": part, "location": location,
                    "body_paragraph": body_i if body_i is not None else "",
                    "xml_paragraph": ordinal,
                    "word_para_id": paragraph.get("{" + NS["w14"] + "}paraId", ""),
                    "style": style.get(W + "val", "") if style is not None else "",
                    "text": text_of(paragraph),
                    "inside_table": cell is not None,
                })
        rels = etree.fromstring(archive.read("word/_rels/document.xml.rels"))
        relationships = {e.get("Id"): e.get("Target") for e in rels}
        for i, paragraph in enumerate(body_paragraphs, 1):
            for rid in paragraph.xpath(".//a:blip/@r:embed", namespaces=NS):
                target = relationships[rid]
                label = ""
                for following in range(i + 1, min(i + 6, len(body_paragraphs) + 1)):
                    candidate = body_text[following].strip()
                    if re.fullmatch(r"Figure\s+(?:S)?\d+", candidate):
                        label = candidate.replace(" ", "")
                        break
                record = {"location": f"body.p{i:04d}", "relationship": rid,
                          "figure": label, "source_part": "word/" + target}
                if label in {"Figure4", "Figure5", "FigureS4", "FigureS5"}:
                    destination = report / "qa" / "manuscript_source_media" / f"{label}_{Path(target).name}"
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    destination.write_bytes(archive.read("word/" + target))
                    image = Image.open(destination)
                    record.update(path=str(destination), width=image.width, height=image.height,
                                  sha256=sha256(destination))
                    # White-background copies are inspection aids only; original
                    # embedded image bytes remain separately preserved above.
                    rgba = image.convert("RGBA")
                    white = Image.new("RGBA", rgba.size, "white")
                    white.alpha_composite(rgba)
                    white.convert("RGB").save(destination.with_name(destination.stem + "_white.png"))
                images.append(record)
    return locations, body_text, images


def inventory(body: dict[int, str], document: Path) -> list[dict]:
    rows = []

    def add(p, claim, old, updater, status="PENDING_NEW_RESULTS", note="", panel="", source=None):
        rows.append({
            "claim_key": claim, "paragraph": p,
            "location": f"body.p{p:04d}" + (f" / {panel}" if panel else ""),
            "verbatim_claim": body[p], "old_value": old,
            "source": source or f"{document}::word/document.xml::body.p{p:04d}",
            "expected_updater": updater, "status": status,
            "new_value": "", "new_result_source": "", "proposed_replacement": "",
            "note": note,
        })

    composition = "ROOT: per-nucleus revised assignment/score table and per-mouse hepatocyte denominators"
    fcg = "FCG: new paired zonation model and outcome-agnostic threshold/selection summaries"
    network = "FCG: new four-context gene-lipid matrices and class membership ledgers"
    validation = "FCG: new male-periportal predictions joined to unchanged Mtmr11 KD lipid statistics"
    missing = "Deferred—James"
    preserve = "PRESERVE_INDEPENDENT_MEASUREMENT"
    for claim, old in [("fcg_nuclei_pericentral", "102,337"), ("fcg_nuclei_periportal", "86,439"),
                       ("composition_p_pericentral_text", "3.7e-12"), ("composition_p_periportal_text", "7.9e-7")]:
        add(91, claim, old, composition, missing,
            "Revised aggregate counts/UMIs are not nucleus numbers or composition denominators. No new value may be inferred from them.", "Figure4A / S4B-C")
    add(91, "assignment_description_results", "Cyp2e1 and Cyp2f2 mutually exclusive domains", "ROOT: confirmed expanded-marker rule and source provenance", "PENDING_MARKER_RULE_CONFIRMATION")
    for claim, old in [("fcg_dual_bias_total", "2,783"), ("fcg_gonad_bias_denominator", "4,276"),
                       ("fcg_dual_bias_percentage", "65%"), ("fcg_gonad_and_PC", "1,441"), ("fcg_gonad_and_PP", "1,342")]:
        add(93, claim, old, fcg, panel="Figure4B / TableS10",
            note="Historical replay establishes the WF05 zonation-model gonad-FDR set:4276 genes,2783 also zone-FDR,1441PC/1342PP. The unchanged WF02 alternative4103/2614 does not match the manuscript. Recompute this WF05 definition; keep the plot's |t|>3 categories separate.")
    add(93, "xenium_zonation_counts", "37/81 (46%) PC;10/81 (12%) PP;logFC>0.1", "Existing Xenium analysis; unchanged by revised sci-RNA zonation", preserve, panel="Figure4C-E / S4D / TableS11")
    add(95, "maldi_detected_and_category_counts", "373 detected;54 femalePC;50 femalePP;55 malePC;42 malePP;two inversions", "Existing MALDI/FISH analysis; independently defined spatial measurements", preserve, panel="Figure4F-I / S4E / TableS12")
    for claim, old in [("female_male_zonation_r", "0.80"), ("female_male_zonation_n_genes", "20,123"),
                       ("female_male_zonation_p", "<1e-300")]:
        add(99, claim, old, fcg, panel="Figure5A / TableS13")
    for claim, old in [("interaction_threshold", "5.4"), ("sex_specific_genes_total", "78"),
                       ("genes_male_PC", "30"), ("genes_male_PP", "25"),
                       ("genes_female_PC", "11"), ("genes_female_PP", "12")]:
        add(101, claim, old, fcg, panel="Figure5A-B / S5A / TableS13")
    add(101, "opposite_sex_gate", "|t|<1", "Retain prespecified descriptive selection gate", "PRESERVE_ANALYSIS_PARAMETER",
        "A weak opposite-sex statistic is not proof of absence/equivalence. Do not describe selected genes as proven biologically exclusive.")
    add(101, "male_gene_count_interpretation", "many more male-polarized genes;testosterone amplification", fcg, "REASSESS_AFTER_NEW_COUNTS",
        "Retain only if new category counts and interpretation support it; the mechanism remains a hypothesis.")
    add(101, "Aox3_inversion_and_Mtmr11_classification", "Aox3 malePP/femalePC;Mtmr11 malePP despite female-biased overall expression", fcg, "REASSESS_GENE_EXAMPLES",
        "Updated zonation statistics must support the spatial classification; overall whole-hepatocyte expression is a separate unchanged source.")
    for claim, old in [("network_gene_total_repeated", "78"), ("network_lipid_universe", "1,254"),
                       ("unique_associated_lipids", "442"), ("associated_lipid_percent", "35.2%"),
                       ("male_PC_gene_fraction", "38.5%"), ("male_PC_association_fraction", "54.7%"),
                       ("associated_male_PC", "287"), ("associated_male_PP", "137"),
                       ("associated_female_PC", "69"), ("associated_female_PP", "77"),
                       ("single_category_total", "332"), ("single_male_PC", "194"),
                       ("single_male_PP", "60"), ("single_female_PC", "49"),
                       ("single_female_PP", "29"), ("shared_category_total", "110")]:
        add(103, claim, old, network, panel="Figure5B-C / TableS14")
    add(103, "association_fraction_denominator_definition", "54.7% of gene-lipid associations", network,
        "RESOLVED_DISTINCT_EDGE_AND_SPECIES_DENOMINATORS",
        "Historical54.7%=387TPC gene-lipid edges/708all-context gene-lipid edges. The287/137/69/77 values count unique category-linked species, not edges. Keep these distinct denominators explicit.")
    add(103, "lipids_per_gene_approximation", "~10 malePC;~6 malePP;~6 femalePC;~6 femalePP", network,
        "RESOLVE_DENOMINATOR_BEFORE_REPLACEMENT",
        "Unique category-associated lipid species divided by gene count is not generally mean gene degree; use actual above-threshold gene-lipid edge counts if reporting mean lipids per gene.")
    add(103, "class_specific_exclusivity", "DG/PA/LPI male-only;sterols/HexCer female-only", network,
        "REASSESS_FROM_NEW_CLASS_FRACTIONS", "Correlations with zonated genes are not direct spatial lipid measurements.", "Figure5C")
    add(105, "network_representatives", "Fmo3 femalePP;Sult3a1 femalePC;Mtmr11 malePP;Kitl malePC", network,
        "REASSESS_GENE_EXAMPLES", "Check whether the historical examples still qualify; distinguish biological retention from the claim that they are now the highest-degree genes.", "Figure5D")
    add(105, "Mtmr11_six_display_partners", "PE(O-24:2_22:6);PE(O-22:1_22:6);PS(20:4_22:6);PE(O-20:1_22:6);PC(14:0_18:3);PC(32:2)", network,
        "PENDING_EDGE_REFRESH", "The embedded figure encodes r by edge width, without printed edge-r values. Its24 historical edge coefficients are separately inventoried from saved SIFs.", "Figure5D")
    add(107, "knockdown_time_hours", "48", "Original experimental protocol", "PRESERVE_EXPERIMENTAL_PARAMETER", panel="Figure5E-F / S5C")
    add(107, "Mtmr11_ULM", "-7.21", validation, panel="Figure5F / TableS15")
    add(107, "Mtmr11_ULM_p", "1.53e-12", validation, panel="Figure5F / TableS15")
    for lipid, corr, kd in [("PI(18:0_18:1)", "0.37", "-19.99"), ("PI(36:1)", "0.42", "-15.07"),
                            ("PS(18:0_18:1)", "0.26", "-18.20"), ("PC(18:0_18:1)", "0.54", "-7.22"),
                            ("PE(O-18:1_22:6)", "0.63", "-4.91"), ("BMP(22:6_22:6)", "-0.45", "8.70"),
                            ("DG(40:7)", "-0.47", "7.44")]:
        add(107, "validation_pair_" + lipid, f"{lipid}:r={corr};KD t={kd}", validation,
            "UPDATE_R_PRESERVE_EXPERIMENTAL_KD_T", "Only the in vivo r changes; retain the experimental KD t after verifying its unchanged source.", "Figure5F")
    add(122, "discussion_composition_interpretation", "fraction dedicated to PC functions vastly greater in gonadal females", composition,
        missing, "This qualitative conclusion also depends on unavailable revised nucleus-composition metadata.")
    add(145, "methods_old_two_marker_rule", ">1 Cyp2e1 and no Cyp2f2=>PC;>1 Cyp2f2 and no Cyp2e1=>PP", "ROOT: confirmed expanded-marker score/classification protocol", "PENDING_MARKER_RULE_CONFIRMATION",
        "Do not invent expanded-marker lists, thresholds, scaling, ambiguous-cell handling or provenance. Clarify whether old thresholds referred to normalized expression rather than raw counts.")
    add(147, "methods_composition_test", "two-sided Wilcoxon rank-sum", composition, missing,
        "The stated test needs per-mouse revised zone counts and appropriate hepatocyte denominators; preserve the test only if the planned estimand agrees.")
    add(167, "methods_zonation_input_definition", "PP/PC from Cyp2f2 and Cyp2e1 marker expression", "ROOT: confirmed expanded-marker protocol", "PENDING_MARKER_RULE_CONFIRMATION")
    add(167, "methods_entropy_selection", "entropy/stability/coverage/specificity sweep", fcg,
        "VERIFY_RETAINED_METHOD_UPDATE_OPERATING_POINT", "Retain the algorithm if it is the specified rerun; update the observed threshold rather than forcing5.4.")
    add(171, "methods_validation_consistency", "predicted correlation × KD t × -1", validation,
        "PRESERVE_FORMULA_UPDATE_PREDICTIONS", "The experimental lipid t-values stay fixed; all derived validation/consistency statistics use revised in vivo predictions.")
    add(280, "caption_4A_composition_and_strategy", "old isolation strategy;sample proportions;median/IQR;two-sided Wilcoxon", composition,
        missing, "Update the classification schematic only after marker provenance is confirmed. Xenium example images and100µm scale bar are independent.", "Figure4A")
    add(281, "caption_4B_axes_and_signs", "male vs female contrast but positive values labelled female", fcg,
        "VERIFY_AXIS_SIGN_CONVENTION", "Ensure the actual contrast orientation and prose agree; do not infer signs from label wording.", "Figure4B")
    add(291, "caption_5A_gene_total", "78", fcg, panel="Figure5A")
    add(292, "caption_5B_correlation_cutoff", "|r|>=0.65", "Retain prespecified display/association rule", "PRESERVE_ANALYSIS_PARAMETER", panel="Figure5B")
    add(292, "caption_5B_strongest_category_definition", "assign lipids to strongest category", network,
        "RESOLVE_CATEGORY_DEFINITION", "Caption strongest-category wording conflicts with explicit shared/multiple category assignments in body paragraph103. Also, the faithful5C butterfly groups lipids by abundance gonad_bias (lipid logFC sign), not by sex of correlated genes.", "Figure5B-C")
    add(294, "caption_5D_partner_count", "four subnetworks;six lipids each;24 edges", network,
        "VERIFY_NEW_DISPLAY_SELECTION", "Maintain selection rule or explicitly label historically retained examples; coefficients are not printed in the embedded network.", "Figure5D")
    add(295, "caption_5E_normalization", "per-gene maximum=1;two-sided Wilcoxon", "Unchanged KD RNA experiment", preserve, panel="Figure5E")
    add(296, "caption_5F_FCG_mouse_n", "43", validation, "CORRECT_SEX_SPECIFIC_SAMPLE_UNIT",
        "A within-male-periportal correlation uses the male cohort, not all43 FCG mice; expected22 if the unchanged mouse roster confirms it.", "Figure5F")
    add(296, "caption_5F_KD_replicates", "4", "Unchanged experimental KD design", "PRESERVE_KD_VERIFY_SOURCE", panel="Figure5F")
    add(296, "caption_5F_regression_interval", "95%", validation, "PRESERVE_ANALYSIS_PARAMETER_UPDATE_FIT", panel="Figure5F")
    add(337, "caption_S4B_old_assignment_scatter", "Cyp2e1 versus Cyp2f2;assignment thresholds", composition, missing,
        "A pseudobulk scatter cannot replace nucleus-level assignment/score QC.", "FigureS4B")
    add(338, "caption_S4C_marker_dotplot", "dot size=proportion expressing nuclei;colour=average scaled expression", composition, missing,
        "Neither fraction of expressing nuclei nor per-nucleus marker distributions can be reconstructed from aggregate gene counts.", "FigureS4C")
    add(344, "caption_S5A_threshold", "selected entropy threshold", fcg, panel="FigureS5A")
    add(346, "caption_S5C_cultured_marker_expression", "Control/Mtmr11 KD marker expression", "Unchanged independent primary-hepatocyte RNA measurements", preserve, panel="FigureS5C")
    for p, claim, old, owner in [(359, "supplementary_table_S10", "zonation FCG genes", fcg),
                                (360, "supplementary_table_S11", "zonation Xenium", "Existing independent Xenium data"),
                                (361, "supplementary_table_S12", "zonation MALDI", "Existing independent MALDI data"),
                                (362, "supplementary_table_S13", "sex-polarized zonation genes", fcg),
                                (363, "supplementary_table_S14", "gene-lipid-TF networks", network),
                                (364, "supplementary_table_S15", "network validation scores", validation)]:
        add(p, claim, old, owner, preserve if p in {360, 361} else "REGENERATE_RESULT_TABLE")
    # Image-only annotations were visually inspected on a white background.
    add(277, "figure4A_image_female_p", "p=3.7e-12", composition, missing,
        "Printed beside the female PP/PC rows; confirm the intended contrast before reconciling it with prose.", "Figure4A embedded image")
    add(277, "figure4A_image_male_p", "p=9.9e-05", composition, missing,
        "Image differs from paragraph91's7.9e-7. Resolve source/contrast; no new proportion estimate is recoverable from pseudobulk alone.", "Figure4A embedded image")
    add(277, "figure4B_additional_correlation_audit", "No r or p annotation printed in the embedded panel", fcg,
        "NEW_AUDIT_STATISTIC_NOT_LITERAL_REPLACEMENT", "The new4B correlation can be reported as a supplementary audit statistic; do not fabricate an old printed r/p.", "Figure4B embedded image")
    add(288, "figure5B_image_counts", "femalePP12;femalePC11;malePC30;malePP25", fcg, panel="Figure5B embedded image")
    add(288, "figure5C_image_class_fractions", "class-specific coloured/grey/white bar fractions;no per-class numeric labels", network,
        "REBUILD_FROM_EXPLICIT_CLASS_NUMERATORS_AND_DENOMINATORS", "Export every class fraction numerically in the replacement ledger; do not estimate values from pixel heights.", "Figure5C embedded image")
    add(288, "figure5E_image_KD_p", "p=0.03", "Unchanged independent KD RNA experiment", preserve, panel="Figure5E embedded image")
    add(288, "figure5F_image_ULM", "ULM score=-7.21", validation, panel="Figure5F embedded image")
    add(341, "figureS5A_image_threshold", "Optimal:5.4", fcg, panel="FigureS5A embedded image")
    return rows


def old_edge_inventory(root: Path) -> list[dict]:
    # These are source coefficients for the24 lipids pictured in the manuscript,
    # not numbers printed on its edges. Keep exact source values and locations.
    displayed = {
        ("Fmo3", "Ovaries_portal"): ["SM(16:2;O2/20:0)", "Cer(18:0;O3/22:0)", "LPE(O-20:0)", "Cer(17:0;O2/23:0;O)", "PI(17:0_18:2)", "Cer(18:0;O3/23:0)"],
        ("Mtmr11", "testes_portal"): ["PE(O-24:2_22:6)", "PE(O-22:1_22:6)", "PS(20:4_22:6)", "PE(O-20:1_22:6)", "PC(14:0_18:3)", "PC(32:2)"],
        ("Sult3a1", "ova_central"): ["PG(20:4_22:4)", "PG(22:4_22:4)", "PE(O-16:1_22:6)", "LPC(20:5)", "Cer(18:2;O2/25:1)", "FA(28:6;O)"],
        ("Kitl", "testes_central"): ["PE(18:1_22:1)", "PE(18:1_24:1)", "SM(18:0;O2/20:0)", "PE(24:1_18:2)", "PC(24:1_18:2)", "PE(22:1_20:4)"],
    }
    rows = []
    for (gene, context), lipids in displayed.items():
        path = analysis_path("historical_fcg05", root) / "networks" / f"cormat_{context}_full_SIF.csv"
        with path.open(newline="", encoding="utf-8-sig") as handle:
            source = list(csv.DictReader(handle))
        for lipid in lipids:
            match = [(i + 2, r) for i, r in enumerate(source) if r["node1"] == gene and r["node2"] == lipid]
            rows.append({"panel": "Figure5D", "gene": gene, "context": context, "lipid": lipid,
                "old_r": match[0][1]["weight"] if len(match) == 1 else "",
                "old_source": str(path), "old_source_csv_line": match[0][0] if len(match) == 1 else "",
                "status": "PENDING_NEW_R" if len(match) == 1 else "OLD_ANNOTATION_REQUIRES_EXACT_MATCH_REVIEW",
                "new_r": "", "new_source": "", "note": "Pictured partner visually transcribed from embedded Figure5; no coefficient is printed on the old edge."})
    return rows


def populate_saved_results(report: Path, rows: list[dict], edge_rows: list[dict]) -> None:
    """Read existing result ledgers; never fit models or change result files."""
    base = report.parent

    def read(path):
        if not path.exists():
            return []
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))

    by_key = {r["claim_key"]: r for r in rows}

    def put(key, value, source, old=None, status="Updated", note=None, prose=None):
        row = by_key[key]
        row.update(new_value=str(value), new_result_source=str(source), status=status)
        if old is not None:
            row["old_pipeline_value"] = str(old)
        if note is not None:
            row["note"] = note
        if prose is not None:
            row["proposed_replacement"] = prose

    preprocessing = base / "inputs/preprocessing_reference/preprocessing_02_fcg.R"
    if preprocessing.exists():
        draft = ("Hepatocyte zonation was assigned using periportal and pericentral marker signatures scored with UCell. "
                 "The provided preprocessing scripts exclude midlobular markers and retain markers detected in 5–95% of hepatocyte nuclei. "
                 "UCell maxRank is max(1000, round(median detected genes per nucleus)). PC and PP scores are standardized across all hepatocyte nuclei within each cohort; "
                 "their difference (PC z score minus PP z score) defines zonation. Nuclei with at least one detected retained marker at each pole are assigned PC when this difference is >0 and PP otherwise; other nuclei remain unassigned.")
        qualification = ("Draft supported by the staged preprocessing scripts and Zonation_Genes.csv, whose hashes are recorded in provenance_manifest.csv. "
                         "James must confirm the exact executed HPC inputs, final retained marker lists and runtime settings; these cannot be recovered from the supplied pseudobulk matrices. "
                         "Replace the old two-marker gates. This methods confirmation does not block the completed pseudobulk numerical update.")
        for key in ["assignment_description_results", "methods_old_two_marker_rule", "methods_zonation_input_definition"]:
            put(key, "Expanded-marker UCell assignment; final retained-marker/runtime confirmation pending", preprocessing,
                status="Methods draft; James confirmation required", note=qualification, prose=draft)

    put("caption_4B_axes_and_signs", "Current regenerated panel: x positive=PP; y positive=Testes (male), negative=Ovaries (female)",
        ROOT / "scripts/workflows/05_zonation_analysis.R", status="Correct caption to actual contrast signs",
        note="The regenerated panel directly plots portal_vs_Central and testes_vs_ova t statistics. No r/p was printed in the old embedded panel. The additional correlation audit uses these uncapped raw axes.",
        prose="The x-axis shows periportal minus pericentral bias and the y-axis shows testes minus ovaries bias; positive values indicate periportal and testes-biased expression, respectively.")
    put("caption_5B_strongest_category_definition", "Allow multiple sex/zone RNA-programme associations at the retained |r|>=0.65 cutoff", base / "fcg05/panel_source_data/Fig5C_lipid_species_membership_old_vs_new.csv",
        status="Correct category definition",
        note="Do not claim a unique strongest category: the body and exported memberships explicitly include shared lipids. Figure5C butterfly gonad grouping is lipid abundance bias, not sex of correlated RNA genes.")

    source = base / "comparisons/Fig4B_FDR_counts_old_vs_new.csv"
    counts = read(source)
    if counts:
        old, new = (next(r for r in counts if r["version"] == v) for v in ["historical", "revised"])
        for key, field in [("fcg_dual_bias_total", "dual_bias"), ("fcg_gonad_bias_denominator", "sex_biased"),
                           ("fcg_gonad_and_PC", "pericentral"), ("fcg_gonad_and_PP", "periportal")]:
            put(key, f"{int(new[field]):,}", source, old[field],
                note="WF05 zonation-model FDR sets reproduce the historical manuscript exactly; do not substitute WF02 or the plot's |t|>3 categories.")
        percentage = 100 * int(new["dual_bias"]) / int(new["sex_biased"])
        put("fcg_dual_bias_percentage", f"{percentage:.1f}%", source,
            old=f"{100*int(old['dual_bias'])/int(old['sex_biased']):.9f}%",
            prose=f"{int(new['dual_bias']):,} of {int(new['sex_biased']):,} ({percentage:.1f}%) gonadal sex-biased genes were also zonation-biased: {int(new['pericentral']):,} pericentral and {int(new['periportal']):,} periportal (FDR<=0.05).")
    source = base / "comparisons/fcg_selected_gene_counts.csv"
    selected = read(source)
    new = {r["context"]: r for r in selected if r["version"] == "revised_reoptimized"}
    if new:
        threshold = float(next(iter(new.values()))["threshold"])
        total = sum(int(r["n_genes"]) for r in new.values())
        for key in ["interaction_threshold", "figureS5A_image_threshold", "caption_S5A_threshold"]:
            put(key, f"{threshold:g}", source, "5.4")
        for key in ["sex_specific_genes_total", "network_gene_total_repeated", "caption_5A_gene_total"]:
            put(key, total, source, "78")
        for key, context in [("genes_male_PC", "TPC"), ("genes_male_PP", "TPP"),
                             ("genes_female_PC", "OPC"), ("genes_female_PP", "OPP")]:
            put(key, new[context]["n_genes"], source)
        put("figure5B_image_counts", ";".join(f"{c}={new[c]['n_genes']}" for c in ["OPP", "OPC", "TPC", "TPP"]), source)
        put("methods_entropy_selection", f"Retain the same objective; new operating threshold={threshold:g}", source,
            status="Retain method; update operating point")
        put("male_gene_count_interpretation", f"{sum(int(new[c]['n_genes']) for c in ['TPC','TPP'])} male vs {sum(int(new[c]['n_genes']) for c in ['OPC','OPP'])} female candidates", source,
            status="Reassess interpretation", note="Male predominance remains a descriptive count result; testosterone mechanism remains a hypothesis.")
    gene_table_source = base / "fcg05/t_table_zonation.csv"
    gene_table = {r["gene"]:r for r in read(gene_table_source)}
    if "Aox3" in gene_table and "Mtmr11" in gene_table:
        a, m = gene_table["Aox3"], gene_table["Mtmr11"]
        assert float(a["portalova_vs_centralova"]) < 0 < float(a["portaltestes_vs_centraltestes"])
        assert float(m["portaltestes_vs_centraltestes"]) > 0 and abs(float(m["portalova_vs_centralova"])) < 1 and float(m["testes_vs_ova"]) < 0
        put("Aox3_inversion_and_Mtmr11_classification", f"Retained: Aox3 female zone t={float(a['portalova_vs_centralova']):.3f},male zone t={float(a['portaltestes_vs_centraltestes']):.3f};Mtmr11 male zone t={float(m['portaltestes_vs_centraltestes']):.3f},female zone t={float(m['portalova_vs_centralova']):.3f},overall gonad t={float(m['testes_vs_ova']):.3f}", gene_table_source,
            status="Retained examples verified with revised zonation model",
            note="Zone contrasts are PP−PC; gonad contrast is Testes−Ovaries. Aox3 inversion and Mtmr11 male-PP/female-overall directions persist. Weak other-sex t remains a descriptive gate, not equivalence evidence.")
    membership_source = base / "fcg05/panel_source_data/Fig5C_lipid_species_membership_old_vs_new.csv"
    membership = read(membership_source)
    network_source = base / "comparisons/fcg_network_counts_old_vs_new.csv"
    network_counts = read(network_source)
    if membership and network_counts:
        m = {v: [r for r in membership if r["version"] == v] for v in ["historical", "revised"]}
        edge_counts = {(r["version"], r["context"]): int(r["gene_lipid_edges"]) for r in network_counts if r["type"] == "displayed"}
        selected_counts = {(r["version"], r["context"]): int(r["n_genes"]) for r in selected}
        n_assoc = lambda d: sum(int(r["n_contexts"]) > 0 for r in d)
        for key, fun in [("unique_associated_lipids", n_assoc),
                         ("single_category_total", lambda d: sum(int(r["n_contexts"]) == 1 for r in d)),
                         ("shared_category_total", lambda d: sum(int(r["n_contexts"]) > 1 for r in d)),
                         ("network_lipid_universe", len)]:
            put(key, fun(m["revised"]), membership_source, fun(m["historical"]))
        put("associated_lipid_percent", f"{100*n_assoc(m['revised'])/len(m['revised']):.1f}%", membership_source,
            f"{100*n_assoc(m['historical'])/len(m['historical']):.9f}%")
        denoms = []; new_means = []; new_unions = []
        for context, label in [("TPC", "male_PC"), ("TPP", "male_PP"), ("OPC", "female_PC"), ("OPP", "female_PP")]:
            hits = {v: sum(r[context] == "TRUE" for r in m[v]) for v in m}
            single = {v: sum(r["unique_context"] == context for r in m[v]) for v in m}
            put("associated_" + label, hits["revised"], membership_source, hits["historical"])
            put("single_" + label, single["revised"], membership_source, single["historical"])
            for version, sel_version in [("historical", "historical"), ("revised", "revised_reoptimized")]:
                ng = selected_counts[(sel_version, context)]; ne = edge_counts[(version, context)]
                denoms.append({"version": version, "context": context, "selected_genes": ng,
                    "gene_lipid_edges": ne, "unique_category_lipid_species": hits[version],
                    "mean_gene_degree_edges_per_gene": ne/ng,
                    "union_species_divided_by_genes_not_mean_degree": hits[version]/ng,
                    "edge_fraction_all_contexts": ne/sum(edge_counts[(version, c)] for c in new)})
                if version == "revised":
                    new_means.append(f"{context}:{ne/ng:.4g}"); new_unions.append(f"{context}:{hits[version]/ng:.4g}")
        write_csv(report / "figure5_network_denominator_audit.csv", denoms)
        old_edges = sum(edge_counts[("historical", c)] for c in new)
        new_edges = sum(edge_counts[("revised", c)] for c in new)
        put("male_PC_association_fraction", f"{100*edge_counts[('revised','TPC')]/new_edges:.1f}%", network_source,
            old=f"387/708={100*387/708:.9f}%",
            note=f"EDGE denominator: revised{edge_counts[('revised','TPC')]}/{new_edges}; historical387/{old_edges}. Unique-species counts use different denominators.")
        put("association_fraction_denominator_definition", f"TPC edges{edge_counts[('revised','TPC')]}/{new_edges}", network_source,
            old="TPC edges387/708", status="Resolved distinct edge/species denominators")
        put("male_PC_gene_fraction", f"{100*int(new['TPC']['n_genes'])/sum(int(r['n_genes']) for r in new.values()):.1f}%", source,
            old=f"30/78={100*30/78:.9f}%")
        put("lipids_per_gene_approximation", "True mean degree="+";".join(new_means)+";union-species/gene ratio="+";".join(new_unions),
            report / "figure5_network_denominator_audit.csv", status="Correct metric name and update",
            note="Both ratios are supplied; only edges/gene is mean partner degree. Historical means12.9/6.2/7.0/7.4167 differ from the manuscript's union-species/genes approximations.")
        male_only, female_only = [], []
        for lipid_class in sorted({r["lipid_class"] for r in m["revised"]}):
            d = [r for r in m["revised"] if r["lipid_class"] == lipid_class]
            male = any(r["TPC"] == "TRUE" or r["TPP"] == "TRUE" for r in d)
            female = any(r["OPC"] == "TRUE" or r["OPP"] == "TRUE" for r in d)
            if male and not female: male_only.append(lipid_class)
            if female and not male: female_only.append(lipid_class)
        put("class_specific_exclusivity", "Male-only:"+",".join(male_only)+";female-only:"+",".join(female_only), membership_source,
            status="Updated class association description", note="Associations with sex/zone gene programmes are not measured spatial lipid abundances.")
    examples_source = base / "comparisons/fcg_historical_four_examples_all_lipid_correlations.csv"
    examples = read(examples_source)
    if examples:
        context = {"testes_central":"TPC", "testes_portal":"TPP", "ova_central":"OPC", "Ovaries_portal":"OPP"}
        lookup = {(r["gene"], r["context"], r["lipid"]): r for r in examples}
        for edge in edge_rows:
            row = lookup.get((edge["gene"], context[edge["context"]], edge["lipid"]))
            if row:
                edge.update(new_r=row["r_new"], new_source=str(examples_source), status="Updated matched historical edge",
                            in_new_top6=row["top6_new"], selected_new=row["selected_new"])
        top6 = [r for r in examples if r["top6_new"] == "TRUE"]
        write_csv(report / "figure5D_new_top6_partners.csv", top6)
        put("Mtmr11_six_display_partners", ";".join(r["lipid"] for r in top6 if r["gene"] == "Mtmr11"), examples_source,
            status="Update partners if retaining top-six rank rule")
        put("network_representatives", "Historical Fmo3/Sult3a1/Mtmr11/Kitl contexts all retained" if all(r["selected_new"] == "TRUE" for r in examples) else "See historical-example eligibility audit",
            base / "comparisons/fcg_historical_four_examples_audit.csv", status="Retain prespecified examples; refresh edges",
            note="Do not claim these are now the highest-degree genes without comparing new degree ranks.")
        put("caption_5D_partner_count", "Four examples;six top partners each;24 displayed edges", examples_source,
            status="Retain display size; refresh selected partners")
    fractions_source = base / "fcg05/panel_source_data/Fig5C_faithful_butterfly_proportions_old_vs_new.csv"
    fractions = read(fractions_source)
    if fractions:
        groups = sorted({(r["lipid_class"], r["gonad_bias"]) for r in fractions})
        display = {r["lipid_class"]:r["displayed_class"] for r in read(base / "fcg05/panel_source_data/Fig5C_explicit_gene_context_class_fractions.csv")}
        lookup = {(r["version"],r["lipid_class"],r["gonad_bias"],r["zonation"]): r for r in fractions}
        all_rows = []
        for lipid_class, gonad in groups:
            for zone in ["pericentral","periportal","both","no_association"]:
                out = {"lipid_class":lipid_class,"lipid_abundance_gonad_bias":gonad,"zone_association":zone,
                       "displayed_class":display.get(lipid_class, "")}
                for v, label in [("historical","old"),("revised","new")]:
                    r = lookup.get((v,lipid_class,gonad,zone))
                    total = next(x["total"] for x in fractions if x["version"]==v and x["lipid_class"]==lipid_class and x["gonad_bias"]==gonad)
                    out.update({label+"_count":r["count"] if r else 0,label+"_denominator":total,label+"_fraction":r["proportion"] if r else 0})
                out["source"] = str(fractions_source)
                out["note"] = "Gonad grouping is lipid abundance logFC sign, not sex of correlated RNA programme."
                all_rows.append(out)
        write_csv(report / "figure5C_class_fraction_replacements.csv", all_rows)
        put("figure5C_image_class_fractions", f"See{len(all_rows)} explicit old/new class×abundance-gonad×zone rows", fractions_source,
            status="Updated fractions with explicit denominators",
            note="Faithful plotted bars and four-context RNA-programme memberships are separate ledgers.")
    validation_source = base / "comparisons/Fig5F_prespecified_male_PP_old_vs_new.csv"
    values = read(validation_source)
    if values:
        old, revised = (next(r for r in values if r["version"]==v) for v in ["historical","revised"])
        put("Mtmr11_ULM", f"{float(revised['ULM_score']):.2f}", validation_source, old["ULM_score"])
        put("Mtmr11_ULM_p", f"{float(revised['ULM_p']):.3g}", validation_source, old["ULM_p"])
        put("figure5F_image_ULM", f"ULM score={float(revised['ULM_score']):.2f}", validation_source,
            note="ULM uses all661 KD targets;345 matched lipids appear in the scatter. Temporary matched-only ULM values are superseded.")
        put("caption_5F_FCG_mouse_n", revised["n_reference_mice"], validation_source, old["n_reference_mice"],
            status="Correct caption reference-mouse sample size",
            note=f"22 reference mice, {revised['n_lipids']} matched lipid points and {revised['n_KD_background']} ULM background targets are distinct units.")
        pairs_source = base / "fcg07/Mtmr11_male_PP_prespecified_source.csv"
        lookup = {r["lipid"]:r for r in read(pairs_source)}
        for row in rows:
            if row["claim_key"].startswith("validation_pair_"):
                lipid = row["claim_key"].removeprefix("validation_pair_")
                if lipid in lookup:
                    r = lookup[lipid]
                    put(row["claim_key"],f"{lipid}:r={float(r['r']):.2f};KD t={float(r['KD_t']):.2f}",pairs_source,
                        status="Update in-vivo r; preserve unchanged experimental KD t")
    for value in read(report / "qa/independent_numeric_summary.csv"):
        put(value["metric_key"],value["new_value"],value["source_file"],value["old_value"],
            status="Independently replayed",note=value["note"] or None)
    for row in rows:
        if row["status"].startswith("PRESERVE_") and not row["new_value"]:
            row["new_value"] = row["old_value"]
            row["proposed_replacement"] = "Retain; this independent experiment or analysis parameter was not changed by the pseudobulk update."
    tables = read(base / "tables/table_manifest.csv")
    for table in ["S10", "S13", "S14", "S15"]:
        sources = [r["file"] for r in tables if r["table_id"] == table]
        if sources:
            assert all(Path(p).exists() for p in sources)
            put("supplementary_table_" + table, "Regenerated source available; replace corresponding supplementary-table content", ";".join(sources),
                status="Regenerated source ready for manuscript table assembly",
                note="Source CSVs/networks are provided separately. The Word manuscript and supplementary workbook were not edited.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", type=Path, default=ROOT / "manuscript/Cleland and Dugourd et al.docx")
    parser.add_argument("--output", type=Path, default=ROOT / "results/zonation_marker_update_20261116/report")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    before = sha256(args.manuscript)
    locations, body, images = extract_docx(args.manuscript, args.output)
    rows = inventory(body, args.manuscript.resolve())
    edge_rows = old_edge_inventory(ROOT)
    populate_saved_results(args.output, rows, edge_rows)
    updates_file = args.output / "manuscript_numeric_updates.json"
    if updates_file.exists():
        updates = json.loads(updates_file.read_text(encoding="utf-8"))
        unknown = set(updates) - {r["claim_key"] for r in rows}
        if unknown:
            raise ValueError(f"Unrecognised claim keys in update ledger: {sorted(unknown)}")
        for row in rows:
            if row["claim_key"] in updates:
                row.update(updates[row["claim_key"]])
    for row in rows:
        section = row["location"].split(" / ", 1)
        row.update(paragraph_id=f"body.p{row['paragraph']:04d}",
                   section_or_panel=section[1] if len(section) > 1 else
                   "Methods" if 130 <= row["paragraph"] < 229 else
                   "Main text" if row["paragraph"] < 229 else "Caption or supplementary table",
                   exact_text=row["verbatim_claim"], old_values=row["old_value"],
                   metric_key=row["claim_key"], scope_status=row["status"])
    preferred = ["paragraph_id", "section_or_panel", "exact_text", "old_values", "metric_key", "source", "scope_status"]
    fields = preferred + list(dict.fromkeys(k for row in rows for k in row if k not in preferred))
    write_csv(args.output / "manuscript_numeric_inventory.csv", rows, fields)
    write_csv(args.output / "manuscript_paragraph_locations.csv", locations)
    write_csv(args.output / "manuscript_embedded_figure_sources.csv", images)
    write_csv(args.output / "figure5D_old_edge_coefficients.csv", edge_rows)
    replacements = [{
        "manuscript_location": row["paragraph_id"], "panel": row["section_or_panel"],
        "metric": row["claim_key"], "old_manuscript_value": row["old_value"],
        "old_pipeline_value": row.get("old_pipeline_value", ""),
        "new_value": row["new_value"], "status": row["status"],
        "notes": row["note"], "source_file": row["new_result_source"] or row["source"],
    } for row in rows]
    replacements.extend({
        "manuscript_location": "body.p0288 / embedded Figure5D", "panel": "Figure5D",
        "metric": f"edge_r:{row['context']}:{row['gene']}:{row['lipid']}",
        "old_manuscript_value": "Not printed; encoded by edge width",
        "old_pipeline_value": row["old_r"], "new_value": row["new_r"],
        "status": row["status"], "notes": row["note"],
        "source_file": row["new_source"] or row["old_source"],
    } for row in edge_rows)
    write_csv(args.output / "new_replacements.csv", replacements)
    candidates = []
    for row in locations:
        if re.search(r"zonat|pericentral|periportal|Mtmr11|Figure\s*[45]", row["text"], re.I):
            candidate = dict(row)
            candidate["numeric_tokens"] = ";".join(re.findall(r"(?<![A-Za-z])[-+<>≤≥]?\s*\d+(?:,\d{3})*(?:\.\d+)?(?:e[-+]?\d+)?%?", row["text"], re.I))
            candidate["curated_claim_keys"] = ";".join(r["claim_key"] for r in rows if r["paragraph"] == row["body_paragraph"] and row["source_part"] == "word/document.xml")
            candidate["review_scope"] = "CURATED" if candidate["curated_claim_keys"] else "BACKGROUND_OR_ORTHOGONAL_OR_WORKING_NOTE_REVIEW"
            candidates.append(candidate)
    write_csv(args.output / "manuscript_zonation_context_scan.csv", candidates)
    assert sha256(args.manuscript) == before, "Source manuscript was changed"
    metadata = {"manuscript": str(args.manuscript.resolve()), "sha256_before": before,
                "sha256_after": sha256(args.manuscript), "inventory_rows": len(rows),
                "body_paragraphs": len(body), "OOXML_table_paragraphs": sum(r["inside_table"] for r in locations),
                "locator_definition": "One-based direct w:body paragraph index, including blanks; not a rendered page number",
                "document_unchanged": True, "rendering": "No PDF rendered; source XML and embedded PNGs inspected only"}
    (args.output / "manuscript_audit_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Manuscript inventory: {len(rows)} claims; {len(body)} body paragraphs; source SHA256 unchanged.")
    print(f"Deferred—James composition/nucleus-QC rows: {sum(r['status'] == 'Deferred—James' for r in rows)}")


if __name__ == "__main__":
    main()
