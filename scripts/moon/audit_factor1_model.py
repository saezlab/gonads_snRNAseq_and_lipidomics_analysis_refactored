#!/usr/bin/env python3
"""Read the existing MOFA HDF5 files; never retrain or modify the model.

Requires numpy and h5py. Run after prepare_factor1_inputs.R from the repo root.
The audit aligns saved Factor1 weights with the file's native factor order,
records the corresponding sample coordinate, and reports descriptive consistency
across existing factor-number fits (not a bootstrap or seed-stability test).
"""
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/cosmos_factor1_pilot/inputs"
MODEL = ROOT / "results/mofa/RNA/mofa_res_12factor.hdf5"


def write_csv(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def decode(values):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in values]


def metadata():
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(ROOT / "data/Suppl_Table_1_OverviewData.xlsx") as z:
        strings = ["".join(t.text or "" for t in e.findall(".//s:t", ns))
                   for e in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("s:si", ns)]
        sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in sheet.findall("s:sheetData/s:row", ns):
            values = {}
            for cell in row:
                col = "".join(c for c in cell.attrib["r"] if c.isalpha())
                v = cell.find("s:v", ns)
                text = v.text if v is not None else ""
                values[col] = strings[int(text)] if cell.attrib.get("t") == "s" else text
            rows.append(values)
        header = rows[0]
        return {r["B"]: {header[k]: v for k, v in r.items()} for r in rows[1:] if r.get("B")}


raw = list(csv.DictReader((OUT / "factor1_all_views_raw.csv").open()))
alignments = []
with h5py.File(MODEL) as f:
    for view in ("Hepatocytes", "lipidomic"):
        saved = {r["feature_id"]: float(r["raw_value"]) for r in raw if r["view"] == view}
        features = decode(f[f"features/{view}"][()])
        target = np.array([saved[x] for x in features])
        native = f[f"expectations/W/{view}"][()]
        correlations = np.array([np.corrcoef(target, row)[0, 1] for row in native])
        idx = int(np.argmax(np.abs(correlations)))
        sign = float(np.sign(correlations[idx]))
        err = float(np.max(np.abs(target - native[idx] * sign)))
        assert err < 1e-10, (view, idx, err)
        alignments.append(dict(view=view, native_factor_index_1based=idx + 1,
                               native_to_saved_sign=int(sign), pearson_r=float(correlations[idx]),
                               max_absolute_difference=err))
    assert len({(x["native_factor_index_1based"], x["native_to_saved_sign"]) for x in alignments}) == 1
    idx = alignments[0]["native_factor_index_1based"] - 1
    native_to_saved = alignments[0]["native_to_saved_sign"]
    samples = decode(f["samples/single_group"][()])
    z = f["expectations/Z/single_group"][()][idx] * native_to_saved
    views = decode(f["views/views"][()])
    variance = f["variance_explained/r2_per_factor/single_group"][()]
    total = f["variance_explained/r2_total/single_group"][()]
    factors = f["training_stats/number_factors"][()]
    training_names = decode(f["training_opts"].attrs["names"])
    training = dict(zip(training_names, f["training_opts"][()].tolist()))
    reference_weights = {v: f[f"expectations/W/{v}"][()][idx] * native_to_saved for v in ("Hepatocytes", "lipidomic")}
    reference_features = {v: decode(f[f"features/{v}"][()]) for v in reference_weights}

meta = metadata()
assert all(s in meta for s in samples)
means = {g: float(np.mean([v for s, v in zip(samples, z) if meta[s]["gonadal_type"] == g]))
         for g in ("Ovaries", "Testes")}
orientation_multiplier = 1 if means["Testes"] > means["Ovaries"] else -1
coordinates = [dict(sample=s, group=meta[s]["group"], gonadal_type=meta[s]["gonadal_type"],
                    allosome_type=meta[s]["allosome_type"], factor1_coordinate_saved=float(v),
                    factor1_coordinate_testes_positive=float(v * orientation_multiplier)) for s, v in zip(samples, z)]
write_csv("factor1_sample_coordinates.csv", coordinates)
write_csv("factor1_model_alignment.csv", alignments)
write_csv("factor1_variance_explained.csv", [dict(view=v, factor1_variance_percent=float(variance[i, idx]),
          all_factors_variance_percent=float(total[i])) for i, v in enumerate(views)])
stability = []
for path in sorted(MODEL.parent.glob("mofa_res_*factor.hdf5")):
    with h5py.File(path) as f:
        other_samples = decode(f["samples/single_group"][()])
        if set(other_samples) != set(samples):
            continue
        other_z = f["expectations/Z/single_group"][()][:, [other_samples.index(s) for s in samples]]
        corr = np.array([np.corrcoef(z, row)[0, 1] for row in other_z])
        k = int(np.argmax(np.abs(corr)))
        sign = float(np.sign(corr[k]))
        row = dict(model=path.name, retained_factors=len(other_z), matching_native_factor_1based=k + 1,
                   aligned_factor_coordinate_r=float(abs(corr[k])))
        for view in reference_weights:
            other_features = decode(f[f"features/{view}"][()])
            assert set(other_features) == set(reference_features[view]), (path.name, view)
            feature_index = {x: i for i, x in enumerate(other_features)}
            aligned_weights = f[f"expectations/W/{view}"][()][k][[feature_index[x] for x in reference_features[view]]] * sign
            row[f"{view}_loading_r"] = float(np.corrcoef(reference_weights[view], aligned_weights)[0, 1])
        stability.append(row)
write_csv("factor1_existing_model_consistency.csv", stability)
audit = dict(model=str(MODEL.relative_to(ROOT)), initial_factors=int(factors[0]), retained_factors=int(factors[-1]),
             training_iterations_recorded=len(factors), training_options=training, sample_count=len(samples),
             native_factor_index_1based=idx + 1, native_to_saved_sign=native_to_saved,
             saved_coordinate_group_means=means, saved_to_testes_positive_multiplier=orientation_multiplier,
             alignment=alignments,
             caveat="Existing factor-number fits only; no bootstrap, repeated-seed, or donor-resampling validation.")
(OUT / "factor1_model_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
print(json.dumps(audit, indent=2))
