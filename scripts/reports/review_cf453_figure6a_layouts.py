"""Render the two Figure6A candidates and optionally package checked exports."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "maintenance"))
from repository_paths import analysis_path, resolve_manifest_path
import xml.etree.ElementTree as ET
import zipfile

import pdfplumber
import pypdfium2 as pdfium
from PIL import Image, ImageDraw, ImageFont
from pypdf import PdfReader


def review(root):
    root = Path(root).resolve()
    qa = root / "QA" / "visual_review"
    qa.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((root / "panel_manifest.csv").open()))
    pdfs = [r for r in rows if r["format"] == "PDF"]
    assert len(rows) == 6 and len(pdfs) == 2
    records = []
    for r in pdfs:
        path = resolve_manifest_path(r["file"], root)
        doc = pdfium.PdfDocument(str(path))
        page = doc[0]
        picture = page.render(scale=1.7).to_pil().convert("RGB")
        picture.save(qa / (path.stem + "_pdf.png"))
        text = PdfReader(str(path)).pages[0].extract_text()
        with pdfplumber.open(path) as pdf:
            p = pdf.pages[0]
            outside = [c for c in p.chars if c["text"].strip() and (
                c["x0"] < -1 or c["x1"] > p.width + 1 or c["top"] < -1 or c["bottom"] > p.height + 1)]
        raster = [n for n in ET.parse(path.with_suffix(".svg")).iter() if n.tag.endswith("}image")]
        png = Image.open(path.with_suffix(".png"))
        expected = (round(float(r["width_in"]) * 300), round(float(r["height_in"]) * 300))
        assert len(doc) == 1 and text.strip()
        assert not outside, (path.name, "text outside page", len(outside))
        assert not raster, (path.name, "raster SVG")
        assert all(abs(x - y) <= 1 for x, y in zip(png.size, expected))
        records.append(dict(panel=path.stem, pdf_pages=len(doc), text_characters=len(text),
                            out_of_page_text=0, vector_svg=True, png_pixels=png.size, dpi=300))
        page.close()
        doc.close()
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 25)
    sheet = Image.new("RGB", (1800, 1250), "#eeeeee")
    draw = ImageDraw.Draw(sheet)
    for j, r in enumerate(pdfs):
        im = Image.open(qa / (resolve_manifest_path(r["file"], root).stem + "_pdf.png"))
        im.thumbnail((880, 1180))
        sheet.paste(im, (j * 900 + (900 - im.width) // 2, 55))
        draw.text((j * 900 + 20, 15), r["option"], font=font, fill="black")
    sheet.save(qa / "two_layouts_comparison.png")
    (qa / "render_manifest.json").write_text(json.dumps(dict(panel_groups=2, records=records), indent=2) + "\n")
    print(json.dumps(records, indent=2))
    return rows


def package(root, rows):
    root = Path(root).resolve()
    payload = {}
    portable = []
    for row in rows:
        r = dict(row)
        for key, folder in (("file", "panels"), ("source_table", "source_data")):
            path = resolve_manifest_path(r[key], root); name = f"{folder}/{path.name}"
            payload[name] = path.read_bytes(); r[key] = name
        portable.append(r)
    for name in ("README.md", "common_gene_ledger.csv", "provenance_sha256.csv",
                 "retained_original_sha256.csv", "sessionInfo.txt"):
        payload[name] = (root / name).read_bytes()
    checks = root / "QA" / "independent_checks.csv"
    assert checks.exists(), "Numerical/source verification must precede packaging"
    payload["verification/independent_checks.csv"] = checks.read_bytes()
    payload["verification/render_manifest.json"] = (root / "QA/visual_review/render_manifest.json").read_bytes()
    buf = io.StringIO(); writer = csv.DictWriter(buf, fieldnames=list(portable[0]))
    writer.writeheader(); writer.writerows(portable)
    payload["panel_manifest.csv"] = buf.getvalue().encode()
    buf = io.StringIO(); writer = csv.writer(buf); writer.writerow(("path", "sha256", "bytes"))
    for name, contents in sorted(payload.items()):
        writer.writerow((name, hashlib.sha256(contents).hexdigest(), len(contents)))
    payload["SHA256SUMS.csv"] = buf.getvalue().encode()
    dest = root / "Figure6A_two_layout_options.zip"
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for name, contents in sorted(payload.items()):
            z.writestr(name, contents)
    with zipfile.ZipFile(dest) as z:
        assert z.testzip() is None
        assert all(hashlib.sha256(z.read(name)).digest() == hashlib.sha256(data).digest() for name, data in payload.items())
    summary = dict(path=str(dest), panel_groups=2, graphic_files=6, source_tables=2,
                   bytes=dest.stat().st_size, verification="CRC and SHA256 passed")
    (root / "bundle_verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", default=str(analysis_path("current_root") / "cf453_layouts"))
    p.add_argument("--package", action="store_true")
    a = p.parse_args(); rows = review(a.output)
    if a.package:
        package(a.output, rows)
