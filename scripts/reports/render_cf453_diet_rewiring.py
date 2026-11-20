"""Render and inspect the three independent diet panel exports."""
import argparse
import csv
import importlib.metadata as metadata
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "maintenance"))
from repository_paths import analysis_current, resolve_manifest_path
import xml.etree.ElementTree as ET

import pdfplumber
import pypdfium2 as pdfium
from PIL import Image, ImageDraw
from pypdf import PdfReader


def review(output):
    root = Path(output).resolve()
    qa = root / "QA" / "visual_review"
    qa.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((root / "panel_manifest.csv").open()))
    pdfs = [row for row in rows if row["format"] == "PDF"]
    assert len(pdfs) == 3 and len(rows) == 9
    assert len({row["source_table"] for row in rows}) == 3
    records = []
    for row in pdfs:
        path = resolve_manifest_path(row["file"], root)
        document = pdfium.PdfDocument(str(path))
        page = document[0]
        image = page.render(scale=1.5).to_pil().convert("RGB")
        image.save(qa / (path.stem + "_pdf.png"))
        text = PdfReader(str(path)).pages[0].extract_text()
        with pdfplumber.open(path) as pdf:
            p = pdf.pages[0]
            outside = [char for char in p.chars if char["text"].strip() and (
                char["x0"] < -1 or char["x1"] > p.width + 1
                or char["top"] < -1 or char["bottom"] > p.height + 1)]
        raster_nodes = [node for node in ET.parse(path.with_suffix(".svg")).iter()
                        if node.tag.endswith("}image")]
        png = Image.open(path.with_suffix(".png"))
        expected = (round(float(row["width_in"]) * 300), round(float(row["height_in"]) * 300))
        assert len(document) == 1 and text.strip()
        assert not raster_nodes and not outside, path.stem + " vector/bounds checks"
        assert all(abs(x - y) <= 1 for x, y in zip(png.size, expected))
        assert resolve_manifest_path(row["source_table"], root).exists()
        records.append(dict(panel=path.stem, pdf_pages=len(document), text_characters=len(text),
                            pdf_bytes=path.stat().st_size, svg_contains_raster_images=False,
                            png_pixels=png.size, expected_png_pixels_300dpi=expected,
                            out_of_page_text_characters=len(outside)))
        page.close()
        document.close()
    sheet = Image.new("RGB", (1920, 850), "#eeeeee")
    draw = ImageDraw.Draw(sheet)
    for index, row in enumerate(pdfs):
        path = resolve_manifest_path(row["file"], root)
        image = Image.open(qa / (path.stem + "_pdf.png"))
        image.thumbnail((625, 800))
        sheet.paste(image, (index * 640 + (640 - image.width) // 2, 35))
        draw.text((index * 640 + 12, 10), path.stem, fill="black")
    sheet.save(qa / "contact_sheet.png")
    report = dict(panel_groups=3, graphic_files=9,
                  versions={x: metadata.version(x) for x in ("pypdfium2", "pdfplumber", "pypdf", "Pillow")},
                  records=records)
    (qa / "render_manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(analysis_current("cf453_rewiring")))
    review(parser.parse_args().output)
