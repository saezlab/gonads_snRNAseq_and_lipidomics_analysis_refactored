"""Validate final HTML/PDF local artifact links without changing reports."""
import csv
import hashlib
import json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlparse
from pypdf import PdfReader

report = Path(__file__).resolve().parents[1]
output = report / "qa/layout_final"
output.mkdir(exist_ok=True)
html_path = report / "zonation_update_report.html"
pdf_path = report / "zonation_update_report.pdf"
initial = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [html_path, pdf_path]}

class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.hrefs = []
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for key, value in attrs:
                if key == "href": self.hrefs.append(value)

parser = Links(); parser.feed(html_path.read_text())
records = []
for href in parser.hrefs:
    url = urlparse(href)
    local = not url.scheme or url.scheme == "file"
    target = Path(unquote(url.path)) if url.scheme == "file" else report / unquote(url.path)
    records.append(dict(surface="html", href=href, local=local,
                        relative=not url.scheme and not href.startswith("/"),
                        exists=target.exists() if local else None))
for page_number, page in enumerate(PdfReader(pdf_path).pages, 1):
    for annotation in page.get("/Annots", []):
        item = annotation.get_object()
        if item.get("/Subtype") != "/Link": continue
        action = item.get("/A", {})
        href = str(action.get("/URI", ""))
        url = urlparse(href)
        local = not url.scheme or url.scheme == "file"
        target = Path(unquote(url.path)) if url.scheme == "file" else report / unquote(url.path)
        records.append(dict(surface=f"pdf_page_{page_number}", href=href, local=local,
                            relative=not url.scheme and not href.startswith("/"),
                            exists=target.exists() if local else None))
with (output / "artifact_links.csv").open("w") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(records[0]))
    writer.writeheader(); writer.writerows(records)
html_records = [r for r in records if r["surface"] == "html"]
summary = dict(html_links=len(html_records), html_all_relative=all(r["relative"] for r in html_records),
               pdf_links=sum(r["surface"].startswith("pdf_") for r in records),
               missing_local_links=[r for r in records if r["local"] and not r["exists"]],
               pdf_links_portability="PDF uses absolute local file URIs; HTML uses portable relative links",
               report_files_unchanged=all(hashlib.sha256((report / n).read_bytes()).hexdigest() == sha for n, sha in initial.items()))
(output / "link_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
assert summary["html_all_relative"] and not summary["missing_local_links"]
assert summary["pdf_links"] > 0 and summary["report_files_unchanged"]
print(json.dumps(summary, indent=2))
