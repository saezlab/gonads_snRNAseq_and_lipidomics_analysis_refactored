"""Read-only PDF/trace review; writes audit records and rendered previews only."""
import argparse, csv, hashlib, json, subprocess
from pathlib import Path
import pdfplumber
from pypdf import PdfReader
from PIL import Image, ImageDraw, ImageFont

parser=argparse.ArgumentParser()
parser.add_argument("--pdf",type=Path,default=Path(__file__).resolve().parents[1]/"zonation_update_report.pdf")
parser.add_argument("--label",default="draft")
args=parser.parse_args()
pdf=args.pdf.resolve(); folder=pdf.parent/"qa"/("layout_"+args.label)
folder.mkdir(parents=True,exist_ok=True)
before=hashlib.sha256(pdf.read_bytes()).hexdigest()
reader=PdfReader(pdf); count=len(reader.pages)
metrics=[]; texts=[]
with pdfplumber.open(pdf) as doc:
    for index,page in enumerate(doc.pages,1):
        text=reader.pages[index-1].extract_text() or ""
        chars=[c for c in page.chars if c.get("text","").strip()]
        outside=[c for c in chars if c["x0"]<-.5 or c["x1"]>page.width+.5 or c["top"]<-.5 or c["bottom"]>page.height+.5]
        weird=[c for c in chars if c["text"] in ("\u25a0","\ufffd") or "ZapfDingbats" in c["fontname"]]
        sizes=[float(c["size"]) for c in chars]
        links=[a for a in (reader.pages[index-1].get("/Annots") or []) if a.get_object().get("/Subtype")=="/Link"]
        metrics.append(dict(page=index,text_characters=len(text),visible_characters=len(chars),
            out_of_page_characters=len(outside),suspicious_glyphs=len(weird),min_font_pt=min(sizes),
            max_font_pt=max(sizes),hyperlinks=len(links),width_pt=page.width,height_pt=page.height))
        texts.append(dict(page=index,text=text))
with (folder/"page_metrics.csv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
(folder/"page_text.json").write_text(json.dumps(texts,indent=2))
poppler="/Users/dugourd/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm"
subprocess.run([poppler,"-r","100","-png",str(pdf),str(folder/"page")],check=True,capture_output=True)
images=sorted(folder.glob("page-*.png"))
assert len(images)==count
font=ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf",22)
def sheet(indices,name,cols,width,maxheight):
    opened=[]
    for i in indices:
        image=Image.open(images[i-1]).convert("RGB");image.thumbnail((width-14,maxheight-36))
        opened.append((i,image))
    rows=(len(opened)+cols-1)//cols
    canvas=Image.new("RGB",(cols*width,rows*maxheight),"#e9edf1")
    draw=ImageDraw.Draw(canvas)
    for j,(i,image) in enumerate(opened):
        x=(j%cols)*width;y=(j//cols)*maxheight
        draw.text((x+10,y+5),f"Page {i}",font=font,fill="#202631")
        canvas.paste(image,(x+(width-image.width)//2,y+34))
    canvas.save(folder/name)
sheet(list(range(1,count+1)),"all_pages_contact.png",3,500,740)
representative=list(dict.fromkeys([1,min(2,count),min(3,count),min(5,count),min(9,count),count]))
sheet(representative,"representative_contact.png",2,800,1170)
trace=list(csv.DictReader((pdf.parent/"panel_trace.csv").open()))
cf=[]
for row in trace:
    if row["panel"] in ("Figure6A","Figure6B","Figure6C"):
        for name in filter(None,row["generated_files"].split("; ")):
            path=Path(name)
            if path.suffix.lower()==".pdf":
                cf.append(dict(panel=row["panel"],path=name,exists=path.exists()))
assert len(cf)==5 and len({r["path"] for r in cf})==5 and all(r["exists"] for r in cf)
with (folder/"cf453_panel_trace_check.csv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=list(cf[0]));w.writeheader();w.writerows(cf)
after=hashlib.sha256(pdf.read_bytes()).hexdigest()
summary=dict(pdf=str(pdf),pdf_sha256=before,unchanged_during_review=before==after,pages=count,
    total_out_of_page_characters=sum(r["out_of_page_characters"] for r in metrics),
    total_suspicious_glyphs=sum(r["suspicious_glyphs"] for r in metrics),
    min_font_pt=min(r["min_font_pt"] for r in metrics),
    blank_pages=sum(r["text_characters"]<40 for r in metrics),pdf_hyperlinks=sum(r["hyperlinks"] for r in metrics),
    cf453_panel_paths_checked=len(cf),cf453_paths_all_exist=True,render_dpi=100,
    representative_pages=representative,report_stage=args.label,
    interpretation=("Final completed numerical report: layout and artifact mapping reviewed; independent statistical replay is documented separately."
                    if args.label=="final" else "Layout audit of the draft report; numerical replacement review is documented separately."))
(folder/"audit_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
assert summary["unchanged_during_review"] and summary["total_out_of_page_characters"]==0
assert summary["total_suspicious_glyphs"]==0 and summary["blank_pages"]==0
print(json.dumps(summary,indent=2))
