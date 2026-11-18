from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "maintenance"))
from repository_paths import analysis_path
import argparse, tempfile
import re, html
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, KeepTogether, Preformatted
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas

parser=argparse.ArgumentParser(description="Render the WP2 lipid-overlap and pooled Figure 3 addendum and append its scientific vector figures.")
parser.add_argument('--input',type=Path,default=analysis_path('historical_wf12') / 'followup/README.md')
parser.add_argument('--output',type=Path,default=None)
parser.add_argument('--font-dir',type=Path,default=Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype')
args=parser.parse_args()
src=args.input.resolve(strict=True)
final=(args.output or src.parent/'WP2_lipid_overlap_and_pooled_Figure3.pdf').resolve()
work=tempfile.TemporaryDirectory(prefix='wdsd_pooled_pdf_')
out=Path(work.name)/'report_text.pdf'
fontdir=args.font_dir
for name,filename in [('Body','DejaVuSans.ttf'),('BodyBold','DejaVuSans-Bold.ttf'),('BodyItalic','DejaVuSans-Oblique.ttf'),('BodyBoldItalic','DejaVuSans-BoldOblique.ttf')]:
 pdfmetrics.registerFont(TTFont(name,str(fontdir/filename)))
pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='BodyItalic',boldItalic='BodyBoldItalic')
ink=colors.HexColor('#172F41'); blue=colors.HexColor('#17638B'); grey=colors.HexColor('#5D6B75')
W,H=A4; margin=43; usable=W-2*margin
styles={
 'body':ParagraphStyle('body',fontName='Body',fontSize=9.2,leading=12.8,spaceAfter=6.5,textColor=ink,splitLongWords=True),
 'title':ParagraphStyle('title',fontName='BodyBold',fontSize=20,leading=25,spaceAfter=14,textColor=ink),
 'h2':ParagraphStyle('h2',fontName='BodyBold',fontSize=13.6,leading=18,spaceBefore=13,spaceAfter=8,textColor=ink,keepWithNext=True),
 'h3':ParagraphStyle('h3',fontName='BodyBold',fontSize=10.5,leading=14.5,spaceBefore=10,spaceAfter=6,textColor=blue,keepWithNext=True),
 'bullet':ParagraphStyle('bullet',fontName='Body',fontSize=9.2,leading=12.8,spaceAfter=5.8,leftIndent=15,bulletIndent=1,bulletFontName='Body',bulletFontSize=9.2,textColor=ink,splitLongWords=True),
 'cell':ParagraphStyle('cell',fontName='Body',fontSize=8,leading=11.4,textColor=ink,spaceAfter=0,splitLongWords=True),
 'thead':ParagraphStyle('thead',fontName='BodyBold',fontSize=8.1,leading=11.6,textColor=colors.white,spaceAfter=0,splitLongWords=True),
}
for style in styles.values():
 style.allowWidows=0
 style.allowOrphans=0

def clean(s):
 return s.replace('\u2013','-').replace('\u2014','-').replace('\u2212','-').replace('\u2011','-')
def inline(s,header=False):
 s=clean(s)
 if header: s=s.replace('Null-bootstrap p','Null bootstrap p')
 pattern=re.compile(r'`([^`]+)`|\*\*([^*]+)\*\*|\[([^\]]+)\]\(([^)]+)\)')
 result=[]; last=0
 for m in pattern.finditer(s):
  result.append(html.escape(s[last:m.start()]))
  if m.group(1) is not None:
   result.append('<font color="'+('#FFFFFF' if header else '#374B59')+'">'+html.escape(m.group(1))+'</font>')
  elif m.group(2) is not None:
   result.append('<b>'+inline(m.group(2),header)+'</b>')
  else:
   target=re.sub(r':\d+$','',m.group(4))
   uri=target if target.startswith(('http://','https://')) else (src.parent / target).resolve().as_uri()
   result.append('<a href="'+html.escape(uri,quote=True)+'" color="'+('#FFFFFF' if header else '#17638B')+'">'+html.escape(m.group(3))+'</a>')
  last=m.end()
 result.append(html.escape(s[last:]))
 return ''.join(result)

class NumberedCanvas(canvas.Canvas):
 def __init__(self,*args,**kwargs):
  canvas.Canvas.__init__(self,*args,**kwargs); self.states=[]
 def showPage(self):
  self.states.append(dict(self.__dict__)); self._startPage()
 def save(self):
  total=len(self.states)
  for state in self.states:
   self.__dict__.update(state); self.footer(total); canvas.Canvas.showPage(self)
  canvas.Canvas.save(self)
 def footer(self,total):
  self.saveState(); self.setStrokeColor(colors.HexColor('#D8E2E8')); self.setLineWidth(.5)
  self.line(margin,35,W-margin,35)
  self.setFont('Body',7.3); self.setFillColor(grey)
  self.drawString(margin,23,'WP2 lipid overlap and pooled Figure 3 | 7 October 2026')
  self.drawRightString(W-margin,23,f'Report {self._pageNumber} / {total}')
  self.restoreState()

def colwidths(header):
 n=len(header)
 if n==7 and header[0]=='Sensitivity': return [usable*x for x in [.21,.15,.08,.11,.15,.15,.15]]
 if n==7 and 'Nominal ULM p' in header: return [usable*x for x in [.21,.09,.10,.10,.16,.17,.17]]
 if n==7 and 'Lipid' in header and header[0]!='Lipid': return [usable*x for x in [.15,.09,.27,.10,.10,.10,.19]]
 if n==7 and header[0]=='Lipid': return [usable*x for x in [.30,.115,.115,.115,.115,.12,.12]]
 if n==6 and 'Matched lipids' in header: return [usable*x for x in [.265,.115,.10,.15,.185,.185]]
 if n==6: return [usable*x for x in [.25,.15,.15,.15,.15,.15]]
 if n==8: return [usable*x for x in [.19,.08,.08,.09,.14,.14,.14,.14]]
 if n==7: return [usable*x for x in [.22,.13,.13,.13,.13,.13,.13]]
 if n==5: return [usable*x for x in [.27,.19,.17,.20,.17]]
 if n==4: return [usable*.40,usable*.20,usable*.20,usable*.20]
 if n==4 and header[0]=='Priority': return [usable*.16,usable*.45,usable*.25,usable*.14]
 if n==3 and header[0]=='Dataset': return [usable*.16,usable*.36,usable*.48]
 if n==3 and header[0]=='Cohort': return [usable*.26,usable*.53,usable*.21]
 if n==3 and header[0]=='Action or information': return [usable*.4,usable*.27,usable*.33]
 if n==3: return [usable*.25,usable*.42,usable*.33]
 if n==2 and header[0] in ['Existing analysis','Deliverable']: return [usable*.28,usable*.72]
 if n==2 and header[0]=='Recording interval': return [usable*.2,usable*.8]
 return [usable/n]*n

story=[]; lines=src.read_text().splitlines(); i=0
while i<len(lines):
 line=lines[i]
 if not line.strip(): i+=1; continue
 if line.startswith('```'):
  code=[]; i+=1
  while i<len(lines) and not lines[i].startswith('```'):
   code.append(lines[i]); i+=1
  i+=1
  story.append(Preformatted('\n'.join(code),ParagraphStyle('code',fontName='Courier',fontSize=7.0,leading=10,spaceBefore=4,spaceAfter=8),maxLineLength=112))
  continue
 if line.startswith('|'):
  block=[]
  while i<len(lines) and lines[i].startswith('|'):
   row=[x.strip() for x in lines[i].split('|')[1:-1]]
   if not all(re.fullmatch(r':?-+:?',x) for x in row): block.append(row)
   i+=1
  header=block[0]; cells=[]
  for k,row in enumerate(block): cells.append([Paragraph(inline(x,k==0),styles['thead' if k==0 else 'cell']) for x in row])
  t=Table(cells,colWidths=colwidths(header),repeatRows=1,hAlign='LEFT')
  t.setStyle(TableStyle([
   ('BACKGROUND',(0,0),(-1,0),ink),('VALIGN',(0,0),(-1,-1),'TOP'),
   ('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),
   ('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),
   ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#F1F5F7'),colors.white]),
   ('LINEBELOW',(0,0),(-1,0),.6,ink),('LINEBELOW',(0,1),(-1,-1),.35,colors.HexColor('#D8E2E8')),
  ]))
  if len(block)<=8:
   table_block=[t,Spacer(1,8)]
   while story and isinstance(story[-1],Paragraph) and getattr(story[-1].style,'keepWithNext',False):
    table_block.insert(0,story.pop())
   story.append(KeepTogether(table_block))
  else: story.append(t)
  continue
 if line.startswith('# '): story.append(Paragraph(inline(line[2:]),styles['title'])); i+=1; continue
 if line.startswith('## '): story.append(Paragraph(inline(line[3:]),styles['h2'])); i+=1; continue
 if line.startswith('### '): story.append(Paragraph(inline(line[4:]),styles['h3'])); i+=1; continue
 m=re.match(r'^(- \[ \]|- |\d+\. )(.+)$',line)
 if m:
  marker=m.group(1); bullet='\u2610' if marker=='- [ ]' else ('\u2022' if marker=='- ' else marker.strip())
  story.append(Paragraph(inline(m.group(2)),styles['bullet'],bulletText=bullet)); i+=1; continue
 para=[line]; i+=1
 while i<len(lines) and lines[i].strip() and not re.match(r'^(#|\||- |\d+\. )',lines[i]): para.append(lines[i]); i+=1
 paragraph_style=styles['body']
 next_content=i
 while next_content<len(lines) and not lines[next_content].strip(): next_content+=1
 if next_content<len(lines) and lines[next_content].startswith('|'):
  paragraph_style=ParagraphStyle('before_table',parent=styles['body'],keepWithNext=True)
 story.append(Paragraph(inline(' '.join(para)),paragraph_style))

doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=margin,leftMargin=margin,topMargin=42,bottomMargin=47,
 title='WP2 lipid overlap and pooled Figure 3 conservation',author='FCGL analysis',subject='Lipid annotation audit, frozen pooled Figure 3 pairs and whole-mouse null checks',allowSplitting=True)
story.append(Paragraph('Figure appendix',styles['h2']))
story.append(Paragraph('The following pages show the molecular-overlap audit, frozen Figure 3 coefficient comparisons, whole-mouse null distributions, predefined Cyp summaries and individual SD animals. Pooled replication and excess alignment beyond shared sex structure are separate questions. Nominal ULM p-values are descriptive under dependent edges.',styles['body']))
doc.build(story,canvasmaker=NumberedCanvas)
print(out)

from pypdf import PdfReader, PdfWriter
writer=PdfWriter()
appendices=['lipid_overlap_audit.pdf','pooled_conservation.pdf','pooled_null_checks.pdf','Cyp_pooled_conservation.pdf','selected_pooled_pairs.pdf']
for file in [out]+[src.parent/'plots'/name for name in appendices]:
 start_page=len(writer.pages)
 if not file.exists(): raise FileNotFoundError(file)
 reader=PdfReader(file)
 for page in reader.pages: writer.add_page(page)
 writer.add_outline_item('Report' if file==out else file.stem.replace('_',' '),start_page)
writer.add_metadata({'/Title':'WP2 lipid overlap and pooled Figure 3 conservation','/Author':'FCGL analysis','/Subject':'Lipid overlap, pooled replication and conditional sex-stratified whole-mouse checks with vector figures'})
with final.open('wb') as f: writer.write(f)
print(final)
work.cleanup()
