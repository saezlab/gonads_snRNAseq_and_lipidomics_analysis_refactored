from pathlib import Path
import argparse, tempfile
import re, html
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, KeepTogether, Preformatted
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas

parser=argparse.ArgumentParser(description="Render the amended meeting research plan as a standalone PDF.")
parser.add_argument('--input',type=Path,default=Path('docs/research_plans/2026-10-05_Western_diet_and_Cyp2c38_analysis_plan.md'))
parser.add_argument('--output',type=Path,default=None)
parser.add_argument('--appendix-bottom-padding',type=float,default=4.0,
                    help='Additional points below each unchanged scientific vector page to keep footer glyph extents inside the displayed page.')
parser.add_argument('--font-dir',type=Path,default=Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype')
args=parser.parse_args()
if not 0 <= args.appendix_bottom_padding <= 72:
 parser.error('--appendix-bottom-padding must be between 0 and 72 points')
src=args.input.resolve(strict=True)
final=(args.output or src.with_suffix('.pdf')).resolve()
final.parent.mkdir(parents=True, exist_ok=True)
work=tempfile.TemporaryDirectory(prefix='wp3_package_pdf_')
out=final
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
 'h4':ParagraphStyle('h4',fontName='BodyBold',fontSize=9.4,leading=13.2,spaceBefore=8,spaceAfter=5,textColor=ink,keepWithNext=True),
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
  self.drawString(margin,23,'FCGL research plan | Meeting 5 October 2026 | Amended 14 October 2026')
  self.drawRightString(W-margin,23,f'{self._pageNumber} / {total}')
  self.restoreState()

def colwidths(header):
 n=len(header)
 if n==4: ratios=[.17,.36,.31,.16]
 elif n==3 and header[0]=='Dataset': ratios=[.17,.36,.47]
 elif n==3: ratios=[.49,.17,.34]
 elif n==2: ratios=[.32,.68]
 else: ratios=[1/n]*n
 return [usable*x for x in ratios]

story=[]; lines=src.read_text().splitlines(); i=0
date_line=next((line for line in lines if line.startswith('Report date: ')), '')
report_date=date_line.split(': ',1)[-1].split('. ',1)[0] if date_line else 'Recorded analysis'
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
 if line.startswith('#### '): story.append(Paragraph(inline(line[5:]),styles['h4'])); i+=1; continue
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
 title='Western diet network comparisons and Cyp2c38 validation: amended research plan',author='FCGL analysis',subject='Meeting 5 October 2026, with requested three-cohort WP3 comparison and corrected cohort ages, 14 October 2026',allowSplitting=True)
doc.build(story,canvasmaker=NumberedCanvas)
print(out)

work.cleanup()
