#!/usr/bin/env python3
"""Render the reviewed KO report without changing its scientific content."""
from __future__ import annotations
import argparse
import hashlib
import html
import json
from pathlib import Path
import re
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, PageTemplate, Paragraph,
                              Table, TableStyle, Spacer, PageBreak,
                              NextPageTemplate, Preformatted)

ROOT = Path(__file__).resolve().parents[2]
FONTDIR = Path('/Users/dugourd/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype')
INK = colors.HexColor('#172C42')
ACCENT = colors.HexColor('#246C77')
MUTED = colors.HexColor('#586674')
RULE = colors.HexColor('#D9E1E5')


def normalize(s):
    for char in ['\u2010', '\u2011', '\u2012', '\u2013', '\u2014']:
        s = s.replace(char, '-')
    return s


def inline(text):
    text = normalize(text)
    # Escape all input; only this small Markdown subset generates markup.
    parts = re.split(r'(\[[^\]]+\]\([^)]+\)|`[^`]+`|\*\*[^*]+\*\*)', text)
    out = []
    for part in parts:
        if part.startswith('[') and (m := re.fullmatch(r'\[([^\]]+)\]\(([^)]+)\)', part)):
            out.append(f'<link href="{html.escape(m[2], quote=True)}" color="#246C77"><u>{html.escape(m[1])}</u></link>')
        elif part.startswith('`') and part.endswith('`'):
            out.append(f'<font name="Mono" size="8.6">{html.escape(part[1:-1])}</font>')
        elif part.startswith('**') and part.endswith('**'):
            out.append('<b>' + html.escape(part[2:-2]) + '</b>')
        else:
            out.append(html.escape(part))
    return ''.join(out)


def blocks(text):
    lines = text.splitlines(); i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1; continue
        if line.startswith('```'):
            i += 1; code = []
            while i < len(lines) and not lines[i].startswith('```'):
                code.append(lines[i]); i += 1
            i += 1; yield 'code', '\n'.join(code); continue
        if line.startswith('# '):
            yield 'title', line[2:]; i += 1; continue
        if line.startswith('## '):
            yield 'heading', line[3:]; i += 1; continue
        if line.startswith('|'):
            rows = []
            while i < len(lines) and lines[i].startswith('|'):
                cells = [x.strip() for x in lines[i].strip('|').split('|')]
                if not all(re.fullmatch(r':?-+:?', x) for x in cells):
                    rows.append(cells)
                i += 1
            yield 'table', rows; continue
        if line.startswith('- '):
            paragraph = [line[2:]]; i += 1
            while i < len(lines) and lines[i].startswith('  '):
                paragraph.append(lines[i].strip()); i += 1
            yield 'bullet', ' '.join(paragraph); continue
        paragraph = [line.strip()]; i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(('#', '|', '- ', '```')):
            paragraph.append(lines[i].strip()); i += 1
        yield 'paragraph', ' '.join(paragraph)


def page_furniture(canvas, doc):
    w, h = canvas._pagesize
    canvas.saveState()
    canvas.setStrokeColor(RULE); canvas.setLineWidth(.65)
    canvas.line(48, h - 40, w - 48, h - 40)
    canvas.setFont('Body', 8); canvas.setFillColor(MUTED)
    canvas.drawString(48, h - 29, 'FCG liver study | Factor 1 TF-to-lipid validation')
    canvas.drawRightString(w - 48, h - 29, 'Reviewed 13 October 2026')
    canvas.line(48, 35, w - 48, 35)
    canvas.setFont('Body', 7.5)
    canvas.drawString(48, 22, 'Experimental proposal | Conditional mechanistic hypotheses')
    canvas.drawRightString(w - 48, 22, str(doc.page))
    canvas.restoreState()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'docs/cosmos_factor1_ko_validation_targets.md')
    parser.add_argument('--output', type=Path, default=ROOT / 'output/pdf/cosmos_factor1_ko_validation_targets.pdf')
    args = parser.parse_args()
    original = args.source.read_bytes(); text = original.decode()
    for name, fname in [('Body','DejaVuSans.ttf'), ('Body-Bold','DejaVuSans-Bold.ttf'),
                        ('Body-Italic','DejaVuSans-Oblique.ttf'), ('Body-BoldItalic','DejaVuSans-BoldOblique.ttf'),
                        ('Mono','DejaVuSansMono.ttf')]:
        pdfmetrics.registerFont(TTFont(name, str(FONTDIR/fname)))
    pdfmetrics.registerFontFamily('Body', normal='Body', bold='Body-Bold', italic='Body-Italic', boldItalic='Body-BoldItalic')
    styles = {
        'body': ParagraphStyle('BodyText', fontName='Body', fontSize=10, leading=14.3, textColor=INK,
                               spaceAfter=9, splitLongWords=True, allowWidows=0, allowOrphans=0),
        'title': ParagraphStyle('Title', fontName='Body-Bold', fontSize=22, leading=28, textColor=INK, spaceAfter=17),
        'heading': ParagraphStyle('Heading', fontName='Body-Bold', fontSize=13, leading=18, textColor=ACCENT,
                                  spaceBefore=10, spaceAfter=9, keepWithNext=True),
        'bullet': ParagraphStyle('Bullet', fontName='Body', fontSize=10, leading=14.3, textColor=INK,
                                 leftIndent=13, firstLineIndent=0, bulletIndent=0, spaceAfter=8,
                                 bulletFontName='Body', bulletFontSize=10,
                                 splitLongWords=True, allowWidows=0, allowOrphans=0),
        'cell': ParagraphStyle('Cell', fontName='Body', fontSize=9.1, leading=12.1, textColor=INK,
                               splitLongWords=True),
        'rankcell': ParagraphStyle('RankCell', fontName='Body', fontSize=9, leading=11.8, textColor=INK,
                                   splitLongWords=True),
        'header': ParagraphStyle('TableHead', fontName='Body-Bold', fontSize=9.1, leading=12.1,
                                 textColor=colors.white),
        'code': ParagraphStyle('Code', fontName='Mono', fontSize=8.3, leading=12,
                               textColor=INK, backColor=colors.HexColor('#F3F6F8'),
                               borderPadding=9, spaceBefore=3, spaceAfter=12),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    doc = BaseDocTemplate(str(args.output), pagesize=A4, leftMargin=48, rightMargin=48,
                          topMargin=54, bottomMargin=46, title='Ranked knockout and knockdown targets for TF to lipid validation',
                          author='Gonadal hormone coordination of sex-specific adult liver physiology',
                          subject='Reviewed alternating enzyme and TF validation panel with compensation rationale')
    portrait_frame = Frame(48, 46, A4[0]-96, A4[1]-100, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    wide = landscape(A4)
    landscape_frame = Frame(48, 46, wide[0]-96, wide[1]-100, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates([PageTemplate(id='portrait', frames=[portrait_frame], onPage=page_furniture, pagesize=A4),
                          PageTemplate(id='landscape', frames=[landscape_frame], onPage=page_furniture, pagesize=wide)])
    story = []; ranking = False; route_page = False; table_heights = []
    for kind, content in blocks(text):
        if kind == 'title':
            story.append(Paragraph(inline(content), styles['title']))
        elif kind == 'heading':
            if content == 'Revised alternating ranking':
                story.extend([NextPageTemplate('landscape'), PageBreak()]); ranking = True
            elif content == 'Mechanistic route tested by each KO set':
                if ranking:
                    story.append(PageBreak())
                else:
                    story.extend([NextPageTemplate('landscape'), PageBreak()])
                ranking = True; route_page = True
            elif ranking:
                story.extend([NextPageTemplate('portrait'), PageBreak()]); ranking = False; route_page = False
            elif content in ['Interpretation of the TF arm', 'Expression, controls and experimental sequence',
                             'Provenance and reproduction']:
                story.append(PageBreak())
            story.append(Paragraph(inline(content), styles['heading']))
        elif kind in ['paragraph','bullet']:
            style = styles['bullet'] if kind == 'bullet' else styles['body']
            story.append(Paragraph(inline(content), style, bulletText='•' if kind=='bullet' else None))
        elif kind == 'code':
            story.append(Preformatted(normalize(content), styles['code']))
        elif kind == 'table':
            if route_page:
                widths = [105, 355, wide[0]-96-460]
            elif ranking:
                widths = [40, 90, 83, 151, wide[0]-96-364]
            else:
                widths = [49, 259, A4[0]-96-308]
            cellstyle = styles['rankcell'] if ranking else styles['cell']
            data = [[Paragraph(inline(cell), styles['header'] if i==0 else cellstyle)
                     for cell in row] for i,row in enumerate(content)]
            table = Table(data, colWidths=widths, repeatRows=1, hAlign='LEFT', splitByRow=True)
            commands = [('BACKGROUND',(0,0),(-1,0),ACCENT), ('VALIGN',(0,0),(-1,-1),'TOP'),
                        ('LEFTPADDING',(0,0),(-1,-1),7), ('RIGHTPADDING',(0,0),(-1,-1),7),
                        ('TOPPADDING',(0,0),(-1,-1),6), ('BOTTOMPADDING',(0,0),(-1,-1),6),
                        ('LINEBELOW',(0,0),(-1,0),.7,ACCENT),
                        ('LINEBELOW',(0,1),(-1,-1),.5,RULE)]
            if ranking:
                for i,row in enumerate(content[1:],1):
                    commands.append(('BACKGROUND',(0,i),(-1,i),colors.HexColor('#EDF5F4' if i%2==0 else '#F7F9FB')))
            else:
                commands.append(('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#F7F9FB'),colors.white]))
            table.setStyle(TableStyle(commands))
            table_heights.append({'ranking':ranking,'width':sum(widths),'height':table.wrap(sum(widths),1000)[1]})
            story.extend([table,Spacer(1,10)])
    doc.build(story)
    assert args.source.read_bytes() == original, 'Source Markdown was changed'
    print(json.dumps({'output':str(args.output), 'source_sha256':hashlib.sha256(original).hexdigest(),
                      'tables':table_heights,'source_modified':False},indent=2))


if __name__=='__main__':
    main()
