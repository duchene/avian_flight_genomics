"""14c_table_s1_generation.py

Adds generation length to Table S2 (sheet "Table S2. Morphometrics") of
../data/Table_S1.xlsx: two columns after the last, "Generation length (years)" and
"Generation length source", from ../output/traits_per_taxon.csv (script 14).

The sheet holds only values and inline strings, so the two cells are appended to each row
of its XML directly. No other part of the workbook is read or rewritten, which keeps the
formulas of the other sheets, and their cached values, intact. Rerunning replaces the two
columns rather than adding them again.
"""
import csv, re, shutil, sys, tempfile, zipfile
from xml.sax.saxutils import escape

XLSX, SHEET = "../data/Table_S1.xlsx", "xl/worksheets/sheet4.xml"
COLS = ("AS", "AT")
HEAD = ("Generation length (years)", "Generation length source")

rows = {r["code"]: r for r in csv.DictReader(open("../output/traits_per_taxon.csv"))}

def label(r):
    if r["gen_source"] == "imputed":
        return ("Imputed with phylo.impute from ln body mass (Bird et al. 2020 value is "
                "a family-level extrapolation)")
    return "Bird et al. 2020, Table S4"

def cell(ref, text=None, num=None):
    if num is not None:
        return f'<c r="{ref}"><v>{num}</v></c>'
    return (f'<c r="{ref}" t="inlineStr"><is><t xml:space="preserve">'
            f'{escape(text)}</t></is></c>')

z = zipfile.ZipFile(XLSX)
xml = z.read(SHEET).decode()
xml = re.sub(r'<c r="(?:%s)\d+"[^>]*>.*?</c>' % "|".join(COLS), "", xml)   # idempotent

def add(m):
    row, n = m.group(0), m.group(1)
    if n == "1":
        new = cell(f"AS{n}", HEAD[0]) + cell(f"AT{n}", HEAD[1])
    else:
        code = re.search(r'<c r="C%s"[^>]*><is><t[^>]*>([^<]*)</t>' % n, row).group(1)
        r = rows[code]
        new = cell(f"AS{n}", num=f'{float(r["gen_years"]):.4f}') + cell(f"AT{n}", label(r))
    return row.replace("</row>", new + "</row>")

xml, k = re.subn(r'<row r="(\d+)"[^>]*>.*?</row>', add, xml, flags=re.S)
assert k == 85, k

tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx").name
with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
    for item in z.infolist():
        out.writestr(item, xml.encode() if item.filename == SHEET else z.read(item.filename))
z.close()
shutil.move(tmp, XLSX)
print("added generation length to", k - 1, "taxa in", XLSX)
