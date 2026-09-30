"""Why do entries differ from the official database? Sample large intensity outliers and classify them.

python dev/outlier_causes.py OURS.hsrdb [n]   (needs COD24_HS4x.hsrdb and the COD mirror in cod-mirror/cif)
"""
import collections, random, re, sqlite3, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_databases import groups  # noqa: E402

ours = sqlite3.connect(f"file:{sys.argv[1]}?mode=ro", uri=True)
ref = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
out = Path(sys.argv[1]).with_suffix(".compare.tsv")
codes = sorted({l.split("\t")[0] for l in out.read_text().splitlines()[1:] if "strong-line difference" in l})
random.seed(3)
sample = random.sample(codes, min(int(sys.argv[2]) if len(sys.argv) > 2 else 200, len(codes)))
Q = "select s.Lines from general_indexed g join general_stored s on s.pid=g.id where g.ProductID=?"
causes = collections.Counter()
for code in sample:
    cod = int(code[2:]) - 1
    f = Path(f"cod-mirror/cif/{str(cod)[0]}/{str(cod)[1:3]}/{str(cod)[3:5]}/{cod}.cif")
    text = f.read_text(errors="replace")
    m = re.search(r"\$Date: (\d{4}-\d\d-\d\d)", text)
    changed = bool(m and m.group(1) > "2024-05-05")
    aniso = "_atom_site_aniso_U_11" in text or "_atom_site_aniso_B_11" in text
    hydro = bool(re.search(r"^\s*H\w*\s+H\b", text, re.M))
    partial = bool(re.search(r"_atom_site_occupancy", text))
    key = ("CIF changed after May 2024" if changed else "CIF unchanged",
           "aniso" if aniso else "iso")
    causes[key] += 1
for k, v in causes.most_common():
    print(v, k)
