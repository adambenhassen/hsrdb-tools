"""Line-by-line comparison of our pattern and COD24's for outlier entries whose CIF has not changed."""
import random, re, sqlite3, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_databases import groups  # noqa: E402

ours = sqlite3.connect(f"file:{sys.argv[1]}?mode=ro", uri=True)
ref = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
out = Path(sys.argv[1]).with_suffix(".compare.tsv")
codes = sorted({l.split("\t")[0] for l in out.read_text().splitlines()[1:] if "strong-line difference" in l})
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
Q = "select s.Lines, g.XTLASPECTN, s.SPGR from general_indexed g join general_stored s on s.pid=g.id where g.ProductID=?"
shown = 0
for code in random.sample(codes, 200):
    cod = int(code[2:]) - 1
    f = Path(f"cod-mirror/cif/{str(cod)[0]}/{str(cod)[1:3]}/{str(cod)[3:5]}/{cod}.cif")
    text = f.read_text(errors="replace")
    m = re.search(r"\$Date: (\d{4}-\d\d-\d\d)", text)
    if m and m.group(1) > "2024-05-05":
        continue
    a, b = ours.execute(Q, (code,)).fetchone(), ref.execute(Q, (code,)).fetchone()
    go, gr = groups(a[0]), groups(b[0])
    print(f"== COD {cod} sg {a[2]} (official {b[2]}), CIF date {m.group(1) if m else '?'}, aniso={'_atom_site_aniso_U_11' in text}")
    for d, i in gr[:14]:
        mm = min(go, key=lambda x: abs(x[0] - d))
        flag = " <--" if abs(mm[1] - i) > 100 else ""
        print(f"   d {d:8.4f}  official {i:7.1f}  ours {mm[1]:7.1f}{flag}")
    shown += 1
    if shown == 3:
        break
