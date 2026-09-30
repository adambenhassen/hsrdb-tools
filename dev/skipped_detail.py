"""For skipped CIFs present in the official database: what the official entry holds and what the CIF says now."""
import re, sqlite3, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import phase

ref = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
for line in Path(sys.argv[1]).read_text().splitlines()[1:]:
    path, reason = line.split("\t", 1)
    if reason.startswith("note:") or "no valid unit cell" in reason:
        continue
    cod = int(Path(path).stem)
    code = str(int(f"96{cod + 1}"))
    row = ref.execute("select s.SPGR, t.chemicalformula from general_indexed g join general_stored s on s.pid=g.id "
                      "join general_text t on t.docid=g.id where g.ProductID=?", (code,)).fetchone()
    if not row:
        continue
    text = Path(path).read_text(errors="replace")
    date = re.search(r"\$Date: (\d{4}-\d\d-\d\d)", text)
    sg = re.findall(r"^_(?:symmetry_space_group_name_H-M|space_group_name_H-M_alt|space_group_name_Hall|symmetry_space_group_name_Hall)\s+(.+)$", text, re.M)
    print(f"{cod} [{reason.split(': ')[-1]}] CIF date {date.group(1) if date else '?'} | CIF symbols {sg} | COD24: {row[0]!r} {row[1]}")
