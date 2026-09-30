"""Official phase record vs ours: same crystal structure? Compare powder patterns computed from both."""
import sys, sqlite3, random
from pathlib import Path
import gemmi, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry, phase, structure, pattern
db = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
cache = sorted(Path("scratch/cif_cache").glob("*.cif")); random.seed(8); random.shuffle(cache)

def to_small(ph, hm):
    st = gemmi.SmallStructure()
    st.cell = gemmi.UnitCell(*[t.values[0] for t in ph.cell])
    st.spacegroup_hm = gemmi.find_spacegroup_by_name(hm).xhm(); st.determine_and_set_spacegroup("2")
    for a in ph.atoms:
        s = gemmi.SmallStructure.Site()
        s.label = a.label; s.type_symbol = a.element; s.element = gemmi.Element(a.element)
        s.fract = gemmi.Fractional(*[t.values[0] for t in a.pos]); s.occ = a.sof.values[0]
        s.u_iso = a.biso.values[0] / (8 * np.pi**2)
        st.add_site(s)
    return st

def key(st):
    lines, _ = pattern.calculate(st, 50)
    return [(round(x.d, 3), round(x.intensity)) for x in lines if x.intensity > 5]

same = diff = 0
for p in cache:
    if same + diff >= int(sys.argv[1] if len(sys.argv) > 1 else 60): break
    cod = int(p.stem); code = str(cod + 1)
    row = db.execute("select p.Phase, s.SPGR from Phases p join general_indexed g on g.ProductID=? join general_stored s on s.pid=g.id where p.Name=?", (code if False else str(int("96" + code)), f"PDF:96-{code[:3]}-{code[3:]}")).fetchone()
    if not row: continue
    ref = phase.decode_payload(phase.decode_blob(row[0]))
    if "TIDY" in ref.trailer["comment"] or len(ref.atoms) > 60: continue
    try:
        e = entry.from_cif(p)
        mine = structure.build(e)
        a, b = key(to_small(ref, row[1])), key(to_small(mine, e.sg_hm))
    except Exception as exc:
        print(cod, "ERR", exc); continue
    if a == b: same += 1
    else:
        diff += 1; print(cod, row[1], e.sg_hm, a[:4], b[:4])
print("same pattern", same, "different", diff)
