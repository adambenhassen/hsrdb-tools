"""Build phase records for entries already in the official database and diff them field by field."""

import collections
import random
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry, phase, structure  # noqa: E402

db = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
cache = sorted(Path("scratch/cif_cache").glob("*.cif"))
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
random.shuffle(cache)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
bad = collections.Counter()
ex = {}
n = 0


def diff(key, want, got, tol=None):
    same = abs(want - got) <= tol if tol is not None else want == got
    if not same:
        bad[key] += 1
        ex.setdefault(key, (cod, want, got))


for p in cache:
    if n >= N:
        break
    cod = int(p.stem)
    code = str(cod + 1)
    row = db.execute("select Phase from Phases where Name=?", (f"PDF:96-{code[:3]}-{code[3:]}",)).fetchone()
    if not row:
        continue
    ref = phase.decode_payload(phase.decode_blob(row[0]))
    if "TIDY" in ref.trailer["comment"]:
        continue  # the official entry was re-standardised; atoms and cell are not comparable
    try:
        e = entry.from_cif(p)
        mine = structure.build(e)
        payload = phase.encode_payload(mine)
        assert phase.encode_payload(phase.decode_payload(payload)) == payload
    except Exception as exc:
        bad["ERROR"] += 1
        ex.setdefault("ERROR", (cod, repr(exc)))
        continue
    n += 1
    diff("prefix", ref.prefix, mine.prefix)
    diff("crystal_system", ref.crystal_system, mine.crystal_system)
    for a, b in zip(ref.cell, mine.cell):
        diff("cell", a.values[0], b.values[0], 1e-3)
        diff("cell_esd", a.values[1], b.values[1], 1e-6)
    diff("natoms", len(ref.atoms), len(mine.atoms))
    if len(ref.atoms) == len(mine.atoms):
        for a, b in zip(ref.atoms, mine.atoms):
            diff("label", a.label, b.label)
            diff("element", a.element, b.element)
            diff("charge", a.charge, b.charge)
            diff("mult", a.multiplicity, b.multiplicity)
            diff("wyckoff", a.wyckoff, b.wyckoff)
            diff("atom_flag", a.tail, b.tail)
            for t1, t2 in zip(a.pos, b.pos):
                diff("xyz", t1.values[0], t2.values[0], 1e-4)
                diff("xyz_esd", t1.values[1], t2.values[1], 1e-6)
            diff("biso", a.biso.values[0], b.biso.values[0], 1e-3)
            diff("biso_esd", a.biso.values[1], b.biso.values[1], 1e-4)
            diff("sof", a.sof.values[0], b.sof.values[0], 1e-4)
            for t1, t2 in zip(a.aniso, b.aniso):
                diff("aniso", t1.values[0], t2.values[0], 1e-3)
                diff("aniso_esd", t1.values[1], t2.values[1], 1e-4)
    for k in ref.trailer:
        if k != "rest":
            diff("trailer." + k, ref.trailer[k], mine.trailer[k])
    rest_a, rest_b = bytearray(ref.trailer["rest"]), bytearray(mine.trailer["rest"])
    rest_a[61:65] = rest_b[61:65] = bytes(4)  # display colour
    diff("trailer.rest", bytes(rest_a), bytes(rest_b))
    for part in ("after_atoms", "mid", "harmonics_head"):
        diff(part, getattr(ref, part), getattr(mine, part))
    diff("harmonics", [(t.name, t.values) for t in ref.harmonics], [(t.name, t.values) for t in mine.harmonics])
    diff("profile", [(t.name, t.values, t.flags) for t in ref.profile + ref.corrections + [ref.scale]],
         [(t.name, t.values, t.flags) for t in mine.profile + mine.corrections + [mine.scale]])
print("compared", n)
for k, v in bad.most_common():
    print(f"{k}: {v}  e.g. {str(ex[k])[:300]}")
