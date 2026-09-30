"""Rebuild entries already in the official database from their CIFs and diff every metadata field."""

import random
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry  # noqa: E402

db = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
db.row_factory = sqlite3.Row
cache = sorted(Path("scratch/cif_cache").glob("*.cif"))
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 11)
random.shuffle(cache)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
bad = Counter()
examples = {}
n = 0
for p in cache:
    if n >= N:
        break
    row = db.execute("select * from general where ProductID=?", (str(int(f"96{int(p.stem) + 1}")),)).fetchone()
    if row is None or row["XTLVOL"] > 20000:
        continue
    try:
        e = entry.from_cif(p)
    except Exception as exc:
        bad["ERROR"] += 1
        examples.setdefault("ERROR", (p.stem, str(exc)))
        continue
    n += 1
    lit = db.execute("select * from literature where pid=?", (row["id"],)).fetchone()
    subs = sorted(r[0] for r in db.execute("select SubFileID from SubFileRef where pid=?", (row["id"],)))
    e1, e2 = entry.element_masks(e.elements)
    checks = {
        "XTSLSYS": (row["XTSLSYS"], e.crystal_system),
        "SPGR": (row["SPGR"], e.sg_hm),
        "XTLASPECTN": (row["XTLASPECTN"], e.sg_number),
        "cell_a": (round(row["A"], 3), round(e.cell[0], 3)),
        "DX~2%": (True, abs(row["DX"] - e.density) / row["DX"] < 0.02),
        "IIC~10%": (True, abs(row["IIC"] - e.iic) / max(row["IIC"], 1e-9) < 0.10),
        "formula": (row["chemicalformula"], e.formula),
        "compound": (row["compoundname"], e.compound_name),
        "mineral": (row["mineralname"], e.mineral_name),
        "common": (row["commonname"], e.compound_name),
        "comment": (row["Comment"].split("Publication title")[-1], e.comment.split("Publication title")[-1]),
        "nelem": (row["NumberOfElements"], len(e.elements)),
        "e1e2": ((row["e1"], row["e2"]), (e1, e2)),
        "subfiles": (subs, sorted(e.subfiles)),
        "lit_volume": (lit["volume"] if lit else None, e.literature["volume"]),
        "lit_pages": (lit["pages"] if lit else None, e.literature["pages"]),
        "lit_year": (lit["year"] if lit else None, e.literature["year"]),
        "lit_journal": (lit["journal"] if lit else None, e.literature["journal"]),
        "lit_authors": (lit["authors"] if lit else None, e.literature["authors"]),
    }
    for k, (want, got) in checks.items():
        if want != got:
            bad[k] += 1
            examples.setdefault(k, (p.stem, want, got))
print(f"checked {n}")
for k, v in bad.most_common():
    print(f"{k}: {v} differ, e.g. {examples[k]}")
