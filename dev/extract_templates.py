"""Derive the constant parts of phase records from an official .hsrdb and write hsrdb_tools/templates.json.

Reports how many distinct variants each part has; anything that is not unique per key needs a closer look.
"""

import collections
import json
import random
import re
import sqlite3
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import phase  # noqa: E402

N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
DB = sys.argv[2] if len(sys.argv) > 2 else "COD24_HS4x.hsrdb"
OUT = sys.argv[3] if len(sys.argv) > 3 else "hsrdb_tools/templates.json"
db = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
random.seed(9)
ids = [r[0] for r in db.execute("select p.ID from Phases p")]
random.shuffle(ids)

trv_kinds = collections.defaultdict(collections.Counter)  # kind -> (refine, flags, limits) counts
laue = collections.defaultdict(collections.Counter)  # laue class -> (head, harmonics) counts
fixed = collections.defaultdict(collections.Counter)


def key_of(t, keep_value):
    v = list(t.values)
    if not keep_value:
        v[0] = v[1] = None  # value and esd come from the CIF
    return json.dumps([t.refine, t.flags.hex(), 0, v])


for pid in ids[:N]:
    blob, name = db.execute("select Phase, Name from Phases where ID=?", (pid,)).fetchone()
    code = name.split(":")[1].replace("-", "")
    sg = db.execute("select g.XTLASPECTN, s.SPGR from general_indexed g join general_stored s on s.pid=g.id "
                    "where g.ProductID=?", (code,)).fetchone()
    p = phase.decode_payload(phase.decode_blob(blob))
    version = p.fmt.version
    profile_names = [t.name for t in p.profile]
    correction_names = [t.name for t in p.corrections]
    trv_kinds["Scale Factor"][key_of(p.scale, True)] += 1
    for t in p.cell:
        trv_kinds[t.name][key_of(t, False)] += 1
    for a in p.atoms:
        for t in a.pos + [a.biso, a.sof] + a.aniso:
            kind = "ATOM " + re.sub(r"^.* ", "", t.name)
            trv_kinds[kind][key_of(t, False)] += 1
    for t in p.profile + p.corrections:
        trv_kinds[t.name][key_of(t, True)] += 1
    group = gemmi.find_spacegroup_by_name(sg[1]) if sg else None
    lc = group.laue_str() if group else "?"
    laue[lc][json.dumps([p.harmonics_head.hex(), [[t.name, key_of(t, True)] for t in p.harmonics]])] += 1
    fixed["after_atoms"][p.after_atoms.hex()] += 1
    fixed["mid"][p.mid.hex()] += 1
    rest = bytearray(p.trailer["rest"])
    rest[24:28] = bytes(4)  # journal issue
    fixed["trailer_rest"][bytes(rest).hex()] += 1

templates = {"trv": {}, "laue": {}, "after_atoms": None, "mid": None, "trailer_rest": None,
             "profile_names": profile_names, "correction_names": correction_names, "record_version": version}
for kind, c in sorted(trv_kinds.items()):
    print(f"{kind}: {len(c)} variants; top {c.most_common(1)[0][1]}/{sum(c.values())}")
    for v, n in c.most_common(3)[1:]:
        print("      also", n, v[:160])
    templates["trv"][kind] = json.loads(c.most_common(1)[0][0])
for lc, c in sorted(laue.items()):
    print(f"laue {lc}: {len(c)} variants; top {c.most_common(1)[0][1]}/{sum(c.values())}")
    head, harm = json.loads(c.most_common(1)[0][0])
    templates["laue"][lc] = {"head": head, "harmonics": [[n, json.loads(k)] for n, k in harm]}
for k in ("after_atoms", "mid", "trailer_rest"):
    print(k, len(fixed[k]), "variants")
    templates[k] = fixed[k].most_common(1)[0][0]
Path(OUT).write_text(json.dumps(templates, indent=1, ensure_ascii=False))
