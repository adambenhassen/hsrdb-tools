"""Decode and re-encode phase records from an official .hsrdb; report mismatches and field statistics."""

import collections
import random
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import phase  # noqa: E402

db = sqlite3.connect(f"file:{sys.argv[3] if len(sys.argv) > 3 else 'COD24_HS4x.hsrdb'}?mode=ro", uri=True)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
maxid = db.execute("select max(ID) from Phases").fetchone()[0]
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
ok = bad = 0
stats = collections.defaultdict(collections.Counter)
for _ in range(N):
    pid = random.randint(1, maxid)
    row = db.execute("select Phase from Phases where ID=?", (pid,)).fetchone()
    if not row:
        continue
    payload = phase.decode_blob(row[0])
    try:
        p = phase.decode_payload(payload)
        out = phase.encode_payload(p)
    except Exception as e:
        bad += 1
        print(pid, "decode error:", e)
        continue
    if out != payload:
        bad += 1
        i = next(k for k in range(min(len(out), len(payload))) if out[k] != payload[k]) if out[:len(payload)] != payload[:len(out)] else min(len(out), len(payload))
        print(pid, "mismatch at", i, len(out), len(payload))
        continue
    ok += 1
    stats["after_atoms"][p.after_atoms.hex()] += 1
    stats["mid"][p.mid.hex()] += 1
    stats["harmonics_head_len"][len(p.harmonics_head)] += 1
    stats["profile"][tuple(t.name for t in p.profile)] += 1
    stats["corrections"][tuple(t.name for t in p.corrections)] += 1
    for a in p.atoms:
        stats["atom_tail"][a.tail.hex()] += 1
        stats["charge"][a.charge] += 1
    t = p.trailer
    for k in ("pre", "source", "i0", "s1", "s2", "s3", "s4"):
        stats["trailer_" + k][t[k] if not isinstance(t[k], bytes) else t[k].hex()] += 1
    rest = bytearray(t["rest"])
    rest[61:65] = b"CCCC"  # display colour
    stats["rest"][bytes(rest).hex()] += 1
print(f"round-trip ok {ok}, failed {bad}")
for k, c in stats.items():
    print(f"== {k}: {len(c)} distinct")
    for v, n in c.most_common(4):
        print(f"   {n:6d} {str(v)[:180]}")
