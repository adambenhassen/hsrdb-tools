"""Check every entry of a built .hsrdb: pattern blobs, metadata sanity and structure records.

python dev/audit_database.py DB.hsrdb [workers]
Prints violation counts per check with example reference codes. Exit status 1 if anything failed.
"""

import collections
import math
import multiprocessing
import re
import sqlite3
import struct
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import phase  # noqa: E402

DB = sys.argv[1]
CS_CODES = {"A": {1}, "M": {3}, "O": {5}, "T": {6}, "H": {8, 10}, "C": {9}}
CS_NAMES = {"A": "triclinic", "M": "monoclinic", "O": "orthorhombic", "T": "tetragonal", "C": "cubic"}
ELEMENT = re.compile(r"([A-Z][a-z]?)(\d+\.\d\d)")

Q = """select g.id, g.ProductID, g.XTSLSYS, g.A, g.B, g.C, g.ALPHA, g.BETA, g.GAMMA, g.DX, g.IIC,
              g.NumberOfElements, g.XTLASPECTN, g.XTLVOL, s.SPGR, s.Lines, s.HKL, {lines_i}, t.chemicalformula,
              t.compoundname, e.e1, e.e2, sl.Lines
       from general_indexed g join general_stored s on s.pid=g.id join general_text t on t.docid=g.id
       join pattern_elements e on e.pid=g.id join pattern_strongestlines sl on sl.pid=g.id
       where g.id between ? and ?"""


def check_range(bounds):
    lo, hi = bounds
    db = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    has_i = "LinesI" in [r[1] for r in db.execute("pragma table_info(general_stored)")]
    bad = collections.defaultdict(list)
    n = 0

    def fail(check, code):
        bad[check].append(code)

    for row in db.execute(Q.format(lines_i="s.LinesI" if has_i else "NULL"), (lo, hi)):
        (pid, code, cs, a, b, c, al, be, ga, dx, iic, nel, sgn, vol, spgr, L, H, I, formula, name,
         e1, e2, strongest) = row
        n += 1
        # pattern blobs
        nl = len(L) // 6
        if not 1 <= nl <= 203 or len(L) != 6 * nl or len(H or b"") != 6 * nl:
            fail("line blob sizes", code)
            continue
        lines = [struct.unpack_from("<fH", L, 6 * i) for i in range(nl)]
        ds = [d for d, _ in lines]
        if not all(math.isfinite(d) and d > 0 for d in ds):
            fail("d not finite/positive", code)
        u16 = [i for _, i in lines]
        if u16 != sorted(u16, reverse=True):
            fail("lines not strongest first", code)
        if u16[0] != 1000:
            fail("strongest line not 1000", code)
        if has_i:
            if len(I) != 4 * nl:
                fail("LinesI size", code)
            else:
                fi = struct.unpack(f"<{nl}f", I)
                if not all(math.isfinite(x) and x >= 0 for x in fi):
                    fail("LinesI not finite", code)
                elif any(abs(x - i) > 0.5001 for x, i in zip(fi, u16)):
                    fail("LinesI vs Lines mismatch", code)
        if strongest != L[:60]:
            fail("strongest-lines copy", code)
        # metadata
        if not all(math.isfinite(x) and x > 0 for x in (a, b, c, vol, dx)) or not all(0 < x < 180 for x in (al, be, ga)):
            fail("cell/volume/density invalid", code)
        cell = gemmi.UnitCell(a, b, c, al, be, ga)
        if abs(cell.volume - vol) / vol > 1e-3:
            fail("volume inconsistent with cell", code)
        if not (math.isfinite(iic) and iic >= 0):
            fail("I/Ic invalid", code)
        sg = gemmi.find_spacegroup_by_name(spgr)
        if sg is None or sg.number != sgn:
            fail("SPGR unparsable or number mismatch", code)
        elif cs in CS_NAMES and sg.crystal_system_str() != CS_NAMES[cs]:
            fail("crystal system mismatch", code)
        els = ELEMENT.findall(formula or "")
        if not els or len(els) != nel:
            fail("formula/element count", code)
        else:
            m1 = m2 = 0
            for el, _ in els:
                z = gemmi.Element(el).atomic_number
                if z < 64:
                    m1 |= 1 << z
                else:
                    m2 |= 1 << (z - 64)
            m1 = m1 - (1 << 64) if m1 >= 1 << 63 else m1
            if (m1, m2) != (e1, e2):
                fail("element bitmask mismatch", code)
        if not name:
            fail("empty compound name", code)
        if db.execute("select count(*) from SubFileRef where pid=?", (pid,)).fetchone()[0] == 0:
            fail("no subfile", code)
        if db.execute("select count(*) from literature where pid=?", (pid,)).fetchone()[0] != 1:
            fail("literature row count", code)
        # structure record
        pname = f"PDF:{code[:2]}-{code[2:5]}-{code[5:]}"
        link = db.execute("select PhaseName from PatternPhase where ReferenceCode=?", (code,)).fetchone()
        prow = db.execute("select Phase from Phases where Name=?", (pname,)).fetchone()
        if link is None or link[0] != pname or prow is None:
            fail("phase missing or unlinked", code)
            continue
        try:
            payload = phase.decode_blob(prow[0])
            p = phase.decode_payload(payload)
        except Exception:
            fail("phase does not decode", code)
            continue
        if phase.encode_payload(p) != payload:
            fail("phase round trip", code)
        if p.fmt.version != (6 if has_i else 3):
            fail("phase version vs database", code)
        if p.crystal_system not in CS_CODES.get(cs, set()):
            fail("phase crystal system code", code)
        if not p.atoms:
            fail("phase without atoms", code)
        trvs = [p.scale] + p.cell + p.profile + p.corrections + p.harmonics
        for at in p.atoms:
            trvs += at.pos + [at.biso, at.sof] + at.aniso
            if not 0 <= at.sof.value <= 1:
                fail("occupancy outside [0,1]", code)
            if at.multiplicity < 1 or not at.label or not at.element:
                fail("atom label/element/multiplicity", code)
        if not all(math.isfinite(v) for t in trvs for v in t.values):
            fail("phase value not finite", code)
        for t, ref in zip(p.cell, (a, b, c, al, be, ga)):
            if abs(t.value - ref) > 1e-4 * max(1, abs(ref)):
                fail("phase cell vs pattern cell", code)
                break
    return n, bad


def main():
    db = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    lo, hi = db.execute("select min(id), max(id) from general_indexed").fetchone()
    step = 5000
    ranges = [(s, min(s + step - 1, hi)) for s in range(lo, hi + 1, step)]
    total = 0
    bad = collections.defaultdict(list)
    with multiprocessing.Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as pool:
        for n, b in pool.imap_unordered(check_range, ranges):
            total += n
            for k, v in b.items():
                bad[k] += v
    print(f"{DB}: {total} entries checked")
    for k, v in sorted(bad.items(), key=lambda kv: -len(kv[1])):
        print(f"  FAIL {k}: {len(v)} e.g. {v[:5]}")
    if not bad:
        print("  all checks passed")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
