"""Compare calculated patterns with the ones stored in an official COD .hsrdb for a sample of entries."""

import math
import random
import sqlite3
import struct
import sys
import urllib.request
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import pattern  # noqa: E402

DB = sys.argv[1] if len(sys.argv) > 1 else "COD24_HS4x.hsrdb"
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
CACHE = Path("scratch/cif_cache")
CACHE.mkdir(parents=True, exist_ok=True)


def stored(db, pid):
    L, _, I = db.execute("select Lines, HKL, LinesI from general_stored where pid=?", (pid,)).fetchone()
    n = len(I) // 4
    ds = [struct.unpack_from("<f", L, 6 * i)[0] for i in range(n)]
    return sorted(zip(ds, struct.unpack(f"<{n}f", I)), key=lambda t: -t[0])


def cif(codid):
    p = CACHE / f"{codid}.cif"
    if not p.exists():
        p.write_bytes(urllib.request.urlopen(f"https://www.crystallography.net/cod/{codid}.cif", timeout=60).read())
    return p


db = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
maxid = db.execute("select max(id) from general_indexed").fetchone()[0]
done = 0
errs = []
failed = 0
while done < N:
    if failed > N:
        sys.exit("too many errors")
    pid = random.randint(1, maxid)
    row = db.execute("select ProductID from general_indexed where id=?", (pid,)).fetchone()
    if not row:
        continue
    codid = int(str(row[0])[2:]) - 1
    try:
        st = gemmi.read_small_structure(str(cif(codid)))
        ref = stored(db, pid)
        tt = 2 * math.degrees(math.asin(min(1, pattern.CU_KA1 / (2 * min(d for d, _ in ref)))))
        calc, _ = pattern.calculate(st, two_theta_max=min(tt + 0.01, 140))
    except Exception as e:  # report and continue with the next sample
        print(codid, "ERROR", e)
        failed += 1
        continue
    done += 1
    # compare summed intensity per d value (lines with identical d may be split differently)
    def groups(pairs):
        g = []
        for d, i in sorted(pairs, key=lambda t: -t[0]):
            if g and abs(g[-1][0] - d) / d < 2e-4:
                g[-1][1] += i
            else:
                g.append([d, i])
        return g
    gref = groups(ref)
    gcalc = groups([(x.d, x.intensity) for x in calc])
    diffs = []
    miss = 0
    for d, i in gref:
        best = min(gcalc, key=lambda x: abs(x[0] - d), default=None)
        if best is None or abs(best[0] - d) / d > 2e-3:
            miss += 1
            continue
        if i >= 50:
            diffs.append(abs(best[1] - i))
    worst = max(diffs, default=0)
    errs.append(worst)
    print(f"{codid} n_ref={len(ref)} n_calc={len(calc)} unmatched={miss} worst_dI(I>=50)={worst:.1f}")
errs.sort()
print("median worst dI", errs[len(errs) // 2], "p90", errs[int(len(errs) * 0.9)])
