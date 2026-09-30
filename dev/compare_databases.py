"""Compare a built database with an official one: coverage, and pattern/metadata agreement on shared entries.

python dev/compare_databases.py OURS.hsrdb OFFICIAL.hsrdb [samples|all] [workers]
With "all", every shared entry is compared and outliers are written to OURS.compare.tsv.
"""

import collections
import multiprocessing
import random
import sqlite3
import struct
import sys

OURS, REF = (sys.argv[1:3] + [None, None])[:2]  # also importable for its helpers
Q = ("select g.XTLASPECTN, g.IIC, e.e1, e.e2, s.Lines, g.XTLVOL from general_indexed g "
     "join general_stored s on s.pid=g.id join pattern_elements e on e.pid=g.id where g.ProductID=?")


def groups(blob):
    """Lines summed per d-spacing, scaled so the strongest is 1000 (older databases use other maxima)."""
    lines = sorted((struct.unpack_from("<fH", blob, 6 * i) for i in range(len(blob) // 6)), key=lambda t: -t[0])
    top = max(i for _, i in lines) or 1
    g = []
    for d, i in lines:
        i = 1000 * i / top
        if g and abs(g[-1][0] - d) / d < 2e-4:
            g[-1][1] += i
        else:
            g.append([d, i])
    return g


def compare(codes):
    ours = sqlite3.connect(f"file:{OURS}?mode=ro", uri=True)
    ref = sqlite3.connect(f"file:{REF}?mode=ro", uri=True)
    stats = collections.Counter()
    worst = []
    outliers = []
    for code in codes:
        a, b = ours.execute(Q, (code,)).fetchone(), ref.execute(Q, (code,)).fetchone()
        checks = {
            "space group same": a[0] == b[0],
            "volume within 0.5%": b[5] is None or abs(a[5] - b[5]) / b[5] < 5e-3,  # older DBs lack some
            "elements same": (a[2], a[3]) == (b[2], b[3]),
            "I/Ic within 10%": b[1] is None or abs(a[1] - b[1]) <= 0.1 * max(b[1], 1e-9),
        }
        go, gr = groups(a[4]), groups(b[4])
        dmin = min(d for d, _ in go)
        strong = [(d, i) for d, i in gr if i >= 100 and d >= dmin]
        diffs = []
        found = 0
        for d, i in strong:
            m = min(go, key=lambda x: abs(x[0] - d))
            if abs(m[0] - d) / d < 2e-3:
                found += 1
                diffs.append(abs(m[1] - i))
        checks["strong lines all found"] = found == len(strong)
        w = max(diffs, default=0)
        checks["largest strong-line difference <= 150"] = w <= 150
        for k, ok in checks.items():
            stats[k] += ok
            if not ok:
                outliers.append((code, k))
        worst.append(w)
    return len(codes), stats, worst, outliers


def main():
    ours = sqlite3.connect(f"file:{OURS}?mode=ro", uri=True)
    ref = sqlite3.connect(f"file:{REF}?mode=ro", uri=True)
    codes_ours = {c for (c,) in ours.execute("select ProductID from general_indexed")}
    codes_ref = {c for (c,) in ref.execute("select ProductID from general_indexed")}
    shared = sorted(codes_ours & codes_ref)
    print(f"ours {len(codes_ours)}, official {len(codes_ref)}, shared {len(shared)}, "
          f"only official {len(codes_ref - codes_ours)}, only ours {len(codes_ours - codes_ref)}")
    mode = sys.argv[3] if len(sys.argv) > 3 else "2000"
    if mode != "all":
        random.seed(1)
        shared = random.sample(shared, min(int(mode), len(shared)))
    chunks = [shared[i:i + 2000] for i in range(0, len(shared), 2000)]
    total = 0
    stats = collections.Counter()
    worst = []
    outliers = []
    with multiprocessing.Pool(int(sys.argv[4]) if len(sys.argv) > 4 else 4) as pool:
        for n, s, w, o in pool.imap_unordered(compare, chunks):
            total += n
            stats.update(s)
            worst += w
            outliers += o
    for k, v in stats.items():
        print(f"{k}: {v}/{total} ({100 * v / total:.2f}%)")
    worst.sort()
    print(f"largest strong-line intensity difference per entry: median {worst[total // 2]:.1f}, "
          f"90th percentile {worst[int(total * 0.9)]:.1f}, 99th {worst[int(total * 0.99)]:.1f} (scale 0-1000)")
    if mode == "all":
        with open(OURS.rsplit(".", 1)[0] + ".compare.tsv", "w") as f:
            f.write("reference_code\tfailed_check\n")
            f.writelines(f"{c}\t{k}\n" for c, k in sorted(outliers))
        print(f"{len(outliers)} outliers written to {OURS.rsplit('.', 1)[0]}.compare.tsv")


if __name__ == "__main__":
    main()
