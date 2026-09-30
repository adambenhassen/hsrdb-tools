"""Rate of strong-line intensity outliers (from compare_databases.py all) per Laue class."""
import collections, sqlite3, sys
import gemmi

db = sqlite3.connect(f"file:{sys.argv[1]}?mode=ro", uri=True)
out = [l.split("\t") for l in open(sys.argv[1].rsplit(".", 1)[0] + ".compare.tsv").read().splitlines()[1:]]
bad = {c for c, k in out if "strong-line" in k}
laue = {}
tot, hit = collections.Counter(), collections.Counter()
for code, sgn, beta in db.execute("select ProductID, XTLASPECTN, BETA from general_indexed"):
    lc = laue.get(sgn) or laue.setdefault(sgn, gemmi.find_spacegroup_by_number(sgn).laue_str())
    if lc == "2/m" and abs(beta - 90) < 0.5:
        lc = "2/m, beta within 0.5 deg of 90"
    tot[lc] += 1
    hit[lc] += code in bad
print(f"{'Laue class':32} entries  outliers   rate")
for lc, n in sorted(tot.items(), key=lambda kv: -hit[kv[0]] / kv[1]):
    print(f"{lc:32} {n:7d} {hit[lc]:8d}  {100 * hit[lc] / n:5.1f}%")
