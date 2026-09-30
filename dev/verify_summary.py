"""Summarise dev/verify_cctbx.py output: python dev/verify_summary.py OUT.tsv [N examples]."""
import collections
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
n_ex = int(sys.argv[2]) if len(sys.argv) > 2 else 5
status = collections.Counter(r["status"].split(":")[0] for r in rows)
print("rows", len(rows), dict(status))
for r in rows:
    if r["status"].startswith("error"):
        print("  ERROR", r["path"], r["status"])
ok = [r for r in rows if r["status"] == "ok"]


def dist(name, key, rows, limits):
    vals = [(float(r[key]), r) for r in rows if r[key] not in ("", "missing")]
    vals.sort(key=lambda v: v[0])
    print(f"{name}: n={len(vals)}", " ".join(f"<={lim:g}:{sum(v <= lim for v, _ in vals)}" for lim in limits),
          f"max={vals[-1][0]:g}" if vals else "")
    return vals


def worst(vals, k=n_ex, cols=()):
    for v, r in vals[::-1][:k]:
        print(f"    {v:.4g} {r['cod_id']}", " ".join(f"{c}={r[c]}" for c in cols))


cols = ("special_shift", "split_sites", "occ_merged", "invalid_adps", "dropped_sites", "cif_status", "cif_dI", "setting_changed")
worst(dist("calc_rel (|F|^2, >=1% of max)", "calc_rel", ok, [1e-12, 1e-6, 1e-4, 1e-3, 1e-2]), cols=cols)
worst(dist("calc_abs (/max)", "calc_abs", ok, [1e-12, 1e-6, 1e-4, 1e-3]), cols=cols)
print("calc_mult_bad>0:", sum(r["calc_mult_bad"] not in ("", "0") for r in ok))
cif_status = collections.Counter(r["cif_status"].split(":")[0] for r in ok)
print("cif_status", dict(cif_status))
for r in ok:
    if r["cif_status"] != "ok":
        print("   ", r["cod_id"], r["cif_status"][:150])
cif_ok = [r for r in ok if r["cif_status"] == "ok"]
worst(dist("cif_dI (0-1000)", "cif_dI", cif_ok, [0.001, 0.01, 0.1, 0.5, 1, 5]), k=3 * n_ex, cols=cols + ("cif_miss_ours", "cif_miss_cctbx"))
worst(dist("cif_miss_ours", "cif_miss_ours", cif_ok, [0, 0.5, 5]), cols=cols)
worst(dist("cif_miss_cctbx", "cif_miss_cctbx", cif_ok, [0, 0.5, 5]), cols=cols)
dev = [(abs(float(r["iic_ratio"]) - 1), r) for r in cif_ok]
dev.sort(key=lambda v: v[0])
print("iic |ratio-1|:", " ".join(f"<={lim:g}:{sum(v <= lim for v, _ in dev)}" for lim in (1e-4, 5e-4, 1e-3, 1e-2)))
worst(dev, cols=cols)
dev = [(abs(float(r["density_ratio"]) - 1), r) for r in cif_ok]
dev.sort(key=lambda v: v[0])
print("density |ratio-1|:", " ".join(f"<={lim:g}:{sum(v <= lim for v, _ in dev)}" for lim in (1e-4, 5e-4, 1e-3, 1e-2)))
worst(dev, cols=cols)
worst(dist("db_dI", "db_dI", cif_ok, [0.001, 0.01, 0.5, 1]), cols=cols)
mb = [r for r in ok if r["mult_bad"] not in ("", "0")]
print("mult_bad>0:", len(mb), "of", sum(r["mult_bad"] != "" for r in ok), [r["cod_id"] for r in mb[:10]])
ns = [r for r in ok if r["occ_merged"] not in ("", "0")]
print("occ_merged>0:", len(ns), [r["cod_id"] for r in ns[:10]])
worst(dist("formula_dev (non-H)", "formula_dev", ok, [0.001, 0.01, 0.05, 0.2]), k=2 * n_ex, cols=cols)
