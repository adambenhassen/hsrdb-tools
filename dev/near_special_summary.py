"""What CIF authors state about sites with nearby images (dev/near_special.py output), by image distance.

python dev/near_special_summary.py NEAR_SPECIAL.tsv
A site counts as stated special when _atom_site_site_symmetry_order > 1, or _atom_site_symmetry_multiplicity
equals the merged image count (SHELXL writes the site-symmetry order) or the general multiplicity divided by it;
as stated general when that item says 1 or the general multiplicity. Occupancy is compared per distinct position
merged (n_within_0.4 / n_same), as sites.merge_tolerances does.
"""
import collections
import math
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
table = collections.defaultdict(collections.Counter)
for r in rows:
    if r["d_max_0.4"] is None or not r["path"].endswith(".cif"):  # text fields with line breaks split a row
        continue
    d = float(r["d_max_0.4"])
    if d == 0:
        continue
    n, n_ops = int(r["n_within_0.4"]), int(r["n_ops"])
    if n == int(r["n_same"]):
        continue  # nothing to merge beyond the exact site symmetry
    occ_n = float(r["occ_at_position"]) * n / int(r["n_same"])
    band = next(b for b in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.41, math.inf) if d < b)  # groups can span > 0.4 A
    stated = ""
    try:
        if r["cif_order"]:
            k = int(float(r["cif_order"]))
            stated = "special" if k > 1 else "general"
        elif r["cif_mult"]:
            m = int(float(r["cif_mult"]))
            if m in (n, n_ops // n) and m not in (1, n_ops):
                stated = "special"
            elif m in (1, n_ops):
                stated = "general"
    except ValueError:
        pass
    table[(band, "occ*n>1" if occ_n > 1.01 else "occ*n<=1")][stated or "no statement"] += 1
print(f"{'d_max <':>8} {'occupancy':>9} {'special':>8} {'general':>8} {'none':>8}")
for (band, occ), c in sorted(table.items()):
    print(f"{band:>8} {occ:>9} {c['special']:>8} {c['general']:>8} {c['no statement']:>8}")
