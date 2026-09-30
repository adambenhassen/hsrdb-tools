"""Assign a cause to every discrepancy in dev/verify_cctbx.py output, or list it as unexplained.

python dev/verify_explain.py OUT.tsv... [--limit dI] [--metric metric_all.tsv] [--list CAUSE]

Discrepancies, by kind:
  calc  same structure in both programs: calc_rel or calc_abs > 1e-6, or a multiplicity that differs. Nothing
        explains these; every one is a bug in hsrdb-tools or in the check.
  cif   cctbx's own reading of the CIF: a matched line differs by more than --limit (0-1000 scale, default 0.01), an
        unmatched line is stronger than --limit, or I/Ic or density differ by more than 1e-6 (relative).
  db    the stored lines differ from cctbx's by more than --limit (only with verify_cctbx.py --db).
  mult  the CIF states a site-symmetry order (or multiplicity) other than ours. Listed, not explained.
Known causes of cif discrepancies, where cctbx reads a CIF differently by design or by limitation; first match wins:
  cctbx_error       cctbx could not read the CIF
  modulated         a modulated structure: the stored pattern has satellites and modulated main lines, cctbx
                    calculates the basic structure (checked against the CIF's own structure factors by
                    dev/verify_modulated.py; the calc check still covers the basic structure)
  corrected         built with a correction from the publication (corrections.json)
  choices_not_applied  cctbx's reading has not one scatterer per atom_site row, so dummy atoms, invalid tensors,
                    the default B and merged images could not be applied to it (verify_cctbx.our_choices)
  cell_constrained  cell does not fit the space group and the setting changed: constrained here as in COD24,
                    published cell in cctbx (with an unchanged setting, the check constrains cctbx's cell too)
  reread            atom types that are not elements, or two capitals read again (cctbx reads the CIF's spelling)
  dropped_rows      atom_site rows without coordinates, left out here
  aniso_beta        beta_ij tensors with repeated labels (cctbx reads only U_ij and B_ij; the check adds beta_ij
                    for labels that occur once)
Every cause is printed with its largest difference; --list CAUSE prints all rows of one cause, largest first.
"""

import argparse
import collections
import csv
import re
from pathlib import Path

BETA_TAG = re.compile(r"^\s*_atom_site_aniso_beta_11\b", re.M)


CONSTRAINED = set()  # CIFs whose cell does not fit their space group (dev/metric_scan.py output, --metric)


def positive(r, key):
    return r.get(key) not in ("", "0", None)


def cause(r):
    if r["cif_status"] != "ok":
        return "cctbx_error"
    if positive(r, "modulated"):
        return "modulated"
    if positive(r, "corrected"):
        return "corrected"
    if r.get("choices") == "0":
        return "choices_not_applied"
    if r["path"] in CONSTRAINED and r["setting_changed"] == "1":
        return "cell_constrained"
    if positive(r, "reread"):
        return "reread"
    if positive(r, "dropped_sites"):
        return "dropped_rows"
    text = Path(r["path"]).read_text(errors="replace")
    if BETA_TAG.search(text):
        return "aniso_beta"
    return "unexplained"


def num(r, key, default=0.0):
    try:
        return float(r[key])
    except (KeyError, ValueError):
        return default


def kinds(r, limit):
    """Discrepancies of a row, as {kind: size}."""
    out = {}
    calc = max(num(r, "calc_rel"), num(r, "calc_abs"))
    if calc > 1e-6 or positive(r, "calc_mult_bad"):
        out["calc"] = calc
    if r["cif_status"] != "ok":
        out["cif"] = 0.0
    else:
        size = max(num(r, "cif_dI"), num(r, "cif_miss_ours"), num(r, "cif_miss_cctbx"))
        if (size > limit or abs(num(r, "iic_ratio", 1.0) - 1) > 1e-6
                or abs(num(r, "density_ratio", 1.0) - 1) > 1e-6):
            out["cif"] = size
    if r.get("db_dI") == "missing" or num(r, "db_dI") > limit:
        out["db"] = num(r, "db_dI", float("inf"))
    if positive(r, "mult_bad"):
        out["mult"] = num(r, "mult_bad")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--limit", type=float, default=0.01)
    ap.add_argument("--metric", help="dev/metric_scan.py output: cells constrained here, published in cctbx")
    ap.add_argument("--list", help="print every row of this cause (e.g. unexplained, calc), largest first")
    args = ap.parse_args()
    if args.metric:
        CONSTRAINED.update(r["path"] for r in csv.DictReader(open(args.metric), delimiter="\t"))
    groups = collections.defaultdict(list)  # (kind, cause) -> [(size, row)]
    total, errors, skipped = 0, [], collections.Counter()
    for f in args.files:
        for r in csv.DictReader(open(f), delimiter="\t"):
            total += 1
            if r["status"].startswith("error"):
                errors.append((r["path"], r["status"]))
                continue
            if r["status"] != "ok":
                skipped[re.sub(r"\d+", "N", r["status"])[:90]] += 1
                continue
            for kind, size in kinds(r, args.limit).items():
                c = cause(r) if kind == "cif" else "unexplained" if kind in ("calc", "db") else "stated_order"
                groups[(kind, c)].append((size, r))
    print(f"rows {total}, ok {total - len(errors) - sum(skipped.values())}, skipped {sum(skipped.values())}, "
          f"errors {len(errors)}")
    for p, s in errors:
        print("  ERROR", p, s)
    for reason, n in skipped.most_common():
        print(f"  skipped {n:6}  {reason}")
    for (kind, c), rows in sorted(groups.items()):
        rows.sort(key=lambda x: -x[0])
        print(f"{kind:4} {c:20} {len(rows):6}  max {rows[0][0]:.3g}  e.g. "
              + " ".join(f"{r['cod_id']}({s:.3g})" for s, r in rows[:6]))
    if args.list:
        for (kind, c), rows in sorted(groups.items()):
            if args.list in (c, kind):
                for s, r in rows:
                    print(kind, c, r["cod_id"], f"{s:.4g}", r["path"], sep="\t")


if __name__ == "__main__":
    main()
