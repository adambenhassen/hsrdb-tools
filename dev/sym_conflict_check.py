"""For CIFs whose listed operations (or Hall symbol) and Hermann-Mauguin symbol describe different groups: which
reading agrees with the CIF's own composition and density.

python dev/sym_conflict_check.py SYM_SOURCES.tsv [--pkg DIR] [-j N] > out.tsv
SYM_SOURCES.tsv is dev/ops_fallback_scan.py output; the CIFs checked are those gemmi cannot match whose operations
or Hall symbol form a group while the H-M symbol is readable. --pkg: directory holding the hsrdb_tools package to
use (e.g. an older version, for comparison). Per CIF: space group, unit-cell formula, largest relative deviation
of a non-H element count from _chemical_formula_sum x Z, and density / _exptl_crystal_density_diffrn.
"""

import argparse
import csv
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

FORMULA = re.compile(r"([A-Z][a-z]?)([0-9.]*)")


def init(pkg):
    sys.path.insert(0, pkg)


def check(path):
    from hsrdb_tools import entry
    try:
        e = entry.from_cif(path)
    except Exception as exc:  # skipped (EntryError) or failed: reported as such
        return f"{path}\t{type(exc).__name__}: {exc}"[:300]
    block = e.block
    fs, z = entry._text(block, "_chemical_formula_sum"), entry._number(block, "_cell_formula_units_Z")
    dev = ""
    if fs and z:
        want = {}
        for el, n in FORMULA.findall(fs):
            want[el] = want.get(el, 0) + (float(n) if n else 1.0) * z
        have = dict(zip(e.elements, (float(t[len(el):]) for el, t in zip(e.elements, e.formula.split()))))
        dev = max((abs(have.get(el, 0) - want.get(el, 0)) / max(have.get(el, 0), want.get(el, 0))
                   for el in set(want) | set(have) if el not in ("H", "D") and max(have.get(el, 0), want.get(el, 0)) > 0),
                  default=0.0)
        dev = f"{dev:.3f}"
    measured = entry._number(block, "_exptl_crystal_density_diffrn")
    ratio = f"{e.density / measured:.3f}" if measured else ""
    return f"{path}\tok\t{e.structure.spacegroup.xhm()}\t{e.formula}\t{dev}\t{ratio}\t{getattr(e, 'symmetry_note', '')}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sources")
    ap.add_argument("--pkg", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("-j", type=int, default=1)
    ap.add_argument("--all", action="store_true", help="every CIF of SYM_SOURCES.tsv")
    args = ap.parse_args()
    paths = [r["path"] for r in csv.DictReader(open(args.sources), delimiter="\t")
             if args.all or (r["gemmi"] == "" and (r["S"].isdigit() or (not r["S"] and r["H"].isdigit()))
                             and r["M"].isdigit())]
    print("path\tstatus\tspace_group\tformula\tformula_dev\tdensity_ratio\tnote", flush=True)
    with ProcessPoolExecutor(args.j, initializer=init, initargs=(args.pkg,)) as ex:
        for row in ex.map(check, paths):
            print(row, flush=True)


if __name__ == "__main__":
    main()
