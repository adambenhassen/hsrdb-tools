"""Atoms whose element gemmi reads from a label in a CIF without _atom_site_type_symbol, where the label is not
plainly element + number: e.g. ring centroids `Cg1` (read as C), `Cent` (Ce), `Cp`, or dummies `DUM1`.

python dev/label_element_scan.py LIST [-j N] > label_element.tsv
One row per such atom: label, element read, calc_flag, occupancy, and whether the element is in the formula.
"""

import argparse
import math
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi


def scan(path):
    try:
        block = gemmi.cif.read(path).sole_block()
        if block.find_values("_atom_site_type_symbol"):
            return None
        st = gemmi.make_small_structure_from_block(block)
    except Exception:
        return None
    flags = [gemmi.cif.as_string(v).strip().lower() for v in block.find_values("_atom_site_calc_flag")]
    formula = gemmi.cif.as_string(block.find_value("_chemical_formula_sum") or "").strip()
    out = []
    for row, s in enumerate(st.sites):
        if not all(map(math.isfinite, (s.fract.x, s.fract.y, s.fract.z))):
            continue
        el = s.element.name
        # the element as the label spells it: leading symbol exactly as written
        m = re.match(r"[A-Z][a-z]?", s.label)
        spelled = m.group(0) if m else ""
        if spelled == el and not s.label[len(el):len(el) + 1].islower():
            continue  # O1, Ca2, C12A, H1'
        flag = flags[row] if row < len(flags) else ""
        in_formula = bool(re.search(rf"\b{el}(?![a-z])", formula))
        out.append(f"{path}\t{s.label}\t{el}\t{flag}\t{s.occ:g}\t{int(in_formula)}\t{formula}")
    return "\n".join(out) or None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("list")
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()
    print("path\tlabel\telement\tcalc_flag\tocc\tin_formula\tformula")
    with ProcessPoolExecutor(args.j) as ex:
        for row in ex.map(scan, Path(args.list).read_text().split(), chunksize=256):
            if row:
                print(row, flush=True)


if __name__ == "__main__":
    main()
