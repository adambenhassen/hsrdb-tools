"""CIFs whose atoms are read as elements that their own chemical formula does not contain, e.g. a type symbol
`CE` (a CH2 group in a Rietveld refinement) read as cerium.

python dev/element_formula_scan.py LIST [-j N] > element_formula.tsv
One row per CIF with such atoms: the elements, their type symbols and labels, and the formula. CIFs without a
formula are listed when a type symbol is an element symbol written in capitals (`CE`, `CO`), which cannot be
checked against a formula.
"""

import argparse
import math
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi

SYMBOL = re.compile(r"([A-Z][a-z]?)")


def formula_elements(text):
    return {s for s in SYMBOL.findall(text) if gemmi.Element(s).atomic_number > 0}


def scan(path):
    try:
        block = gemmi.cif.read(path).sole_block()
        st = gemmi.make_small_structure_from_block(block)
    except Exception:
        return None
    formula = gemmi.cif.as_string(block.find_value("_chemical_formula_sum") or "").strip()
    wanted = formula_elements(formula) if formula not in ("", "?", ".") else None
    odd = {}
    for s in st.sites:
        if s.element == gemmi.Element("X") or not all(map(math.isfinite, (s.fract.x, s.fract.y, s.fract.z))):
            continue
        el = s.element.name
        if wanted is not None:
            if el not in wanted and not (el == "D" and "H" in wanted):
                odd.setdefault(el, set()).add((s.type_symbol, s.label))
        elif re.fullmatch(r"[A-Z]{2}[0-9+-]*", s.type_symbol or ""):
            odd.setdefault(el, set()).add((s.type_symbol, s.label))
    if not odd:
        return None
    kind = "not in formula" if wanted is not None else "capitals, no formula"
    details = "; ".join(f"{el}: " + ", ".join(sorted(f"{t or '-'}/{lab}" for t, lab in v)[:4]) for el, v in odd.items())
    return f"{path}\t{kind}\t{','.join(sorted(odd))}\t{details}\t{formula}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("list")
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()
    print("path\tkind\telements\ttypes/labels\tformula")
    with ProcessPoolExecutor(args.j) as ex:
        for row in ex.map(scan, Path(args.list).read_text().split(), chunksize=256):
            if row:
                print(row, flush=True)


if __name__ == "__main__":
    main()
