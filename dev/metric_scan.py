"""CIFs whose unit cell does not fit their space group: the largest change of the metric tensor under a
symmetry operation, relative to its largest element (0 for a consistent cell).

python dev/metric_scan.py LIST [-j N] > metric.tsv   (rows with a deviation above 1e-6)
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi
import numpy as np


def scan(path):
    try:
        block = gemmi.cif.read(path).sole_block()
        st = gemmi.make_small_structure_from_block(block)
        if st.spacegroup is None:
            st.determine_and_set_spacegroup("SH12")
        if st.spacegroup is None or not st.cell.volume > 0:
            return None
    except Exception:
        return None
    orth = np.array(st.cell.orth.mat.tolist())
    g = orth.T @ orth
    dev = max(np.abs(r.T @ g @ r - g).max() for r in
              (np.array(op.rot, dtype=float) / gemmi.Op.DEN for op in st.spacegroup.operations().sym_ops))
    dev /= np.abs(g).max()
    if dev > 1e-6:
        c = st.cell
        return f"{path}\t{dev:.3e}\t{st.spacegroup.xhm()}\t{c.a} {c.b} {c.c} {c.alpha} {c.beta} {c.gamma}"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("list")
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()
    print("path\tdeviation\tspace_group\tcell")
    with ProcessPoolExecutor(args.j) as ex:
        for row in ex.map(scan, Path(args.list).read_text().split(), chunksize=256):
            if row:
                print(row, flush=True)


if __name__ == "__main__":
    main()
