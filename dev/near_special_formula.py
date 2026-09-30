"""Which reading of near-special sites matches the CIF's own formula?

python dev/near_special_formula.py NEAR_SPECIAL.tsv [-j N] > out.tsv
For each CIF listed by dev/near_special.py with images between 0.001 and 0.4 A apart: the largest relative
difference of a non-H element count in the unit cell from _chemical_formula_sum x Z, for each merge rule in
POLICIES (images closer than the tolerance are one atom; "occ" rules also merge images up to 0.4 A apart when
the occupancy at the position times the number of images exceeds 1).
"""

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from hsrdb_tools import entry  # noqa: E402
from verify_cctbx import FORMULA  # noqa: E402

ROUND = 0.05
POLICIES = {"merge_0.4": (0.4,), "merge_0.05": (0.05,), "merge_0.01": (0.01,), "merge_0.001": (0.001,),
            "occ_0.01": (0.01, True), "occ_0.001": (0.001, True)}


def counts(st, tol, occ_rule=False):
    """Element counts in the cell when images closer than tol are one atom. With occ_rule, images up to 0.4 A
    apart are also one atom when the occupancy at the position times the number of images exceeds 1."""
    orth = np.array(st.cell.orth.mat.tolist())
    metric = orth.T @ orth
    ops = list(st.spacegroup.operations())
    rot = np.array([op.rot for op in ops], dtype=float) / gemmi.Op.DEN
    tran = np.array([op.tran for op in ops], dtype=float) / gemmi.Op.DEN
    xyz = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in st.sites])
    out = {}
    for i, s in enumerate(st.sites):
        x = (rot @ xyz[i] + tran) % 1.0
        diff = (x[:, None, :] - x[None, :, :] + 0.5) % 1.0 - 0.5
        dist2 = np.einsum("abi,ij,abj->ab", diff, metric, diff)
        t = tol
        if occ_rule:
            d = (xyz - xyz[i] + 0.5) % 1.0 - 0.5
            occ_here = sum(st.sites[k].occ for k in np.flatnonzero(np.einsum("ai,ij,aj->a", d, metric, d) < ROUND**2))
            if occ_here * np.count_nonzero(dist2[0] < 0.4**2) > 1.01:
                t = 0.4
        n = int(np.count_nonzero(~np.triu(dist2 < t**2, 1).any(axis=0)))
        out[s.element.name] = out.get(s.element.name, 0.0) + s.occ * n
    return out


def dev(have, want):
    d = 0.0
    for el in set(want) | set(have):
        if el in ("H", "D"):
            continue
        w, h = want.get(el, 0.0), have.get(el, 0.0)
        if max(w, h) > 0:
            d = max(d, abs(h - w) / max(w, h))
    return d


def check(path):
    try:
        block = gemmi.cif.read(path).sole_block()
        st = gemmi.make_small_structure_from_block(block)
        if st.spacegroup is None:
            st.determine_and_set_spacegroup("SH12")
        st, _, _, _ = entry._clean_sites(st, st.spacegroup)
        fs, z = entry._text(block, "_chemical_formula_sum"), entry._number(block, "_cell_formula_units_Z")
        if not fs or not z:
            return [path] + [""] * len(POLICIES) + ["no formula or Z"]
        want = {}
        for el, n in FORMULA.findall(fs):
            want[el] = want.get(el, 0) + (float(n) if n else 1.0) * z
        return [path] + [f"{dev(counts(st, *p), want):.4f}" for p in POLICIES.values()] + [""]
    except Exception as e:
        return [path] + [""] * len(POLICIES) + [f"error {e}"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("near_special")
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()
    rows = list(csv.DictReader(open(args.near_special), delimiter="\t"))
    ambiguous = sorted({r["path"] for r in rows if float(r["d_max_0.4"]) > 0.001})
    print("\t".join(["path", *POLICIES, "note"]))
    with ProcessPoolExecutor(args.j) as ex:
        for r in ex.map(check, ambiguous, chunksize=16):
            print("\t".join(r), flush=True)


if __name__ == "__main__":
    main()
