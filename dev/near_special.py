"""Sites whose symmetry images lie close to each other but do not coincide: on a special position with rounded
coordinates, or a split (disordered) atom next to one.

python dev/near_special.py LIST [-j N] > near_special.tsv
LIST: file with one CIF path per line. One row per such site: the shortest and longest image distance within
0.6 A, the number of images in the site's group within 0.4 A and within 0.01 A (sites._groups: the operations
mapping the site onto the positions sites.unit_cell_images would merge, and its exact site-symmetry order), the
largest distance in that group, the occupancy, the summed
occupancy of all sites at that position, and what the CIF states about the site (symmetry multiplicity,
site-symmetry order, Wyckoff letter).
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry, sites  # noqa: E402

COLUMNS = ["path", "label", "element", "occ", "occ_at_position", "n_ops", "n_within_0.4", "n_same", "n_within_0.6",
           "d_min", "d_max_0.4", "cif_mult", "cif_order", "cif_wyckoff", "formula_sum", "Z"]


def value(block, tag, row):
    col = block.find_values(tag)
    if not col or row >= len(col):
        return ""
    return gemmi.cif.as_string(col[row])


def scan(path):
    out = []
    try:
        block = gemmi.cif.read(path).sole_block()
        st = gemmi.make_small_structure_from_block(block)
        if st.spacegroup is None:
            st.determine_and_set_spacegroup("SH12")
        if st.spacegroup is None or not st.cell.volume > 0:
            return out
        st, rows, _, _ = entry._clean_sites(st, st.spacegroup)
    except Exception:
        return out
    orth = np.array(st.cell.orth.mat.tolist())
    ops = list(st.spacegroup.operations())
    rot = np.array([op.rot for op in ops], dtype=float) / gemmi.Op.DEN
    tran = np.array([op.tran for op in ops], dtype=float) / gemmi.Op.DEN
    metric = orth.T @ orth
    xyz = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in st.sites])
    for i, s in enumerate(st.sites):
        x = (rot @ xyz[i] + tran) % 1.0
        diff = (x[1:] - x[0] + 0.5) % 1.0 - 0.5  # images relative to the site itself
        dist = np.sqrt(np.einsum("ai,ij,aj->a", diff, metric, diff))
        near = dist[(dist > 1e-3) & (dist < 0.6)]
        if not len(near):
            continue
        d_all = (xyz - xyz[i] + 0.5) % 1.0 - 0.5
        colocated = np.sqrt(np.einsum("ai,ij,aj->a", d_all, metric, d_all)) < 0.05
        occ_here = sum(st.sites[k].occ for k in np.flatnonzero(colocated))
        row = rows[i]
        xi, diff_i, dist2 = sites._image_distances(st, i, rot, tran, metric)
        group = sites._groups(xi, diff_i, dist2, metric, sites.SPECIAL_POSITION_TOLERANCE) == 0
        n_same = int(np.count_nonzero(sites._groups(xi, diff_i, dist2, metric, sites.SAME_POSITION) == 0))
        out.append([path, s.label, s.element.name, f"{s.occ:.4f}", f"{occ_here:.4f}", len(ops),
                    int(np.count_nonzero(group)), n_same,
                    1 + int(np.count_nonzero(dist < 0.6)), f"{near.min():.4f}",
                    f"{np.sqrt(dist2[0, group].max()):.4f}",
                    value(block, "_atom_site_symmetry_multiplicity", row),
                    value(block, "_atom_site_site_symmetry_order", row),
                    value(block, "_atom_site_Wyckoff_symbol", row),
                    entry._text(block, "_chemical_formula_sum"), entry._text(block, "_cell_formula_units_Z")])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("list")
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()
    paths = Path(args.list).read_text().split()
    print("\t".join(COLUMNS), flush=True)
    with ProcessPoolExecutor(args.j) as ex:
        for rows in ex.map(scan, paths, chunksize=64):
            for r in rows:
                print("\t".join(map(str, r)), flush=True)


if __name__ == "__main__":
    main()
