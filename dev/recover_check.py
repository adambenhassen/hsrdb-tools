"""Check CIFs recovered by setting.from_operations: the standardised pattern must equal the pattern of the CIF's
own sites expanded with its listed operations in its own cell.

python dev/recover_check.py UNKNOWN_SG.tsv   (from dev/unknown_sg.py)
"""

import csv
import sys
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry, pattern, setting, sites  # noqa: E402


def direct_pattern(path):
    """Pattern of the CIF expanded with its own operations, as a P 1 structure in its own cell."""
    block = gemmi.cif.read(path).sole_block()
    st = gemmi.make_small_structure_from_block(block)
    ops = setting.listed_operations(st, block)
    st, _, _, _ = entry._clean_sites(st, gemmi.SpaceGroup("P 1"))
    images = sites.unit_cell_images(st, ops)
    p1 = gemmi.SmallStructure()
    p1.cell = st.cell
    p1.spacegroup_hm = "P 1"
    p1.determine_and_set_spacegroup("2")
    orth = np.array(st.cell.orth.mat.tolist())
    for i, x, u in zip(images.site, images.xyz, images.ucart):
        s = st.sites[i].clone()
        s.fract = gemmi.Fractional(*x)
        s.aniso = setting._cart_to_ucif(u, st.cell) if st.sites[i].aniso.nonzero() else gemmi.SMat33d(0, 0, 0, 0, 0, 0)
        p1.add_site(s)
    del orth
    return pattern.calculate(p1, 90)[0]


def compare(a, b):
    """Largest intensity difference of lines matched within 0.006 deg 2-theta, strongest unmatched line."""
    tt = lambda d: 2 * np.degrees(np.arcsin(pattern.CU_KA1 / (2 * d)))  # noqa: E731
    ta, tb = [tt(x.d) for x in a], [tt(x.d) for x in b]
    used, di, unmatched = set(), 0.0, 0.0
    for t, x in zip(ta, a):
        j = min(range(len(tb)), key=lambda k: abs(tb[k] - t), default=None)
        if j is not None and abs(tb[j] - t) <= 0.006:
            used.add(j)
            di = max(di, abs(x.intensity - b[j].intensity))
        else:
            unmatched = max(unmatched, x.intensity)
    unmatched = max([unmatched] + [b[j].intensity for j in range(len(b)) if j not in used])
    return di, unmatched


rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
for r in rows:
    path = r["path"]
    name = Path(path).stem
    symbol = r["_symmetry_space_group_name_H-M"] or r["_space_group_name_H-M_alt"] or r["_space_group_name_Hall"] \
        or r["_symmetry_space_group_name_Hall"]
    try:
        e = entry.from_cif(path)
    except (entry.EntryError, ValueError, RuntimeError) as exc:
        print(f"{name}\tSKIP\t{symbol!r}\tops={r['n_ops']}\t{exc}")
        continue
    ours = pattern.calculate(e.structure, 90)[0]
    if not setting.listed_operations(gemmi.make_small_structure_from_block(gemmi.cif.read(path).sole_block()),
                                     gemmi.cif.read(path).sole_block()):
        print(f"{name}\tOK\t{symbol!r}\t-> {e.sg_hm} ({e.sg_number})\tcell {tuple(round(v, 4) for v in e.cell)}"
              f"\tsymbol only, no direct check\t{e.formula}")
        continue
    ref = direct_pattern(path)
    di, unmatched = compare(ours, ref)
    print(f"{name}\tOK\t{symbol!r}\t-> {e.sg_hm} ({e.sg_number})\tcell {tuple(round(v, 4) for v in e.cell)}"
          f"\tmax dI {di:.2e}, unmatched {unmatched:.2e}\t{e.formula}")
