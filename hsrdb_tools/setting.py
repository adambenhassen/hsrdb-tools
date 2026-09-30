"""Move a structure into the space-group setting HighScore assumes for a bare Hermann-Mauguin symbol.

That is gemmi's reference setting: standard axes, origin choice 2 and hexagonal axes for rhombohedral groups,
the convention of the official COD databases (whose atom positions may differ by symmetry-equivalent choices).
"""

import re
from functools import lru_cache

import gemmi
import numpy as np
import spglib

from .sites import ucif_to_cart, unit_cell_images


@lru_cache(maxsize=None)
def reference_spacegroup(number):
    for g in gemmi.spacegroup_table():
        if g.number == number and g.basisop == gemmi.Op("x,y,z"):
            return g
    raise ValueError(f"no reference setting for space group {number}")


def identify_shifted(symops):
    """Space group whose operations match the listed ones after an origin shift s.

    For CIFs whose symbol gemmi does not know (e.g. 'P 2yb (x+1/4,y,z)'). Returns (group, s) with s in
    eighths, or None; coordinates in the group's setting are x - s.
    """
    try:
        ops = [gemmi.Op(o) for o in symops]
    except RuntimeError:
        return None
    den = gemmi.Op.DEN
    rots = {tuple(map(tuple, op.rot)) for op in ops}
    listed = {(tuple(map(tuple, op.rot)), tuple(t % den for t in op.tran)) for op in ops}
    for g in gemmi.spacegroup_table():
        gops = list(gemmi.SpaceGroup(g.xhm()).operations())
        if len(gops) != len(listed) or {tuple(map(tuple, o.rot)) for o in gops} != rots:
            continue
        for sx in range(0, den, den // 8):
            for sy in range(0, den, den // 8):
                for sz in range(0, den, den // 8):
                    s = np.array([sx, sy, sz])
                    shifted = set()
                    for o in gops:
                        r = np.array(o.rot)
                        t = np.array(o.tran) + s - (r @ s) // den  # t + s - R s, all in units of 1/den
                        shifted.add((tuple(map(tuple, o.rot)), tuple(int(v) % den for v in t)))
                    if shifted == listed:
                        return gemmi.SpaceGroup(g.xhm()), s / den
    return None


def _cell_from_orth(m):
    length = [np.linalg.norm(m[:, i]) for i in range(3)]

    def angle(u, v):
        return np.degrees(np.arccos(np.clip(np.dot(m[:, u], m[:, v]) / (length[u] * length[v]), -1, 1)))

    return gemmi.UnitCell(*length, angle(1, 2), angle(0, 2), angle(0, 1))


def _cart_to_ucif(ucart, cell):
    rc = cell.reciprocal()
    ninv = np.diag([1 / rc.a, 1 / rc.b, 1 / rc.c])
    f = np.array(cell.frac.mat.tolist())
    u = ninv @ f @ ucart @ f.T @ ninv
    return gemmi.SMat33d(u[0, 0], u[1, 1], u[2, 2], u[0, 1], u[0, 2], u[1, 2])


def _transformed(st, rot, tran, spacegroup):
    """st with fractional coordinates x' = rot x + tran (the cell changes accordingly) in the given space group."""
    old_orth = np.array(st.cell.orth.mat.tolist())
    new_orth = old_orth @ np.linalg.inv(rot)  # same Cartesian positions: old_orth @ x == new_orth @ (rot @ x)
    cell = _cell_from_orth(new_orth)
    # gemmi rebuilds the orthogonalisation in its own orientation; q rotates our frame into it
    q = np.array(cell.orth.mat.tolist()) @ np.linalg.inv(new_orth)
    out = gemmi.SmallStructure()
    out.name = st.name
    out.cell = cell
    out.spacegroup_hm = spacegroup.xhm()
    out.determine_and_set_spacegroup("2")
    for s in st.sites:
        n = s.clone()
        x = rot @ np.array([s.fract.x, s.fract.y, s.fract.z]) + tran
        n.fract = gemmi.Fractional(*(x % 1.0))
        if s.aniso.nonzero():
            n.aniso = _cart_to_ucif(q @ ucif_to_cart(s.aniso, st.cell) @ q.T, cell)
        out.add_site(n)
    return out


def to_reference(st):
    """The structure expressed in the reference setting of its space group, and whether it changed.

    st must have st.spacegroup set. The input is not modified.
    """
    op = st.spacegroup.basisop.inverse()  # coordinates in this setting -> reference setting
    if op == gemmi.Op("x,y,z"):
        return st, False
    rot = np.array(op.rot, dtype=float) / gemmi.Op.DEN
    tran = np.array(op.tran, dtype=float) / gemmi.Op.DEN
    return _transformed(st, rot, tran, reference_spacegroup(st.spacegroup.number)), True


def _candidate_matrices(system):
    """Integer changes of axes with determinant 1 that can keep a reference setting: for monoclinic (unique axis
    b) matrices mixing a and c with entries -2..2, b kept or reversed; for orthorhombic the axis permutations
    with signs."""
    if system == "monoclinic":
        for p, q, r, s in np.ndindex(5, 5, 5, 5):
            for sb in (1, -1):
                m = np.array([[p - 2, 0, q - 2], [0, sb, 0], [r - 2, 0, s - 2]])
                if round(np.linalg.det(m)) == 1:
                    yield m
    else:
        for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
            for signs in np.ndindex(2, 2, 2):
                m = np.zeros((3, 3), dtype=int)
                for row, col in enumerate(perm):
                    m[row, col] = 1 - 2 * signs[row]
                if round(np.linalg.det(m)) == 1:
                    yield m


@lru_cache(maxsize=None)
def _setting_choices(xhm):
    """Changes of coordinates x' = m x + t (t in eighths) under which the operations of a monoclinic or
    orthorhombic reference setting are again its own, i.e. other cells of the same setting."""
    den = gemmi.Op.DEN
    sg = gemmi.SpaceGroup(xhm)
    ops = list(sg.operations())
    rots = [np.array(o.rot) // den for o in ops]
    trans = [np.array(o.tran) / den for o in ops]
    key = {(tuple(map(tuple, r)), tuple(np.rint(t * 8).astype(int) % 8)) for r, t in zip(rots, trans)}
    shifts = [np.array(t) / 8 for t in np.ndindex(8, 8, 8)]
    out = []
    for m in _candidate_matrices(sg.crystal_system_str()):
        minv = np.rint(np.linalg.inv(m)).astype(int)
        conj = [(m @ r @ minv, m @ t) for r, t in zip(rots, trans)]
        for t0 in shifts:
            if {(tuple(map(tuple, r)), tuple(np.rint((t + t0 - r @ t0) * 8).astype(int) % 8)) for r, t in conj} == key:
                out.append((m, t0))
                break
    return out


def _cell_rank(system, orth):
    """The official COD databases' preference among cells of one setting (smaller is better); None if not
    allowed. Monoclinic: beta >= 90 degrees, then the shortest a + c, then the shortest a; orthorhombic: the
    shortest a, then b."""
    a, b, c = (np.linalg.norm(orth[:, i]) for i in range(3))
    if system == "monoclinic":
        beta = np.degrees(np.arccos(np.clip(orth[:, 0] @ orth[:, 2] / (a * c), -1, 1)))
        return None if beta < 90 - 1e-6 else (round(a + c, 6), round(a, 6), round(beta, 6))
    return (round(a, 6), round(b, 6), round(c, 6))


def reduce_cell(st):
    """A structure in its reference setting moved to the cell the official COD databases store for it: triclinic
    cells Niggli-reduced; monoclinic and orthorhombic cells chosen among those that keep the space-group symbol
    (_setting_choices, _cell_rank): e.g. c' = c + 2a for P 1 21/c 1 (a fixed change of setting alone gives
    beta = 140 deg for P 1 21/n 1, or an acute beta for P 1 21/a 1), c' = c + a for C 1 2/c 1, the axes of
    P 21 21 21 in increasing length, a cyclic permutation for P b c a. Returns the structure, the matrix m of
    x' = m x + t and whether it changed; the input is not modified."""
    total = np.eye(3)
    sg = st.spacegroup
    if sg is None or sg.basisop != gemmi.Op("x,y,z"):
        return st, total, False
    system = sg.crystal_system_str()
    if system == "triclinic":
        gv = gemmi.GruberVector(st.cell, None, True)
        # the tolerance that reproduces COD24's choice best (98% of triclinic entries; the rest are near-ties)
        gv.niggli_reduce(1e-3 * st.cell.volume ** (2 / 3))
        op = gv.change_of_basis.inverse()  # gemmi's op maps reduced coordinates to the given ones
        if op == gemmi.Op("x,y,z"):
            return st, total, False
        m = np.array(op.rot, dtype=float) / gemmi.Op.DEN
        return _transformed(st, m, np.zeros(3), sg), m, True
    if system not in ("monoclinic", "orthorhombic"):
        return st, total, False
    for _ in range(10):  # a monoclinic pass mixes a and c by at most twice the other
        orth = np.array(st.cell.orth.mat.tolist())
        best = None
        for m, t0 in _setting_choices(sg.xhm()):
            rank = _cell_rank(system, orth @ np.linalg.inv(m))
            if rank is not None and (best is None or rank < best[0]):
                best = (rank, m, t0)
        m, t0 = best[1], best[2]
        if np.array_equal(m, np.eye(3)) and not t0.any():
            break
        st = _transformed(st, m.astype(float), t0, sg)
        total = m @ total
    return st, total, not np.array_equal(total, np.eye(3))


def constrain_cell(st):
    """Cell parameters made exact for the crystal system of st's space group (in its reference setting), as the
    official COD databases store them: angles of 90 and 120 degrees, a = b (tetragonal, trigonal, hexagonal)
    and a = b = c (cubic), taking b and c respectively. Some CIFs give cells that do not fit their space group
    (e.g. 89.99 or 89.0 degrees in an orthorhombic cell); unconstrained, symmetry-equivalent reflections would
    fall at different d. Coordinates stay fractional. Returns whether the cell changed.
    """
    c = st.cell
    a, b, cc, al, be, ga = c.a, c.b, c.c, c.alpha, c.beta, c.gamma
    system = st.spacegroup.crystal_system_str()
    if system == "monoclinic":
        al = ga = 90.0
    elif system == "orthorhombic":
        al = be = ga = 90.0
    elif system == "tetragonal":
        a, al, be, ga = b, 90.0, 90.0, 90.0
    elif system in ("trigonal", "hexagonal"):
        a, al, be, ga = b, 90.0, 90.0, 120.0
    elif system == "cubic":
        a, b, al, be, ga = cc, cc, 90.0, 90.0, 90.0
    new = (a, b, cc, al, be, ga)
    if new == (c.a, c.b, c.c, c.alpha, c.beta, c.gamma):
        return False
    st.cell = gemmi.UnitCell(*new)
    return True


def normalise_symbol(hm):
    """Hermann-Mauguin symbol as gemmi reads it, for three notations seen in COD: a cell-choice qualifier on the
    short symbol ('C2:b1'), an origin written out ('F d d d {origin @ -1 @ d d d}', origin choice 2) and a
    superspace-group symbol ('Cmca(00\\g)s00', COD 2105669), whose basic space group precedes the wave vector."""
    hm = re.sub(r"\s*\{\s*origin\s*@\s*-1\b[^}]*\}\s*$", ":2", hm)
    hm = re.sub(r"^\s*([^()]+?)\s*\([^()]*\)\s*[0stqh]+\s*$", r"\1", hm)
    if ":" in hm:
        base, qualifier = (t.strip() for t in hm.split(":", 1))
        g = gemmi.find_spacegroup_by_name(base)
        if g is not None and g.qualifier == qualifier:
            return g.xhm()
    return hm


def listed_operations(st, block=None):
    """Symmetry operations of a structure as read from its CIF: the listed ones, else those of its Hall symbol
    (which may carry a change of basis; gemmi's st.spacegroup_hall, which a correction may have cleared). None if
    there are none or they cannot be read. block is not used (kept for callers)."""
    hall = st.spacegroup_hall.strip()
    try:
        if st.symops:
            return [gemmi.Op(o) for o in st.symops]
        if hall not in ("", "?", "."):
            return list(gemmi.symops_from_hall(hall))
    except RuntimeError:
        return None
    return None


def _is_group(ops):
    den = gemmi.Op.DEN
    key = {(tuple(map(tuple, o.rot)), tuple(t % den for t in o.tran)) for o in ops}
    return len(key) == len(ops) and all(
        (tuple(map(tuple, p.rot)), tuple(t % den for t in p.tran)) in key for a in ops for b in ops for p in [a * b])


def from_operations(st, ops):
    """A structure given with explicit operations in a setting gemmi does not know (unusual axes or centring,
    a shifted origin, a supercell), in the reference setting of its space group, with the matrix P of the change
    (x' = P x + p); None if it cannot be trusted.

    The sites are expanded with the operations and spglib identifies the space group of the result. Accepted only
    if the operations form a group that preserves the cell metric and spglib finds exactly that many operations
    in the cell, each of them one of the listed ones: the standardised structure is then the listed one.
    """
    if not ops or not _is_group(ops):
        return None
    orth = np.array(st.cell.orth.mat.tolist())
    metric = orth.T @ orth
    rots = [np.array(o.rot, dtype=float) / gemmi.Op.DEN for o in ops]
    if not all(np.allclose(r.T @ metric @ r, metric, rtol=0, atol=1e-3 * np.abs(metric).max()) for r in rots):
        return None
    images = unit_cell_images(st, ops)
    listed = {(tuple(map(tuple, np.rint(r).astype(int))), tuple(np.rint(np.array(o.tran) / gemmi.Op.DEN * 24) % 24))
              for r, o in zip(rots, ops)}
    for symprec in (1e-3, 1e-2):  # Å; the looser pass covers coordinates like 0.333 for 1/3
        try:
            ds = spglib.get_symmetry_dataset((orth.T, images.xyz, images.site), symprec=symprec)
        except spglib.error.SpglibError:
            continue
        if ds is None or len(ds.rotations) != len(ops):
            continue
        found = {(tuple(map(tuple, r)), tuple(np.rint(t * 24) % 24)) for r, t in zip(ds.rotations, ds.translations)}
        if found != listed:
            continue
        hall = spglib.get_spacegroup_type(ds.hall_number).hall_symbol
        standard = gemmi.find_spacegroup_by_ops(gemmi.symops_from_hall(hall))
        if standard is None:
            return None
        rot = np.array(ds.transformation_matrix, dtype=float)
        return _transformed(st, rot, np.array(ds.origin_shift), standard), rot
    return None
