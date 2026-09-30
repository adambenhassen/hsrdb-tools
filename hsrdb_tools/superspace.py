"""Diffraction patterns of modulated structures from their superspace description (msCIF, as JANA writes it).

In a modulated structure the atoms are displaced, or their sites partly occupied, by waves whose period does not
fit the basic unit cell (wave vectors q). A CIF describes it in (3+d)-dimensional superspace: the basic structure,
the modulation functions of each atom along its internal coordinate v = t + q.xbar (JANA's x4), and superspace
symmetry operations. The diffraction pattern holds main reflections h and satellites h + m.q.

Each atom is a string in superspace, x_s(v) = (xbar + u(v), v + Q u(v)). An operation (W, w) maps H.x_s to
(H W).x_s + H.w, so the structure factor of H = (h, m) is

    F(H) = sum_atoms f / stab  sum_ops exp(2 pi i H.w) < p(v) T_K(v) exp(2 pi i [k.xbar + S_K.u(v) + n.v]) >_v

with K = (k, n) = H W, S_K = k + n Q the scattering vector of K, stab the number of operations that leave the
atom's basic position in place, and <>_v the mean over v: over one period, or over the crenel interval the atom
occupies (incommensurate), or over the v values that the section t0 selects in the supercell (commensurate).
JANA writes the average occupancy of a crenel atom (within its interval of width w the site holds occ / w), and the
occupancy within the interval for an atom that only a sawtooth function confines.
The structure factors agree with those JANA lists in the CIFs themselves (dev/verify_modulated.py).
"""

import math
import re
from dataclasses import dataclass, field
from fractions import Fraction

import gemmi
import numpy as np

from . import pattern

# coefficients of modulation functions whose definitions are specific to JANA's implementation and are not
# reproduced here (JANA also writes the definitions of functions it does not use: only coefficients count)
_UNSUPPORTED = re.compile(r"(crenel_ortho|xharm|_displace_ortho|_u_ortho|_occ_ortho|zigzag|_b_fourier|"
                          r"_adp_[c-f]_fourier)\w*_param_(coeff|cos|sin)")
_FUNCTIONS = re.compile(r"(_fourier|_legendre)\w*_param_(cos|sin|coeff)|_special_func_(crenel|sawtooth)")
_IJ = {"u11": (0, 0), "u22": (1, 1), "u33": (2, 2), "u12": (0, 1), "u13": (0, 2), "u23": (1, 2)}
# quadrature: points per period (trapezoid, exact for trigonometric polynomials of lower degree) or Gauss-Legendre
# points per crenel interval, chosen for the highest frequency of each integrand (see _grid), within these bounds
MIN_POINTS, MAX_POINTS = {1: 16, 2: 16}, {1: 1024, 2: 160}
SATELLITE_LIMIT = 0.5
MAX_ORDER = {1: 64, 2: 12}
# largest R against the structure factors a CIF lists (agreement) for which it is built: the reproduced CIFs agree
# within 0.9% (rounding of the listed values), the others differ by 2% to 50%
AGREEMENT_LIMIT = 0.01


class ModulationError(ValueError):
    """The CIF's modulation cannot be evaluated (a function this module does not reproduce, missing data)."""


@dataclass
class Atom:
    site: int  # index into the structure the builder reads (st.sites)
    element: str
    xyz: np.ndarray  # basic position, in the CIF's setting
    occ: float  # average occupancy, as the CIF gives it
    uq: np.ndarray  # T = exp(-2 pi^2 S uq S) for fractional S, from the builder's displacement parameters
    stab: int = 1
    disp: list = field(default_factory=list)  # (axis, wave (d,), cos, sin), fractional
    occf: list = field(default_factory=list)  # (wave, cos, sin), absolute
    uf: list = field(default_factory=list)  # ((i, j), wave, cos, sin), U_ij against reciprocal-axis lengths
    crenel: tuple | None = None  # (centre, width) in v: the atom is present only there
    occ_inside: bool = False  # occ is the occupancy within the interval (a sawtooth's), not the average
    sawtooth: np.ndarray | None = None  # amplitude (3,), fractional, across the crenel interval
    legendre_disp: list = field(default_factory=list)  # (axis, order, coefficient), within the crenel interval
    legendre_u: list = field(default_factory=list)
    legendre_occ: list = field(default_factory=list)


@dataclass
class Modulated:
    cell: gemmi.UnitCell  # basic cell, CIF setting
    q: np.ndarray  # (d, 3) wave vectors
    ops: list  # (W (3+d, 3+d), w (3+d,)) superspace operations, centring included
    atoms: list
    ssg: str  # superspace group symbol as the CIF gives it (may be empty)
    t0: np.ndarray | None = None  # commensurate: the section
    supercell: np.ndarray | None = None  # commensurate: basic lattice translations inside the supercell
    f2: dict = field(default_factory=dict, repr=False)  # |F|^2 (Cu Ka1) by symmetry-unique reflection, see _intensities

    @property
    def d(self):
        return len(self.q)


def _num(v):
    if v is None:
        return None
    s = gemmi.cif.as_string(v).strip()
    if s in ("", "?", "."):
        return None
    try:
        return float(re.sub(r"\([^)]*\)", "", s))
    except ValueError:
        raise ModulationError(f"not a number: {s!r}") from None


def _op(text, n):
    """Superspace operation 'x1+1/2,-x2,x3,x4-x1' as (W, w)."""
    names = [f"x{i + 1}" for i in range(n)]
    rows = text.replace(" ", "").lower().split(",")
    if len(rows) != n:
        raise ModulationError(f"superspace operation {text!r} is not {n}-dimensional")
    W, w = np.zeros((n, n)), np.zeros(n)
    for r, expr in enumerate(rows):
        for sign, term in re.findall(r"([+-]?)([^+-]+)", expr):
            s = -1.0 if sign == "-" else 1.0
            m = re.fullmatch(r"(\d*\.?\d*/?\d*)\*?(x\d)", term)
            try:
                if m:
                    W[r, names.index(m.group(2))] += s * (float(Fraction(m.group(1))) if m.group(1) else 1.0)
                else:
                    w[r] += s * float(Fraction(term))
            except (ValueError, ZeroDivisionError):
                raise ModulationError(f"unreadable superspace operation {text!r}") from None
    return W, w


def _loop(block, prefix, cols):
    """Rows of a loop as lists (missing optional columns '?')."""
    t = block.find(prefix, cols)
    return [[t[i][j] if t[i].has(j) else "?" for j in range(len(cols))] for i in range(len(t))]


def _terms(block, prefix, keys, params):
    """Rows [atom_site_label, *keys, *params] of a modulation-function loop, from either msCIF layout: one loop with
    everything, or a loop of ids (with the label and keys) and a loop of their parameters (prefix + 'param_id').
    Rows without a label ('?', '.' or '', JANA's empty loops) are left out."""
    rows = _loop(block, prefix, ["atom_site_label"] + keys + [f"param_{p}" for p in params])
    if not rows:
        by_id = {r[0]: r[1:] for r in _loop(block, prefix, ["id", "atom_site_label"] + keys)}
        rows = []
        for r in _loop(block, prefix + "param_", ["id"] + params):
            if r[0] not in by_id:
                raise ModulationError(f"{prefix}param_id {r[0]} has no {prefix}id")
            rows.append(by_id[r[0]] + r[1:])
    return [r for r in rows if gemmi.cif.as_string(r[0]).strip() not in ("?", ".", "")]


def _used(value):
    """A value other than ?, . or a number that is zero."""
    if value in ("?", "."):
        return False
    try:
        return float(re.sub(r"\([^)]*\)", "", value)) != 0
    except ValueError:
        return True


def _used_tags(block):
    """Lower-case tags with at least one value other than ?, . or zero (JANA writes empty loops for unused
    functions, and zero coefficients)."""
    out = set()
    for item in block:
        if item.pair:
            if _used(item.pair[1]):
                out.add(item.pair[0].lower())
        elif item.loop:
            tags, vals = [t.lower() for t in item.loop.tags], item.loop.values
            for j, t in enumerate(tags):
                if any(_used(v) for v in vals[j::len(tags)]):
                    out.add(t)
    return out


def _lattice_points(M):
    """Basic lattice translations inside the supercell whose axes are the rows of M (in basic units)."""
    Minv = np.linalg.inv(M)
    lo = np.floor(np.minimum(M, 0).sum(axis=0)).astype(int) - 1
    hi = np.ceil(np.maximum(M, 0).sum(axis=0)).astype(int) + 1
    L = np.stack(np.meshgrid(*[np.arange(a, b + 1) for a, b in zip(lo, hi)], indexing="ij"), -1).reshape(-1, 3)
    f = L @ Minv
    pts = L[np.all((f > -1e-9) & (f < 1 - 1e-9), axis=1)]
    if len(pts) != round(abs(np.linalg.det(M))):
        raise ModulationError("inconsistent commensurate supercell")
    return pts


def _waves(block, d):
    """Wave vector of each Fourier term (seq id -> integer coefficients of the q's)."""
    for pre in ("_jana_atom_site_fourier_wave_vector_q", "_atom_site_fourier_wave_vector_q"):
        lp = block.find_loop(pre + "1_coeff").get_loop()
        if lp is None:
            continue
        tags = [t.lower() for t in lp.tags]
        try:
            seq = next(j for j, t in enumerate(tags) if t.endswith("wave_vector_seq_id"))
            cols = [tags.index(f"{pre}{i + 1}_coeff") for i in range(d)]
        except (StopIteration, ValueError):
            raise ModulationError("incomplete Fourier wave vector list") from None
        return {lp[r, seq]: np.array([int(_num(lp[r, c])) for c in cols]) for r in range(lp.length())}
    return None


def _explicit_waves(block, q):
    """Wave vectors listed as components (_atom_site_Fourier_wave_vector_x/y/z), as integer multiples of the q's;
    None if not listed. JANA's older CIFs leave out a component column even where it is not zero, so only the
    listed components are compared: wave i <= d is q_i if they agree (JANA numbers the q's first), any other wave
    must match exactly one combination of the q's with coefficients up to 4."""
    for axis in "xyz":
        lp = block.find_loop(f"_atom_site_fourier_wave_vector_{axis}").get_loop()
        if lp is not None:
            break
    else:
        return None
    tags = [t.lower() for t in lp.tags]
    seq = next((j for j, t in enumerate(tags) if t.endswith("wave_vector_seq_id")), None)
    listed = [i for i, a in enumerate("xyz") if f"_atom_site_fourier_wave_vector_{a}" in tags]
    if seq is None or not listed:
        raise ModulationError("Fourier wave vectors without ids or components")
    d = len(q)
    grid = np.stack(np.meshgrid(*[np.arange(-4, 5)] * d, indexing="ij"), -1).reshape(-1, d)
    out = {}
    for r in range(lp.length()):
        vec = np.array([_num(lp[r, tags.index(f"_atom_site_fourier_wave_vector_{'xyz'[i]}")]) for i in listed])
        key = lp[r, seq]
        i = int(key) - 1 if key.isdigit() else -1
        if 0 <= i < d and np.allclose(q[i, listed], vec, atol=1e-4):
            out[key] = np.eye(d, dtype=int)[i]
            continue
        match = grid[np.all(np.abs(grid @ q[:, listed] - vec) < 1e-4, axis=1)]
        if len(match) != 1:
            raise ModulationError(f"Fourier wave vector {key} ({vec} in the listed components) does not identify "
                                  "one combination of the modulation wave vectors")
        out[key] = match[0]
    return out


def _ssg_operations(block):
    """The listed superspace operations as (W, w), or [] (their dimension from the operations themselves)."""
    texts = [gemmi.cif.as_string(r[0]) for r in block.find("_space_group_symop_ssg_", ["operation_algebraic"])]
    return [_op(t, len(t.split(","))) for t in texts]


def operations_3d(block):
    """The distinct 3D parts of the listed superspace operations, as gemmi triplets: the space group of the basic
    structure, for CIFs that list no other symmetry. [] if there are none or they cannot be read."""
    try:
        ops = _ssg_operations(block)
    except ModulationError:
        return []
    out = {}
    for W, w in ops:
        if W.shape[0] < 3 or not np.allclose(W[:3, 3:], 0):
            return []
        op = gemmi.Op("x,y,z")
        op.rot = np.rint(W[:3, :3] * gemmi.Op.DEN).astype(int).tolist()
        op.tran = (np.rint(w[:3] * gemmi.Op.DEN).astype(int) % gemmi.Op.DEN).tolist()
        out.setdefault(op.triplet(), None)
    return list(out)


def states_modulation(block):
    """Whether the CIF describes its structure as modulated (a modulation dimension or wave vectors)."""
    dim = _num(block.find_value("_cell_modulation_dimension"))
    q = block.find("_cell_wave_vector_", ["x"])
    return bool(dim and dim > 0) or any(_num(r[0]) for r in q)


def composite(block):
    """Whether the CIF describes a composite (intergrowth) crystal: subsystems with different basic lattices."""
    n = _num(block.find_value("_cell_subsystems_number"))
    return n is not None and n > 1


def read(block, st, rows, shift=(0.0, 0.0, 0.0)):
    """The CIF's superspace description, or None when it has no modulation functions (an ordinary structure, or
    only the basic structure of a modulated one).

    st: the builder's structure in the CIF's setting (elements resolved, displacement parameters completed),
    rows: its atom_site row for each site, shift: the origin shift st's coordinates carry. Raises
    ModulationError when the modulation cannot be evaluated."""
    if composite(block):
        raise ModulationError("composite structure: its subsystems have different basic lattices, which one 3D "
                              "cell cannot hold")
    tags = _used_tags(block)
    if not any(_FUNCTIONS.search(t) for t in tags):
        return None
    unsupported = sorted({_UNSUPPORTED.search(t).group(1).strip("_") for t in tags if _UNSUPPORTED.search(t)})
    if unsupported:
        raise ModulationError(f"modulation functions not supported: {', '.join(unsupported)}")
    q = np.array([[_num(v) for v in r] for r in block.find("_cell_wave_vector_", ["x", "y", "z"])], dtype=float)
    q = q.reshape(-1, 3)
    if not len(q) or not np.isfinite(q).all():
        raise ModulationError("modulation functions without wave vectors")
    d = len(q)
    if d > 2:
        raise ModulationError(f"({3}+{d})-dimensional modulation not supported")
    ops = _ssg_operations(block)
    if not ops:
        raise ModulationError("modulation functions without superspace symmetry operations")
    if any(W.shape[0] != 3 + d for W, _ in ops):
        raise ModulationError(f"superspace operations do not have {3 + d} dimensions")
    t0 = [_num(block.find_value(f"_jana_cell_commen_t_section_{i + 1}")) for i in range(d)]
    t0 = np.array(t0) if all(x is not None for x in t0) else None
    supercell = None
    if t0 is not None:
        M = np.array([[_num(block.find_value(f"_jana_cell_commen_supercell_matrix_{i}_{j}")) for j in (1, 2, 3)]
                      for i in (1, 2, 3)], dtype=float)
        if not np.isfinite(M).all():
            raise ModulationError("commensurate section without supercell")
        # a commensurate q is a fraction that CIFs write rounded (13/47 as 0.276596, COD 2310038); rounded, the
        # satellite of highest order misses the lattice vector it falls on and becomes a line at d ~ 10^5 A
        exact = np.array([[float(Fraction(float(x)).limit_denominator(1000)) for x in row] for row in q])
        if np.abs(exact - q).max() < 1e-5:
            q = exact
        supercell = _lattice_points(M)
    waves = _waves(block, d)
    if waves is None:
        waves = _explicit_waves(block, q)

    def wave(seq):
        if waves is not None:
            if seq not in waves:
                raise ModulationError(f"unknown wave vector {seq!r}")
            return waves[seq]
        if d != 1 or not seq.isdigit():
            raise ModulationError("Fourier terms without wave vectors")
        return np.array([int(seq)])

    labels = [gemmi.cif.as_string(v) for v in block.find_values("_atom_site_label")]
    cell = st.cell
    rc = cell.reciprocal()
    astar = np.array([rc.a, rc.b, rc.c])
    gstar = np.array(cell.reciprocal_metric_tensor().as_mat33().tolist())
    atoms, by_label = [], {}
    for i, (s, row) in enumerate(zip(st.sites, rows)):
        if row is None:
            raise ModulationError("atoms added by a correction have no modulation")
        if s.aniso.nonzero():
            u = s.aniso
            uq = np.array([[u.u11, u.u12, u.u13], [u.u12, u.u22, u.u23], [u.u13, u.u23, u.u33]]) * np.outer(astar, astar)
        else:
            uq = s.u_iso * gstar
        a = Atom(i, s.element.name, np.array([s.fract.x, s.fract.y, s.fract.z]) + shift, s.occ, uq)
        atoms.append(a)
        by_label.setdefault(labels[row], []).append(a)

    all_labels = set(labels)

    def each(label):
        """The structure's sites of an atom_site label (none if the builder left the row out, e.g. no coordinates);
        a label the CIF does not list is an error: its modulation would be lost silently."""
        label = gemmi.cif.as_string(label).strip()
        if label in ("?", ".", ""):  # a row of JANA's empty loops
            return []
        if label not in all_labels:
            raise ModulationError(f"modulation function for {label!r}, which is not an atom_site label")
        return by_label.get(label, [])

    def value(v):
        x = _num(v)
        return 0.0 if x is None else x

    for r in _terms(block, "_atom_site_displace_Fourier_", ["axis", "wave_vector_seq_id"], ["cos", "sin"]):
        for a in each(r[0]):
            a.disp.append(("xyz".index(gemmi.cif.as_string(r[1]).lower()), wave(gemmi.cif.as_string(r[2])),
                           value(r[3]), value(r[4])))
    for r in _terms(block, "_atom_site_occ_Fourier_", ["wave_vector_seq_id"], ["cos", "sin"]):
        for a in each(r[0]):
            a.occf.append((wave(gemmi.cif.as_string(r[1])), value(r[2]), value(r[3])))
    for r in _terms(block, "_atom_site_U_Fourier_", ["tens_elem", "wave_vector_seq_id"], ["cos", "sin"]):
        for a in each(r[0]):
            a.uf.append((_IJ[gemmi.cif.as_string(r[1]).lower()], wave(gemmi.cif.as_string(r[2])), value(r[3]),
                         value(r[4])))
    for r in _loop(block, "_atom_site_occ_special_func_", ["atom_site_label", "crenel_c", "crenel_w"]):
        if _num(r[1]) is not None:
            for a in each(r[0]):
                a.crenel = (_num(r[1]), _num(r[2]))
    for r in _loop(block, "_atom_site_displace_special_func_", ["atom_site_label", "sawtooth_ax", "sawtooth_ay",
                                                                  "sawtooth_az", "sawtooth_c", "sawtooth_w"]):
        if _num(r[4]) is not None:
            for a in each(r[0]):
                a.sawtooth = np.array([value(x) for x in r[1:4]])
                interval = (_num(r[4]), _num(r[5]))
                if a.crenel is None:
                    # the atom occupies the sawtooth's interval; JANA then writes the occupancy within it (COD
                    # 4002590: Na and Mn alternate on one site, each at occupancy 1 over half of x4)
                    a.crenel, a.occ_inside = interval, True
                elif not np.allclose(a.crenel, interval):
                    raise ModulationError("sawtooth and crenel intervals of an atom differ")
    for pre in ("_jana_atom_site_displace_legendre_", "_atom_site_displace_legendre_"):
        for r in _loop(block, pre, ["atom_site_label", "axis", "param_order", "param_coeff"]):
            for a in each(r[0]):
                a.legendre_disp.append(("xyz".index(gemmi.cif.as_string(r[1]).lower()), int(value(r[2])),
                                        value(r[3])))
    for r in _loop(block, "_jana_atom_site_u_legendre_", ["atom_site_label", "tens_elem", "param_order",
                                                            "param_coeff"]):
        for a in each(r[0]):
            a.legendre_u.append((_IJ[gemmi.cif.as_string(r[1]).lower()], int(value(r[2])), value(r[3])))
    for r in _loop(block, "_jana_atom_site_occ_legendre_", ["atom_site_label", "param_order", "param_coeff"]):
        for a in each(r[0]):
            a.legendre_occ.append((int(value(r[1])), value(r[2])))
    # a function the CIF gives values for must have reached the atoms (another loop layout would otherwise be
    # ignored silently)
    read = {"displace": any(a.disp for a in atoms), "occ": any(a.occf for a in atoms), "u": any(a.uf for a in atoms),
            "legendre": any(a.legendre_disp or a.legendre_u or a.legendre_occ for a in atoms),
            "special": any(a.crenel is not None for a in atoms)}
    for t in tags:
        m = re.search(r"_(displace|occ|u)_fourier_param_(cos|sin)|_(legendre)_param_coeff|_(special)_func_", t)
        if m and not read[next(g for g in m.groups() if g in read)]:
            raise ModulationError(f"{t} has values that could not be assigned to atoms")
    for a in atoms:
        if (a.legendre_disp or a.legendre_u or a.legendre_occ or a.sawtooth is not None) and a.crenel is None:
            raise ModulationError("Legendre or sawtooth functions without crenel interval")
        if a.crenel is not None and not 0 < a.crenel[1] <= 1:
            raise ModulationError(f"crenel width {a.crenel[1]}")
        if a.crenel is not None and d != 1:
            raise ModulationError("crenel functions in more than one internal dimension")
    ssg = gemmi.cif.as_string(block.find_value("_space_group_ssg_name") or "").strip()
    return Modulated(cell, q, ops, atoms, "" if ssg in ("?", ".") else ssg, t0, supercell)


def check_operations(model, ops3d):
    """ModulationError unless the superspace operations restricted to 3D are the structure's space group."""
    den = gemmi.Op.DEN
    key = lambda r, t: (tuple(np.rint(r).astype(int).ravel()), tuple(np.rint(np.asarray(t) * 24).astype(int) % 24))  # noqa: E731
    ss = {key(W[:3, :3], w[:3]) for W, w in model.ops}
    if ss != {key(np.array(o.rot) / den, np.array(o.tran) / den) for o in ops3d}:
        raise ModulationError("the superspace operations do not match the space group")


def _static(a):
    """Whether atom a has no modulation (its scattering is the basic structure's; no satellites)."""
    return not (a.disp or a.occf or a.uf or a.crenel is not None)


def _bandwidth(model, a, smax):
    """Per internal axis, an upper estimate of the highest frequency (cycles per period of v) in the integrand of
    atom a for scattering vectors up to smax (1/A), besides n.v. A displacement wave k of amplitude A (A) makes
    exp(2 pi i S.u) a phase modulation of index beta = 2 pi smax A, whose Bessel components J_n(beta) fall below
    1e-9 of the total beyond n = beta + 10 + 3 beta^(1/3) (along k); several waves add their beta |k| and need the
    margin once. Occupancy waves add their order, and displacement-parameter waves act through exp(-2 pi^2 S U S)
    like a displacement wave of index 2 pi^2 smax^2 U."""
    lengths = np.array([model.cell.a, model.cell.b, model.cell.c])
    rc = model.cell.reciprocal()
    astar = np.array([rc.a, rc.b, rc.c])
    b = np.zeros(model.d)
    disp, uwaves = {}, {}
    for axis, k, c, s in a.disp:
        disp.setdefault(tuple(k), np.zeros(3))[axis] += math.hypot(c, s) * lengths[axis]
    for (i, j), k, c, s in a.uf:  # U_ij against reciprocal lengths -> A^2 (bound by the smallest a*)
        uwaves[tuple(k)] = uwaves.get(tuple(k), 0.0) + math.hypot(c, s) * astar[i] * astar[j] / astar.min() ** 2
    # a sum of phase-modulation waves: the frequencies of each (index beta, along k) add; the Bessel tail of the
    # product needs its margin once, along the highest wave order
    total, kmax = 0.0, np.zeros(model.d)
    for waves, index in ((disp, lambda x: 2 * math.pi * smax * float(np.abs(x).sum())),
                         (uwaves, lambda x: 2 * math.pi**2 * smax**2 * x)):
        for k, amp in waves.items():
            beta = index(amp)
            b += np.abs(k) * beta
            total += beta
            kmax = np.maximum(kmax, np.abs(k))
    b += kmax * (10 + 3 * total ** (1 / 3))
    for k, _, _ in a.occf:  # occupancy waves multiply: their orders add
        b += np.abs(k)
    return b


def _samples(model, a, W, w, nmax=0, smax=0.0):
    """Internal coordinates v (n, d) of atom a and their weights, for its image under (W, w), for reflections with
    satellite indices up to nmax and scattering vectors up to smax (incommensurate: the quadrature follows them)."""
    if model.t0 is None:
        top = nmax + _bandwidth(model, a, smax)  # highest frequency of the integrand, per internal axis
        lo, hi = MIN_POINTS[model.d], MAX_POINTS[model.d]
        if a.crenel is not None:
            c, width = a.crenel
            n = min(hi, max(lo, math.ceil(4 * top[0] * width) + 16))
            x, gw = np.polynomial.legendre.leggauss(n)
            return (c + width * x / 2)[:, None], gw * width / 2
        # the midpoint rule over a period with n points is exact for frequencies below n
        n = np.ceil(top).astype(int) + 4
        if n.max() > hi:
            raise ModulationError(f"the modulation needs {n.max()} quadrature points per period (more than {hi})")
        axes = [(np.arange(max(lo, k)) + 0.5) / max(lo, k) for k in n]
        v = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, model.d)
        return v, np.full(len(v), 1 / len(v))
    # commensurate: the image's atoms sit at v' = t0 + Q.(xbar' + L) for L in the supercell, and v' = R_I v + R_M xbar + w4
    xp = W[:3, :3] @ a.xyz + w[:3]
    vp = model.t0 + (xp + model.supercell) @ model.q.T
    v = (vp - (W[3:, :3] @ a.xyz + w[3:])) @ np.linalg.inv(W[3:, 3:]).T
    weight = np.full(len(v), 1 / len(v))
    if a.crenel is not None:
        c, width = a.crenel
        dv = np.abs((v[:, 0] - c + 0.5) % 1.0 - 0.5)
        if np.any(np.abs(dv - width / 2) < 1e-6):
            raise ModulationError("a commensurate section falls on a crenel edge")
        weight = weight * (dv < width / 2)
    return v, weight


def _functions(model, a, v):
    """Displacement u (n, 3), occupancy p (n,) and displacement quadratic form (n, 3, 3) of atom a at v."""
    n = len(v)
    u = np.zeros((n, 3))
    p = np.full(n, a.occ / a.crenel[1] if a.crenel is not None and not a.occ_inside else a.occ)
    U = np.broadcast_to(a.uq, (n, 3, 3)).copy()
    rc = model.cell.reciprocal()
    astar = np.array([rc.a, rc.b, rc.c])
    for axis, k, c, s in a.disp:
        arg = 2 * math.pi * (v @ k)
        u[:, axis] += c * np.cos(arg) + s * np.sin(arg)
    for k, c, s in a.occf:
        arg = 2 * math.pi * (v @ k)
        p += c * np.cos(arg) + s * np.sin(arg)
    for (i, j), k, c, s in a.uf:
        val = (c * np.cos(2 * math.pi * (v @ k)) + s * np.sin(2 * math.pi * (v @ k))) * astar[i] * astar[j]
        U[:, i, j] += val
        if i != j:
            U[:, j, i] += val
    if a.crenel is not None:
        c, width = a.crenel
        x = 2 * ((v[:, 0] - c + 0.5) % 1.0 - 0.5) / width  # -1 .. 1 across the interval
        if a.sawtooth is not None:
            u += np.outer(x, a.sawtooth)
        for axis, order, coef in a.legendre_disp:
            u[:, axis] += coef * np.polynomial.legendre.Legendre.basis(order)(x)
        for order, coef in a.legendre_occ:
            p += coef * np.polynomial.legendre.Legendre.basis(order)(x)
        for (i, j), order, coef in a.legendre_u:
            val = coef * np.polynomial.legendre.Legendre.basis(order)(x) * astar[i] * astar[j]
            U[:, i, j] += val
            if i != j:
                U[:, j, i] += val
    return u, p, U


def structure_factors(model, H, wavelength=pattern.CU_KA1, scattering=None, dispersion=None, chunk=4096):
    """Complex F for superspace indices H (n, 3+d), per basic cell. scattering: element -> scattering length
    (neutrons) or function of sin^2(theta)/lambda^2 (electrons) instead of X-ray form factors; dispersion:
    element -> f' + i f'' instead of Cromer-Liberman (both for comparisons with a refinement's own values)."""
    H = np.asarray(H, dtype=float)
    out = np.zeros(len(H), complex)
    if not len(H):
        return out
    gstar = np.array(model.cell.reciprocal_metric_tensor().as_mat33().tolist())
    elements = sorted({a.element for a in model.atoms})
    # the images' indices K = H W reach satellite orders up to about those of H; the largest |K[3:]| bounds them
    nmax = max(int(np.abs(H @ W)[:, 3:].max()) for W, _ in model.ops)
    S = H[:, :3] + H[:, 3:] @ model.q
    smax = float(np.sqrt(np.einsum("ni,ij,nj->n", S, gstar, S).max()))
    samples = [[(*_samples(model, a, W, w, nmax, smax), W, w) for W, w in model.ops] for a in model.atoms]
    funcs = [[_functions(model, a, v) for v, _, _, _ in s] for a, s in zip(model.atoms, samples)]
    for lo in range(0, len(H), chunk):
        h = H[lo:lo + chunk]
        S = h[:, :3] + h[:, 3:] @ model.q
        stol2 = np.einsum("ni,ij,nj->n", S, gstar, S) / 4
        if scattering:
            f = {el: (scattering[el](stol2) if callable(scattering[el]) else np.full(len(h), scattering[el]))
                 .astype(complex) for el in elements}
        else:
            f = {el: pattern._f0(el, stol2) + (dispersion[el] if dispersion and el in dispersion else
                                                complex(*pattern._anomalous(el, wavelength))) for el in elements}
        for a, sa, fa in zip(model.atoms, samples, funcs):
            acc = np.zeros(len(h), complex)
            for (v, weight, W, w), (u, p, U) in zip(sa, fa):
                K = h @ W
                SK = K[:, :3] + K[:, 3:] @ model.q
                if _static(a) and model.t0 is None:  # the mean over v vanishes unless n = 0
                    T = np.exp(-2 * math.pi**2 * np.einsum("ni,ij,nj->n", SK, a.uq, SK))
                    main = ~K[:, 3:].any(axis=1)
                    acc += np.where(main, a.occ * T * np.exp(2j * math.pi * (K[:, :3] @ a.xyz + h @ w)), 0)
                    continue
                phase = (K[:, :3] @ a.xyz)[:, None] + SK @ u.T + K[:, 3:] @ v.T
                if a.uf or a.legendre_u:
                    T = np.exp(-2 * math.pi**2 * np.einsum("ni,gij,nj->ng", SK, U, SK))
                    acc += np.exp(2j * math.pi * (h @ w)) * ((np.exp(2j * math.pi * phase) * T) @ (p * weight))
                else:  # displacement parameters constant along v
                    T = np.exp(-2 * math.pi**2 * np.einsum("ni,ij,nj->n", SK, a.uq, SK))
                    acc += np.exp(2j * math.pi * (h @ w)) * T * (np.exp(2j * math.pi * phase) @ (p * weight))
            out[lo:lo + chunk] += f[a.element] * acc / a.stab
    return out


def _order_indices(d, k):
    """All m (n, d) with max |m_i| == k."""
    if k == 0:
        return np.zeros((1, d), dtype=int)
    g = np.stack(np.meshgrid(*[np.arange(-k, k + 1)] * d, indexing="ij"), -1).reshape(-1, d)
    return g[np.abs(g).max(axis=1) == k]


def _reflections(model, m, dmin, seen=None):
    """(h, k, l, m) for the given satellite indices m (n, d) with d >= dmin in the CIF's cell; for a commensurate
    structure only those whose scattering vector is not in seen (updated)."""
    rc = model.cell.reciprocal()
    astar = np.array([rc.a, rc.b, rc.c])
    gstar = np.array(model.cell.reciprocal_metric_tensor().as_mat33().tolist())
    smax = 1 / dmin
    shift = m @ model.q
    out = []
    for mv, sh in zip(m, shift):
        hmax = np.ceil(smax / astar + np.abs(sh)).astype(int) + 1
        hkl = np.stack(np.meshgrid(*[np.arange(-x, x + 1) for x in hmax], indexing="ij"), -1).reshape(-1, 3)
        S = hkl + sh
        s2 = np.einsum("ni,ij,nj->n", S, gstar, S)
        ok = (s2 <= smax**2) & (s2 > 1e-12)
        if seen is not None:
            for i in np.flatnonzero(ok):
                key = tuple(np.rint(S[i] * 1e6).astype(np.int64))
                if key in seen:
                    ok[i] = False
                else:
                    seen.add(key)
        out.append(np.hstack([hkl[ok], np.broadcast_to(mv, (ok.sum(), model.d))]))
    return np.vstack(out) if out else np.zeros((0, 3 + model.d))


def _canonical(model, H):
    """Per reflection, the key of its symmetry-equivalent set (the smallest of its images H W, packed into an
    integer), or None when the indices are too large to pack. Symmetry-equivalent reflections have equal |F|, with
    anomalous scattering too (an operation of the structure maps its density onto itself)."""
    Hi = np.rint(H).astype(np.int64)
    bits = 60 // Hi.shape[1]
    off = 1 << (bits - 1)
    images = np.stack([Hi @ np.rint(W).astype(np.int64) for W, _ in model.ops])  # (ops, n, 3+d)
    if np.abs(images).max(initial=0) >= off:
        return None
    radix = np.left_shift(np.int64(1), bits * np.arange(Hi.shape[1] - 1, -1, -1, dtype=np.int64))
    return ((images + off) @ radix).min(axis=0)


def _intensities_f2(model, H, wavelength):
    """|F|^2 for superspace indices H, each set of symmetry-equivalent reflections calculated once (and remembered
    for the Cu Ka1 wavelength across calls)."""
    keys = _canonical(model, H)
    if keys is None:
        return np.abs(structure_factors(model, H, wavelength)) ** 2
    uniq, first, inverse = np.unique(keys, return_index=True, return_inverse=True)
    cache = model.f2 if wavelength == pattern.CU_KA1 else {}
    todo = [j for j, k in enumerate(uniq) if int(k) not in cache]
    if todo:
        f2 = np.abs(structure_factors(model, H[first[todo]], wavelength)) ** 2
        cache.update(zip((int(uniq[j]) for j in todo), f2))
    return np.array([cache[int(k)] for k in uniq])[inverse.ravel()]


def _commensurate_order(model):
    """Largest satellite order needed to reach every distinct scattering vector of a commensurate structure."""
    return max(Fraction(float(x)).limit_denominator(1000).denominator for x in model.q.ravel())


def calculate(model, basis, cell, dmin, wavelength=pattern.CU_KA1):
    """All lines with d >= dmin, as pattern.calculate, in the setting of cell (the stored structure's) reached from
    the CIF's by x' = basis x + t. Satellite orders are added until two in a row give no reflection of
    SATELLITE_LIMIT on the 0-1000 scale. Returns (lines, imax, highest order used, strongest reflection beyond
    it if MAX_ORDER stopped the search, else 0)."""
    binv = np.linalg.inv(basis)
    gstar = np.array(cell.reciprocal_metric_tensor().as_mat33().tolist())
    commensurate = model.t0 is not None
    seen = set() if commensurate else None
    last = _commensurate_order(model) if commensurate else MAX_ORDER[model.d]
    hs, inten = [], []
    top, quiet, left = 0.0, 0, 0.0
    order = 0
    dlim = 0.0  # largest d of the main and first-order reflections
    for order in range(last + 1):
        H = _reflections(model, _order_indices(model.d, order), dmin, seen)
        if len(H):
            S = (H[:, :3] + H[:, 3:] @ model.q) @ binv
            d = 1 / np.sqrt(np.einsum("ni,ij,nj->n", S, gstar, S))
            if order > 1 and not commensurate:
                # a satellite of high order m lies near the origin where m q comes close to a lattice vector
                # (mullite, COD 2108042: 10 q ~ (3 0 5), d = 3800 A); such reflections, at angles no powder
                # measurement reaches, are left out: none may have a larger d than the main and first-order ones
                keep = d <= dlim * (1 + 1e-9)
                H, S, d = H[keep], S[keep], d[keep]
        if len(H):
            theta = np.arcsin(wavelength / (2 * d))
            lp = (1 + np.cos(2 * theta) ** 2) / (np.sin(theta) ** 2 * np.cos(theta))
            i = _intensities_f2(model, H, wavelength) * lp
            hs.append((H, S, d))
            inten.append(i)
            strongest = i.max()
            if order <= 1:
                dlim = max(dlim, float(d.max()))
        else:
            strongest = 0.0
        top = max(top, strongest) if order == 0 else top
        if order > 0 and not commensurate:
            quiet = quiet + 1 if strongest < SATELLITE_LIMIT / 1000 * top else 0
            if quiet == 2:
                break
    else:
        if not commensurate and quiet < 2:
            left = strongest / top * 1000
    H = np.vstack([h for h, _, _ in hs])
    S = np.vstack([s for _, s, _ in hs])
    d = np.concatenate([x for _, _, x in hs])
    i = np.concatenate(inten)
    present = i > 1e-12 * i.max()  # e.g. reflections a centred basic cell forbids, stored in its primitive cell
    H, S, d, i = H[present], S[present], d[present], i[present]
    hkl = np.rint(S).astype(np.int64)
    main = ~H[:, 3:].any(axis=1) if not commensurate else np.all(np.abs(S - hkl) < 1e-6, axis=1)
    if not np.allclose(S[main], hkl[main], atol=1e-6):
        raise ModulationError("main reflections are not integral in the stored setting")
    hkl[~main] = 0  # satellites have no index in the basic cell
    lines, imax = pattern.merge_lines(hkl, d, i, wavelength)
    weak = ~main & (i < SATELLITE_LIMIT / 1000 * imax)
    if weak.any() and not commensurate:  # a commensurate structure is a supercell: all its lines count
        lines, imax = pattern.merge_lines(hkl[~weak], d[~weak], i[~weak], wavelength)
    return lines, imax, order, left


def stick_pattern(model, basis, cell, wavelength=pattern.CU_KA1):
    """Lines to store, as pattern.stick_pattern: up to 90 deg 2-theta (140 deg if fewer than 10), at most
    pattern.MAX_LINES by decreasing d. Returns (lines, imax, highest satellite order, strongest reflection
    beyond it if MAX_ORDER was reached)."""
    d90 = pattern.d_for_two_theta(90.0, wavelength)
    dmin = max(d90, (4 * math.pi * cell.volume / (3 * 4 * pattern.MAX_LINES)) ** (1 / 3))
    while True:
        lines, imax, order, left = calculate(model, basis, cell, dmin, wavelength)
        if len(lines) >= pattern.MAX_LINES or dmin <= d90:
            break
        dmin = max(d90, dmin * 0.8)
    if len(lines) < 10:
        lines, imax, order, left = calculate(model, basis, cell, pattern.d_for_two_theta(140.0, wavelength),
                                             wavelength)
    lines = lines[:pattern.MAX_LINES]
    top = max((x.intensity for x in lines), default=0)
    if 0 < top < 1000:
        imax *= top / 1000
        lines = [pattern.Line(x.d, 1000 * x.intensity / top, x.hkl) for x in lines]
    return lines, imax, order, left


def agreement(model, block):
    """How well the structure factors reproduce those the refinement program lists in the CIF (JANA writes F_calc
    or F^2_calc for the observed main and satellite reflections), with the CIF's own f' and f'', neutron
    scattering lengths or electron form factors, and twin domains: R = sum |k F2 - F2_cif| / sum F2_cif, k
    fitted, over all, main and satellite reflections. None if the CIF lists no such values (or its twin law does
    not map the wave vectors)."""
    idx = ["index_h", "index_k", "index_l"] + [f"index_m_{i + 1}" for i in range(model.d)]
    for fc in ("F_squared_calc", "F_calc"):
        t = block.find("_refln_", idx + [fc])
        if len(t):
            break
    else:
        return None
    rows = [[t[i][j] for j in range(len(idx) + 1)] for i in range(len(t))]
    rows = [r for r in rows if _num(r[-1]) is not None]
    if len(rows) < 10:
        return None
    H = np.array([[int(_num(x)) for x in r[:-1]] for r in rows])
    ref = np.array([_num(r[-1]) for r in rows]) ** (2 if fc == "F_calc" else 1)
    wl = _num(block.find_value("_diffrn_radiation_wavelength")) or 0.71073
    dispersion, lengths = {}, {}
    for r in _loop(block, "_atom_type_", ["symbol", "?scat_dispersion_real", "?scat_dispersion_imag",
                                                      "?scat_length_neutron"]):
        el = re.match(r"[A-Za-z]+", gemmi.cif.as_string(r[0])).group(0).capitalize()
        if _num(r[1]) is not None:
            dispersion[el] = complex(_num(r[1]), _num(r[2]) or 0.0)
        if _num(r[3]) is not None:
            lengths[el] = _num(r[3])
    probe = gemmi.cif.as_string(block.find_value("_diffrn_radiation_probe") or "").strip().lower()
    if probe == "electron":  # kinematic electron form factors (International Tables C, Table 4.3.2.2)
        els = {a.element for a in model.atoms}
        coefs = {el: np.array(gemmi.Element(el).c4322.get_coefs()) for el in els}
        lengths = {el: (lambda s2, c=c: np.exp(-np.outer(s2, c[5:])) @ c[:5]) for el, c in coefs.items()}
    twins = _loop(block, "_twin_individual_", ["mass_fraction_refined"] + [
        f"twin_matrix_{i}{j}" for i in (1, 2, 3) for j in (1, 2, 3)])
    twins = [[_num(x) for x in r] for r in twins if all(_num(x) is not None for x in r)] or [[1, 1, 0, 0, 0, 1, 0, 0, 0, 1]]
    F2 = np.zeros(len(H))
    for r in twins:
        T = np.array(r[1:]).reshape(3, 3)
        M = np.linalg.lstsq(model.q.T, (model.q @ T).T, rcond=None)[0].T  # the twin maps each q onto a combination
        if not np.allclose(np.rint(M) @ model.q, model.q @ T, atol=1e-3):
            return None
        Ht = np.hstack([H[:, :3] @ T, H[:, 3:] @ np.rint(M)])
        F2 += r[0] * np.abs(structure_factors(model, Ht, wl, lengths or None, dispersion)) ** 2
    ok = ref > 0
    k = (ref[ok] * F2[ok]).sum() / (F2[ok] ** 2).sum()
    sat = H[:, 3:].any(axis=1)

    def R(s):
        return float(np.abs(k * F2[s] - ref[s]).sum() / ref[s].sum()) if ref[s].sum() > 0 else float("nan")

    return {"n": len(H), "nsat": int(sat.sum()), "scale": round(float(k), 4), "R": round(R(ok), 4),
            "Rmain": round(R(ok & ~sat), 4), "Rsat": round(R(ok & sat), 4),
            "probe": probe if probe in ("electron", "neutron") else "neutron" if lengths else "xray",
            "twins": len(twins)}
