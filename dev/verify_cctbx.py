"""Independent check of calculated patterns against cctbx.

Runs in the cctbx environment of the dev image (docker build --target dev -t hsrdb-dev .):
  docker run --rm -v $PWD:/w hsrdb-dev /opt/cctbx/bin/python dev/verify_cctbx.py LIST [-j N] [--db DB.hsrdb] > out.tsv
LIST is a file with one CIF path per line, or CIF files and folders.

Per CIF, one TSV row with:
  calc_*  our per-reflection powder intensities against cctbx for the same structure (our cleaned sites in the
          reference setting, merged images at their mean position, same f' and f''), over the reflections of the
          stored lines. Checks structure factors, form factors, displacement factors, symmetry expansion,
          multiplicities and Lorentz-polarisation; any difference is a bug. calc_rel: largest relative difference
          over reflections >= 1% of the strongest; calc_abs: largest difference over all reflections, relative to
          the strongest. friedel_rel: same measure between |F(h)|^2 and the Friedel average.
  cif_*   our lines against lines cctbx computes from the raw CIF with its own parser and symmetry, in the CIF's
          setting, with the same f', f'', neutral-atom IT92 form factors, and the deliberate choices of hsrdb-tools
          applied to cctbx's reading (our_choices; choices = 1 if applied). Checks parsing, setting change and site
          cleaning. cif_dI: largest difference of matched lines on the 0-1000 scale; cif_miss_*: strongest line
          without a partner within 0.006 deg 2-theta on each side; both without the last stored line, which ends
          at the range limit or the line cap where the two sides cut a chain of merged reflections differently
          (cif_dI_last, for information); iic_ratio, density_ratio: ours / cctbx; cif_retry: items left out so
          that cctbx could read the CIF (read_structures).
  db_dI   stored lines (LinesI, or Lines for HighScore 3.x) against the cctbx lines, when --db is given.
  mult_bad: atoms whose site-symmetry order is not one the CIF states (entry._site_orders); split_sites /
          occ_merged: sites with images between sites.SAME_POSITION and SPECIAL_POSITION_TOLERANCE apart, kept as
          a split atom / merged (stated by the CIF or by the occupancy rule; our_choices merges the same sites in
          cctbx's reading); special_shift: largest move (A) cctbx made to put one of our
          sites on its special position; dummy_sites: atoms flagged 'dum'; reread: sites whose element differs from
          gemmi's reading of the CIF (cctbx reads the CIF's spelling); corrected: built with a correction from
          corrections.json; formula_dev: largest relative difference of a non-H element count from
          _chemical_formula_sum x Z.
"""

import argparse
import math
import types
import re
import sqlite3
import struct
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi
import numpy as np
from cctbx import adptbx, crystal, miller, sgtbx, xray
from cctbx.array_family import flex
from cctbx.eltbx import xray_scattering
from iotbx import cif as icif
from iotbx.cif import builders

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import build, entry, pattern, sites  # noqa: E402

WL = pattern.CU_KA1
MATCH_TWO_THETA = 0.006
COLUMNS = ["path", "cod_id", "status", "n_refl", "calc_rel", "calc_abs", "friedel_rel", "calc_mult_bad",
           "special_shift", "cif_status", "cif_lines", "cif_dI", "cif_d2t", "cif_miss_ours", "cif_miss_cctbx",
           "iic_ratio", "density_ratio", "db_dI", "mult_bad", "split_sites", "formula_dev", "invalid_adps",
           "dropped_sites", "setting_changed", "occ_merged", "dummy_sites", "choices", "reread", "corrected",
           "cif_dI_last", "cif_retry", "modulated"]


def two_theta(d):
    return np.degrees(2 * np.arcsin(WL / (2 * np.asarray(d))))


def lp(d):
    th = np.arcsin(WL / (2 * np.asarray(d)))
    return (1 + np.cos(2 * th) ** 2) / (np.sin(th) ** 2 * np.cos(th))


def anomalous(el):
    return pattern._anomalous(el, WL)


def neutral(scattering_type):
    el, _ = xray_scattering.get_element_and_charge_symbols(scattering_type, exact=False)
    return el.capitalize()


def set_scattering(xs):
    for sc in xs.scatterers():
        sc.scattering_type = neutral(sc.scattering_type)
        sc.fp, sc.fdp = anomalous(sc.scattering_type)
    xs.scattering_type_registry(table="it1992")


def cif_group(block, e):
    """cctbx space group of the symmetry statement hsrdb-tools used, in the CIF's own setting: the listed
    operations, else the Hall symbol, parsed by cctbx; the H-M symbol when the operations contradict the CIF's
    composition (e.symmetry_note). None to keep cctbx's own reading (only an H-M symbol, or a correction)."""
    if e.correction and {"spacegroup", "operations"} & set(entry._corrections().get(str(e.cod_id), {})):
        return None
    if e.symmetry_note:
        return sgtbx.space_group("Hall: " + entry._symbol_group(block).hall)
    ops = [gemmi.cif.as_string(v) for v in (block.find_values("_space_group_symop_operation_xyz")
                                            or block.find_values("_symmetry_equiv_pos_as_xyz"))]
    if ops:
        g = sgtbx.space_group()
        for op in ops:
            g.expand_smx(sgtbx.rt_mx(gemmi.Op(op).triplet()))
        return g
    for tag in ("_space_group_name_Hall", "_symmetry_space_group_name_Hall"):
        hall = entry._text(block, tag)
        if hall:
            return sgtbx.space_group_info("Hall: " + hall).group()
    return None


SYMBOL_TAGS = ("_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt", "_symmetry_space_group_name_Hall",
               "_space_group_name_Hall", "_symmetry_Int_Tables_number", "_space_group_IT_number")
WAVELENGTH_TAGS = ("_diffrn_radiation_wavelength",)


def read_structures(path):
    """cctbx's reading of the CIF, and what had to be left out for cctbx to read it: the space-group symbols when
    they disagree with the listed operations (cctbx refuses; cif_group then sets the group hsrdb-tools used), and a
    wavelength that is not a number (cctbx parses it; the calculation does not use it)."""
    reader = icif.reader(file_path=path)
    try:
        return reader.build_crystal_structures(), ""
    except Exception as exc:  # retried once below; the original error is raised if the retry fails too
        error, message = exc, str(exc)
    drop = (SYMBOL_TAGS if "Inconsistent symmetry" in message else ()) + \
        (WAVELENGTH_TAGS if "could not convert string to float" in message else ())
    removed = set()
    for block in reader.model().values():
        for tag in drop:
            if tag in block:
                del block[tag]
                removed.add("symbol" if tag in SYMBOL_TAGS else "wavelength")
    if not removed:
        raise error
    try:
        return reader.build_crystal_structures(), ",".join(sorted(removed))
    except Exception:
        raise error


def our_choices(xs, path, e):
    """cctbx's reading of the CIF with the deliberate choices of hsrdb-tools applied, so that they explain no
    difference: aniso rows as given (beta_ij too, which cctbx does not read; repeated labels paired in order, which
    cctbx gets wrong); dummy atoms (calc_flag dum) left out; tensors with all components zero or a negative
    mean-square displacement (judged before site symmetry is imposed) replaced by the isotropic value; B = sites.DEFAULT_B for atoms without any displacement parameter
    (cctbx leaves them at U = 0); images of an atom merged within the distance hsrdb-tools chose for it
    (sites.merge_tolerances), with cctbx's own special-position handling; the space group of the symmetry
    statement hsrdb-tools used (cif_group); the cell constrained to the crystal system when the setting is
    unchanged. Needs one scatterer per atom_site row, in row order; returns
    (structure, applied)."""
    block = gemmi.cif.read(path).sole_block()
    if len(block.find_values("_atom_site_label")) != xs.scatterers().size():
        return xs, False
    tol = {}
    for (t, _), r in zip(sites.merge_tolerances(e.structure, None, e.site_orders, e.separate_stated), e.site_rows):
        if r is not None:
            tol[r] = max(t, tol.get(r, 0.0))
    u_col, b_col, occ_col = (block.find_values(t) for t in ("_atom_site_U_iso_or_equiv", "_atom_site_B_iso_or_equiv",
                                                           "_atom_site_occupancy"))
    flags = [gemmi.cif.as_string(v).strip().lower() for v in block.find_values("_atom_site_calc_flag")]
    uc = xs.unit_cell()
    # the aniso rows as the CIF gives them, converted with cctbx: U_ij, B_ij, and beta_ij (U* = beta / 2 pi^2),
    # which cctbx's reader ignores; a label repeated with one row per atom (disorder parts) pairs the rows in order,
    # where cctbx gives every atom the first row (plans/upstream note); None for a row that is not all numbers
    tensors = {}
    for kind in ("U", "B", "beta"):
        table = block.find("_atom_site_aniso_", ["label"] + [f"{kind}_{ij}" for ij in ("11", "22", "33", "12", "13", "23")])
        for r in table:
            v = [gemmi.cif.as_number(r[k]) for k in range(1, 7)]
            if not all(map(math.isfinite, v)):
                raw = None
            elif kind == "beta":
                raw = tuple(x / (2 * math.pi**2) for x in v)
            else:
                raw = adptbx.u_cif_as_u_star(uc, tuple(x / (8 * math.pi**2) if kind == "B" else x for x in v))
            tensors.setdefault(r.str(0), []).append(raw)
        if tensors:
            break
    labels = [sc.label for sc in xs.scatterers()]
    adp_types = [gemmi.cif.as_string(v).strip().lower() for v in
                 block.find_values("_atom_site_adp_type") or block.find_values("_atom_site_thermal_displace_type")]
    keep = flex.bool()
    for row, sc in enumerate(xs.scatterers()):
        keep.append(not (row < len(flags) and flags[row] == "dum"))
        # cctbx drops a whole occupancy or U column when one value is '?' or '.'; those values default to 1 and
        # (below) to U_iso_or_equiv, B_iso_or_equiv or B = sites.DEFAULT_B
        occ = gemmi.cif.as_number(occ_col[row]) if occ_col and row < len(occ_col) else math.nan
        sc.occupancy = occ if math.isfinite(occ) else 1.0
        rows = tensors.get(sc.label, [])
        if rows:
            same = [r for r, lab in enumerate(labels) if lab == sc.label]
            aniso = [r for r in same if r < len(adp_types) and adp_types[r] in ("uani", "bani")]
            if len(rows) == len(same):  # one row per atom (disorder parts), in order
                raw = rows[same.index(row)]
            elif len(rows) == len(aniso):  # one row per atom the CIF marks anisotropic
                raw = rows[aniso.index(row)] if row in aniso else None
            else:  # one row for the atoms of a mixed site, at the position of the first
                f = [(a - b + 0.5) % 1 - 0.5 for a, b in zip(sc.site, xs.scatterers()[same[0]].site)]
                raw = rows[0] if uc.length(f) < 1e-3 else None
            sc.flags.set_use_u_iso(raw is None)
            sc.flags.set_use_u_aniso(raw is not None)
            sc.u_star = raw if raw is not None else (0, 0, 0, 0, 0, 0)
        u = gemmi.cif.as_number(u_col[row]) if u_col and row < len(u_col) else math.nan
        b = gemmi.cif.as_number(b_col[row]) if b_col and row < len(b_col) else math.nan
        # isotropic value: U_iso, else B_iso, else trace / 3 of an invalid tensor (>= entry.MIN_UEQ), else default
        iso = u if math.isfinite(u) else max(b, 0.0) / (8 * math.pi**2) if math.isfinite(b) else math.nan
        if sc.flags.use_u_aniso():
            ucart = np.array(adptbx.u_star_as_u_cart(uc, sc.u_star))
            m = ucart[[0, 3, 4, 3, 1, 5, 4, 5, 2]].reshape(3, 3)
            if not ucart.any() or np.linalg.eigvalsh(m).min() < -entry.ADP_MARGIN:
                sc.flags.set_use_u_aniso(False)
                sc.flags.set_use_u_iso(True)
                sc.u_star = (0, 0, 0, 0, 0, 0)
                if not math.isfinite(iso) and (ucart[0] + ucart[1] + ucart[2]) / 3 >= entry.MIN_UEQ:
                    iso = (ucart[0] + ucart[1] + ucart[2]) / 3
        if not sc.flags.use_u_aniso():
            sc.u_iso = max(iso, 0.0) if math.isfinite(iso) else sites.DEFAULT_B / (8 * math.pi**2)
    cs = xs.crystal_symmetry()
    group = cif_group(block, e)
    if group is not None:
        cs = crystal.symmetry(unit_cell=cs.unit_cell(), space_group=group)
    if not e.setting_changed:  # the cell as constrained to the crystal system (setting.constrain_cell), U_cif kept
        cs = crystal.symmetry(unit_cell=e.cell, space_group_info=cs.space_group_info())
    out = xray.structure(special_position_settings=cs.special_position_settings(
        min_distance_sym_equiv=sites.SAME_POSITION))
    raw = gemmi.make_small_structure_from_block(block).sites
    for row, sc in enumerate(xs.scatterers()):
        if keep[row]:
            if sc.flags.use_u_aniso():
                sc.u_star = adptbx.u_cif_as_u_star(cs.unit_cell(), adptbx.u_star_as_u_cif(uc, sc.u_star))
            if tol.get(row) == 0.0:  # images kept apart: the CIF's coordinates, which cctbx's reader moved
                sc.site = (raw[row].fract.x, raw[row].fract.y, raw[row].fract.z)
            sps = cs.special_position_settings(min_distance_sym_equiv=tol.get(row, sites.SAME_POSITION))
            out.add_scatterer(sc, sps.site_symmetry(sc.site))
    out.scattering_type_registry(table="it1992")  # a new structure starts with cctbx's default table
    return out, True


def powder(xs, indices, shells=()):
    """cctbx: multiplicity * Friedel-averaged |F|^2 * LP, and |F(h)|^2 alone, for the given indices. shells:
    sites.Shell scatterers, which cctbx has no form factor for: their term sum(n f(q)) * sin(qR)/(qR) *
    exp(-2 pi^2 U / d^2) at each centre, with cctbx's IT92 form factors, is added to cctbx's structure factors."""
    cs = xs.crystal_symmetry()
    idx = flex.miller_index([tuple(int(v) for v in h) for h in indices])
    neg = flex.miller_index([tuple(-int(v) for v in h) for h in indices])
    fh = np.array(miller.set(cs, idx, anomalous_flag=True).structure_factors_from_scatterers(
        xs, algorithm="direct").f_calc().data())
    fm = np.array(miller.set(cs, neg, anomalous_flag=True).structure_factors_from_scatterers(
        xs, algorithm="direct").f_calc().data())
    d = np.array(miller.set(cs, idx).d_spacings().data())
    h = np.asarray(indices, dtype=float)
    for sh in shells:
        x = 2 * math.pi * sh.radius / d
        stol = 1 / (2 * d)
        atoms = sum(n / len(sh.xyz) * (np.array([xray_scattering.it1992(el).fetch().at_stol(v) for v in stol])
                                       + complex(*anomalous(el))) for el, n in sh.per_cell.items())
        term = atoms * np.sin(x) / x * np.exp(-2 * math.pi**2 * sh.u_iso / d**2)
        phase = np.exp(2j * math.pi * (h @ sh.xyz.T)).sum(axis=1)
        fh, fm = fh + term * phase, fm + term * phase.conj()
    ih = np.abs(fh) ** 2
    im = np.abs(fm) ** 2
    mult = np.array(miller.set(cs, idx, anomalous_flag=False).multiplicities().data())
    return mult * (ih + im) / 2 * lp(d), mult * ih * lp(d), mult, d


def merge(d, inten):
    """Our line rule: a reflection within pattern.MERGE_TWO_THETA of the previous one joins its line; reflections
    of zero intensity are left out."""
    present = inten > 1e-12 * inten.max()
    d, inten = d[present], inten[present]
    t = two_theta(d)
    out = []  # [last 2theta, total, strongest, d]
    for i in np.argsort(t, kind="stable"):
        if out and t[i] - out[-1][0] <= pattern.MERGE_TWO_THETA:
            m = out[-1]
            m[0] = t[i]
            m[1] += inten[i]
            if inten[i] > m[2] * (1 + 1e-9):
                m[2], m[3] = inten[i], d[i]
        else:
            out.append([t[i], inten[i], inten[i], d[i]])
    return np.array([m[3] for m in out]), np.array([m[1] for m in out])


def compare_lines(d_a, i_a, d_b, i_b):
    """Largest intensity and 2-theta difference of matched lines; strongest unmatched line on each side."""
    ta, tb = two_theta(d_a), two_theta(d_b)
    order = np.argsort(tb)
    tb_s = tb[order]
    used = np.zeros(len(tb), bool)
    d_i = d_t = miss_a = 0.0
    for k in range(len(ta)):
        j = np.searchsorted(tb_s, ta[k])
        cand = [c for c in (j - 1, j) if 0 <= c < len(tb_s)]
        best = min(cand, key=lambda c: abs(tb_s[c] - ta[k]), default=None)
        if best is not None and abs(tb_s[best] - ta[k]) <= MATCH_TWO_THETA:
            used[order[best]] = True
            d_i = max(d_i, abs(i_a[k] - i_b[order[best]]))
            d_t = max(d_t, abs(tb_s[best] - ta[k]))
        else:
            miss_a = max(miss_a, i_a[k])
    miss_b = max((i_b[j] for j in range(len(tb)) if not used[j]), default=0.0)
    return d_i, d_t, miss_a, miss_b


def _structure_with_our_tolerance(crystal_symmetry, scatterers, wavelength=None):
    sps = crystal_symmetry.special_position_settings(min_distance_sym_equiv=sites.SAME_POSITION)
    return xray.structure(special_position_settings=sps, scatterers=scatterers, wavelength=wavelength)


# cctbx's CIF builder merges symmetry images closer than 0.5 A; ours only those closer than sites.SAME_POSITION
# (plus merges stated by the CIF or by the occupancy rule, counted in occ_merged)
builders.xray = types.SimpleNamespace(**{**vars(xray), "structure": _structure_with_our_tolerance})


def cctbx_from_ours(st, images, separate=()):
    """cctbx structure of our scattering sites. Each site enters at its first image (the identity operation), which
    for images we merged is their mean position, on the special position, with the mean of their rotated tensors:
    cctbx then finds the same site symmetry, so merged sites are compared as well. Sites in separate (indices)
    keep all their images (a general position for cctbx)."""
    cs = crystal.symmetry(unit_cell=(st.cell.a, st.cell.b, st.cell.c, st.cell.alpha, st.cell.beta, st.cell.gamma),
                          space_group_symbol="Hall: " + st.spacegroup.hall)
    uc = cs.unit_cell()
    scs = flex.xray_scatterer()
    for i in sorted(set(images.site)):
        s, k = st.sites[i], int(np.flatnonzero(images.site == i)[0])
        uk = images.ucart[k]
        u = (adptbx.u_cart_as_u_star(uc, (uk[0, 0], uk[1, 1], uk[2, 2], uk[0, 1], uk[0, 2], uk[1, 2]))
             if s.aniso.nonzero() else s.u_iso)
        scs.append(xray.scatterer(label=s.label, site=tuple(float(v) for v in images.xyz[k]), u=u, occupancy=s.occ,
                                  scattering_type=s.element.name))
    sps = cs.special_position_settings(min_distance_sym_equiv=sites.SAME_POSITION)
    xs = xray.structure(special_position_settings=sps)
    for i, sc in zip(sorted(set(images.site)), scs):
        tol = 0.0 if i in separate else sites.SAME_POSITION
        xs.add_scatterer(sc, cs.special_position_settings(min_distance_sym_equiv=tol).site_symmetry(sc.site))
    shift = max((uc.distance(tuple(a.site), tuple(b.site)) for a, b in zip(scs, xs.scatterers())), default=0.0)
    set_scattering(xs)
    return xs, shift


def cell_mass(xs):
    """Unit-cell mass of cctbx's atoms with the IUPAC weights hsrdb-tools uses (gemmi): cctbx's own table differs
    by about 2e-4, which would hide a counting error of that size."""
    return sum(sc.occupancy * sc.multiplicity() * gemmi.Element(sc.scattering_type).weight for sc in xs.scatterers())


def cctbx_lines(xs, d_last):
    """cctbx lines down to d_last (merge rule and scale as ours) and the absolute strongest intensity."""
    dmin = WL / (2 * math.sin(math.radians(two_theta(d_last) + 0.05) / 2))
    ms = miller.build_set(xs.crystal_symmetry(), anomalous_flag=False, d_min=dmin)
    inten, _, _, d = powder(xs, np.array(ms.indices()))
    ld, li = merge(d, inten)
    keep = ld >= d_last * (1 - 1e-7)
    ld, li = ld[keep], li[keep]
    imax = li.max()
    return ld, 1000 * li / imax, imax


_K_CORUNDUM = None


def k_factor(xs, imax):
    return imax / xs.unit_cell().volume() ** 2 / (cell_mass(xs) * 1.66054 / xs.unit_cell().volume())


def corundum_k():
    xs = icif.reader(input_string=entry._CORUNDUM_CIF).build_crystal_structures()["corundum"]
    set_scattering(xs)
    ms = miller.build_set(xs.crystal_symmetry(), anomalous_flag=False, d_min=pattern.d_for_two_theta(90.0))
    inten, _, _, d = powder(xs, np.array(ms.indices()))
    return k_factor(xs, merge(d, inten)[1].max())


FORMULA = re.compile(r"([A-Z][a-z]?)([0-9.]*)")


def formula_dev(block, st, images):
    fs, z = entry._text(block, "_chemical_formula_sum"), entry._number(block, "_cell_formula_units_Z")
    if not fs or not z:
        return ""
    want = {}
    for el, n in FORMULA.findall(fs):
        want[el] = want.get(el, 0) + (float(n) if n else 1.0) * z
    have = {}
    for i in images.site:
        el = st.sites[i].element.name
        have[el] = have.get(el, 0) + st.sites[i].occ
    dev = 0.0
    for el in set(want) | set(have):
        if el in ("H", "D"):
            continue
        w, h = want.get(el, 0.0), have.get(el, 0.0)
        if max(w, h) > 0:  # formulas can list an element with count 0
            dev = max(dev, abs(h - w) / max(w, h))
    return f"{dev:.4f}"


def special_stats(e):
    """Sites kept as split atoms, and sites merged although their images are more than sites.SAME_POSITION apart
    (by the CIF's statement or the occupancy rule; cctbx here merges only below SAME_POSITION, so these differ)."""
    reasons = sites.merge_tolerances(e.structure, None, e.site_orders, e.separate_stated)
    split = sum(why in ("split", "stated") and tol == sites.SAME_POSITION for tol, why in reasons)
    merged = sum(tol > sites.SAME_POSITION for tol, _ in reasons)
    return split, merged


def mult_bad(e):
    """Sites whose site-symmetry order differs from the one the CIF states (entry._site_orders, which reads the
    multiplicity in the CIF's own setting); the order does not depend on the setting."""
    if not any(e.site_orders):
        return ""
    ours = e.images.multiplicities(len(e.structure.sites))
    n_ops = len(e.structure.spacegroup.operations())
    real = set(e.scattering.site)  # dummy atoms (e.g. hydrogens of unknown position at -1,-1,-1) are no atoms
    return str(sum(stated is not None and k in real and n_ops // ours[k] not in stated
                   for k, stated in enumerate(e.site_orders)))


def reread(e):
    """Sites whose element differs from gemmi's reading of the CIF (types that are not elements, two capitals
    read again, labels without a type column): cctbx reads these as the CIF spells them."""
    raw = gemmi.make_small_structure_from_block(e.block).sites
    return sum(r is not None and r < len(raw) and raw[r].element.name != s.element.name
               for s, r in zip(e.structure.sites, e.site_rows))


def db_lines(db, cod_id):
    row = db.execute("select s.Lines, {} from general_stored s join general_indexed g on g.id = s.pid "
                     "where g.ProductID = ?".format("s.LinesI" if HAS_I else "NULL"), (f"96{cod_id + 1}",)).fetchone()
    if row is None:
        return None
    L, I = row
    n = len(L) // 6
    d = np.array([struct.unpack_from("<f", L, 6 * k)[0] for k in range(n)])
    inten = (np.array(struct.unpack(f"<{n}f", I)) if I else
             np.array([struct.unpack_from("<H", L, 6 * k + 4)[0] for k in range(n)], dtype=float))
    return d, inten


DB = None
HAS_I = False


def init(db_path):
    global DB, HAS_I, _K_CORUNDUM
    if db_path:
        DB = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        HAS_I = "LinesI" in [r[1] for r in DB.execute("pragma table_info(general_stored)")]
    _K_CORUNDUM = corundum_k()


def check(path):
    row = dict.fromkeys(COLUMNS, "")
    row["path"] = path
    try:
        try:
            e = entry.from_cif(path)
        except build.DATA_ERRORS as exc:  # what the build skips; anything else is an error row
            row["status"] = f"skipped: {type(exc).__name__}: {exc}"
            return row
        st, images = e.structure, e.scattering
        row.update(cod_id=e.cod_id, status="ok", invalid_adps=e.invalid_adps, dropped_sites=e.dropped_sites,
                   setting_changed=int(e.setting_changed),
                   dummy_sites=len(st.sites) - len(set(images.site)), corrected=int(bool(e.correction)),
                   reread=reread(e), modulated=int(e.modulated is not None))
        d_last = min(x.d for x in e.lines)

        # calc: same structure, reflections of the stored lines
        hkl = np.asarray(gemmi.make_miller_array(st.cell, st.spacegroup, d_last * (1 - 1e-7)), dtype=np.int64)
        frac = np.array(st.cell.frac.mat.tolist())
        d = 1 / np.linalg.norm(hkl @ frac, axis=1)
        ours = pattern._intensities(st, hkl, d, WL, images)
        our_mult = pattern._multiplicities(hkl, st.spacegroup)
        xs, shift = cctbx_from_ours(st, images, {k for k, w in enumerate(e.images.reasons)
                                                  if w == "stated_general" and e.separate_stated})
        c_avg, c_h, c_mult, _ = powder(xs, hkl, images.shells)
        top = c_avg.max()
        strong = c_avg >= 0.01 * top
        row["n_refl"] = len(hkl)
        row["calc_rel"] = f"{np.max(np.abs(ours - c_avg)[strong] / c_avg[strong]):.3e}"
        row["calc_abs"] = f"{np.max(np.abs(ours - c_avg)) / top:.3e}"
        row["friedel_rel"] = f"{np.max(np.abs(c_h - c_avg)[strong] / c_avg[strong]):.3e}"
        row["calc_mult_bad"] = int(np.count_nonzero(our_mult != c_mult))
        row["special_shift"] = f"{shift:.4f}"
        row["mult_bad"] = mult_bad(e)
        row["split_sites"], row["occ_merged"] = special_stats(e)
        row["formula_dev"] = formula_dev(e.block, st, images)

        # cif: cctbx reads the CIF itself
        our_d = np.array([x.d for x in e.lines])
        our_i = np.array([x.intensity for x in e.lines])
        try:
            structures, row["cif_retry"] = read_structures(path)
            if len(structures) != 1:
                raise ValueError(f"{len(structures)} structures")
            xs2 = next(iter(structures.values()))
            set_scattering(xs2)
            xs2, applied = our_choices(xs2, path, e)
            row["choices"] = int(applied)
            cd, ci, imax = cctbx_lines(xs2, d_last)
        except Exception as exc:  # cctbx could not read or use the CIF: reported, not fatal
            row["cif_status"] = f"{type(exc).__name__}: {str(exc).splitlines()[0] if str(exc) else ''}"[:200]
        else:
            row["cif_status"] = "ok"
            # the last stored line ends at the range limit (90 or 140 deg) or the line cap, where the two
            # sides cut a chain of merged reflections differently: compared separately (cif_dI_last)
            cut = two_theta(our_d[-2]) + pattern.MERGE_TWO_THETA if len(our_d) > 1 else 0.0
            body = two_theta(cd) <= cut
            di, dt, ma, mb = compare_lines(our_d[:-1], our_i[:-1], cd[body], ci[body])
            row.update(cif_lines=len(cd), cif_dI=f"{di:.3e}", cif_d2t=f"{dt:.4f}", cif_miss_ours=f"{ma:.3e}",
                       cif_miss_cctbx=f"{mb:.3e}")
            last = cd[~body]
            row["cif_dI_last"] = (f"{abs(our_i[-1] - ci[~body][np.argmin(abs(last - our_d[-1]))]):.3e}"
                                  if len(last) else f"{our_i[-1]:.3e}")
            row["iic_ratio"] = f"{e.iic / (k_factor(xs2, imax) / _K_CORUNDUM):.9f}"
            row["density_ratio"] = f"{e.density / (cell_mass(xs2) * 1.66054 / xs2.unit_cell().volume()):.9f}"
            if DB is not None:
                got = db_lines(DB, e.cod_id)
                if got is None:
                    row["db_dI"] = "missing"
                else:
                    di, _, ma, mb = compare_lines(got[0], got[1], cd, ci)
                    row["db_dI"] = f"{max(di, ma, mb):.3e}"
    except Exception:
        row["status"] = "error: " + traceback.format_exc().strip().splitlines()[-1][:200]
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("-j", type=int, default=1)
    ap.add_argument("--db")
    args = ap.parse_args()
    paths = []
    for p in map(Path, args.inputs):
        if p.is_dir():
            paths += sorted(str(x) for x in p.rglob("*.cif"))
        elif p.suffix == ".cif":
            paths.append(str(p))
        else:
            paths += [s for s in p.read_text().split() if s]
    print("\t".join(COLUMNS), flush=True)
    with ProcessPoolExecutor(args.j, initializer=init, initargs=(args.db,)) as ex:
        for row in ex.map(check, paths, chunksize=16):
            print("\t".join(str(row[c]) for c in COLUMNS), flush=True)


if __name__ == "__main__":
    main()
