"""Build HighScore phase records (crystal structures for Rietveld use) from CIF data."""

import json
import math
import re
import struct
from functools import lru_cache
from pathlib import Path

import gemmi
import numpy as np
import spglib
import spglib.error

from . import phase
from . import sites as sites_mod

spglib.error.OLD_ERROR_HANDLING = False  # raise errors instead of returning None

EIGHT_PI2 = 8 * math.pi**2
ISSUE_OFFSET = 24


TEMPLATE_FILES = {phase.HS4.version: "templates.json", phase.HS3.version: "templates_hs3.json"}


@lru_cache(maxsize=None)
def _templates(version=phase.HS4.version):
    return json.loads(Path(__file__).with_name(TEMPLATE_FILES[version]).read_text())


def _trv(kind, name, value=None, esd=0.0, version=phase.HS4.version):
    refine, flags, _, values = _templates(version)["trv"][kind]
    values = list(values)
    if value is not None:
        values[0] = value
    if values[1] is None:
        values[1] = esd
    return phase.TRV(name, tuple(values), refine, bytes.fromhex(flags))


def _num(block, tag):
    """Value and standard uncertainty of a numeric CIF item; None when absent or not a number."""
    v = block.find_value(tag)
    return None if v is None else _num_str(v)


_NUM = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)(?:\((\d+)\)?)?")


def _num_str(v):
    """Value and standard uncertainty from a CIF number such as '0.1234(5)'; None if not a number."""
    m = _NUM.match(v.strip("'\""))
    if not m:
        return None
    body, su = m.group(1), m.group(2)
    mantissa, _, exponent = body.lower().partition("e")
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    return float(body), (int(su) * 10 ** (int(exponent or 0) - decimals) if su else 0.0)


def _text(block, tag):
    """CIF text as HighScore keeps it: line breaks become spaces."""
    v = block.find_value(tag)
    if v is None or v in ("?", "."):
        return ""
    return gemmi.cif.as_string(v).replace("\r", "").replace("\n", " ")


def _int(block, tag):
    """Integer CIF item, 0 when absent, not a plain number or too large for an int32 field."""
    t = _text(block, tag).strip()
    return int(t) if t.isdigit() and int(t) < 2**31 else 0


def site_info(st, images):
    """Per asymmetric-unit site: (multiplicity, Wyckoff symbol). Independent of the record version.

    images is sites.unit_cell_images(st).
    """
    letters = _wyckoff_letters(st, images)
    return [(int(m), f"{m}{letters[i]}" if i in letters else "")
            for i, m in enumerate(images.multiplicities(len(st.sites)))]


def _wyckoff_letters(st, images):
    """Wyckoff letter per asymmetric-unit site index, from spglib.

    Sites are first told apart by element, occupancy and site; disordered structures can then show lower
    symmetry, so element-only typing is tried next, each at a strict and then a looser tolerance. Sites stay
    unassigned if no attempt finds the CIF's group.
    """
    lattice = np.array(st.cell.orth.mat.tolist()).T
    attempts = [(kind_of, symprec)
                for symprec in (1e-3, 1e-2)  # Å; the looser pass covers coordinates like 0.333 for 1/3
                for kind_of in (lambda i: (st.sites[i].element.name, round(st.sites[i].occ, 3), i),
                                lambda i: st.sites[i].element.name)]
    for kind_of, symprec in attempts:
        kinds = {}
        types = [kinds.setdefault(kind_of(int(i)), len(kinds)) for i in images.site]
        try:
            ds = spglib.get_symmetry_dataset((lattice, images.xyz, types), symprec=symprec)
        except spglib.error.SpglibError:  # spglib rejects some degenerate cells
            continue
        if ds.number == st.spacegroup.number:
            letters = {}
            for i, w in zip(images.site, ds.wyckoffs):
                letters.setdefault(int(i), w)
            return letters
    return {}


def _charge(type_symbol):
    """Charge from a CIF type symbol such as 'O2-', 'O-2', 'Fe3+' or 'Na+'."""
    m = re.search(r"(\d*)([+-])(\d*)$", type_symbol or "")
    if not m:
        return 0
    n = int(m.group(1) or m.group(3) or 1)
    return n if m.group(2) == "+" else -n


def _clamp(kind, v, version):
    """HighScore keeps refinable values inside the parameter's limits (which differ between versions)."""
    lo, hi = _templates(version)["trv"][kind][3][3:5]
    return min(max(v, lo), hi) if hi > lo else v


def build(entry, fmt=phase.HS4, sites=None):
    """Phase record for an entry built in this process (entry.structure is in the reference setting).

    Values come from the cleaned structure; standard uncertainties from the matching CIF rows, and for
    anisotropic tensors only where the setting change leaves them valid. fmt selects the record version:
    phase.HS4 (HighScore 4.x and later) or phase.HS3 (HighScore 3.x). sites is site_info(), passed in when
    building several versions of the same structure.
    """
    st, block = entry.structure, entry.block
    if st is None:
        raise ValueError("entry has no structure (it was built in another process)")
    v = fmt.version
    tpl = _templates(v)
    prefix = entry.compound_name + " "
    c = st.cell
    cell = [_trv(name, name, val, version=v)
            for name, val in zip(phase.CELL_NAMES, (c.a, c.b, c.c, c.alpha, c.beta, c.gamma))]

    # CIF atom_site rows in file order; entry.site_rows maps each structure site to its row
    rows = list(block.find("_atom_site_", ["label", "?U_iso_or_equiv", "?B_iso_or_equiv", "?occupancy",
                                           "?calc_flag"]))
    aniso_rows = list(block.find("_atom_site_aniso_", ["label", "U_11", "U_22", "U_33", "U_12", "U_13", "U_23"]))
    aniso_labels = [r[0] for r in aniso_rows]
    sites = sites or site_info(st, entry.images)

    def su(row, col):
        parsed = _num_str(row[col]) if row is not None and row.has(col) else None
        return parsed[1] if parsed else 0.0

    atoms = []
    for k, (site, (mult, wyckoff), r) in enumerate(zip(st.sites, sites, entry.site_rows, strict=True)):
        if k in entry.shell_sites:  # a shell scatterer has no representation among atoms
            continue
        row = rows[r] if r is not None and r < len(rows) and rows[r][0] == site.label else None
        pos = [_trv("ATOM " + axis, f"{site.label} {axis}", x, version=v)
               for axis, x in zip(phase.ATOM_POS, (site.fract.x, site.fract.y, site.fract.z))]
        u_given = row is not None and row.has(1) and _num_str(row[1]) is not None
        if site.u_iso > 0:  # U_iso_or_equiv, or B_iso_or_equiv converted by gemmi
            biso = site.u_iso * EIGHT_PI2
            bsu = su(row, 1) * EIGHT_PI2 if u_given else su(row, 2)
        elif site.aniso.nonzero():  # equivalent isotropic value of the tensor
            biso, bsu = EIGHT_PI2 * np.trace(sites_mod.ucif_to_cart(site.aniso, c)) / 3, 0.0
        else:  # stated as 0, or negative and set to 0 (entry._default_b gave atoms without one the default)
            biso, bsu = 0.0, 0.0
        u = site.aniso
        uij = (u.u11, u.u22, u.u33, u.u12, u.u13, u.u23) if u.nonzero() else (0.0,) * 6
        uij_su = (0.0,) * 6
        if u.nonzero() and not entry.setting_changed and aniso_labels.count(site.label) == 1:
            arow = aniso_rows[aniso_labels.index(site.label)]
            uij_su = tuple(su(arow, k) for k in range(1, 7))
        aniso = [_trv("ATOM " + comp, f"{site.label} {comp}", _clamp("ATOM " + comp, x * EIGHT_PI2, v),
                      e * EIGHT_PI2, version=v)
                 for comp, x, e in zip(phase.ATOM_ANISO, uij, uij_su)]
        calc = row[4].strip("'\"").lower() if row is not None and row.has(4) else ""
        tail = bytearray(26)
        tail[20] = 1 if calc in ("calc", "c") else 2 if calc == "dum" else 0
        atoms.append(phase.Atom(
            label=site.label, element=site.element.name,
            pos=pos,
            biso=_trv("ATOM Biso", f"{site.label} Biso", _clamp("ATOM Biso", biso, v), bsu, version=v),
            sof=_trv("ATOM sof", f"{site.label} sof", _clamp("ATOM sof", site.occ, v), su(row, 3), version=v),
            aniso=aniso,
            multiplicity=mult,
            wyckoff=wyckoff,
            charge=_charge(site.type_symbol),
            tail=bytes(tail),
        ))

    laue = tpl["laue"][st.spacegroup.laue_str()]
    harmonics = [phase.TRV(n, tuple(v[3]), v[0], bytes.fromhex(v[1])) for n, v in laue["harmonics"]]

    r_all = (_num(block, "_refine_ls_R_factor_all") or (0.0, 0.0))[0]
    rest = bytearray.fromhex(tpl["trailer_rest"])
    struct.pack_into("<i", rest, ISSUE_OFFSET, _int(block, "_journal_issue"))
    comment = []
    if _text(block, "_exptl_crystal_colour"):
        comment.append("Crystal color: " + _text(block, "_exptl_crystal_colour"))
    if _text(block, "_exptl_crystal_description"):
        comment.append("Crystal description: " + _text(block, "_exptl_crystal_description"))
    trailer = {
        "pre": bytes(16), "r_factor": r_all, "name": entry.compound_name, "source": fmt.source,
        "cod_id": entry.cod_id, "i0": 0,
        "date": _text(block, "_audit_creation_date"), "method": _text(block, "_audit_creation_method"),
        "s1": "", "s2": "",
        "systematic": _text(block, "_chemical_name_systematic"), "mineral": _text(block, "_chemical_name_mineral"),
        "common": _text(block, "_chemical_name_common"), "formula_sum": _text(block, "_chemical_formula_sum"),
        "formula_struct": _text(block, "_chemical_formula_structural"), "title": _text(block, "_publ_section_title"),
        "authors": ", ".join(gemmi.cif.as_string(a) for a in block.find_values("_publ_author_name")),
        "s3": "", "journal": _text(block, "_journal_name_full"), "coden": _text(block, "_journal_coden_ASTM"),
        "volume": _int(block, "_journal_volume"), "year": _int(block, "_journal_year"),
        "page_first": _int(block, "_journal_page_first"), "page_last": _int(block, "_journal_page_last"),
        "s4": "", "comment": "\r\n".join(comment), "rest": bytes(rest),
    }
    return phase.Phase(
        prefix=prefix,
        scale=_trv("Scale Factor", "Scale Factor", version=v),
        cell=cell,
        crystal_system=phase.CRYSTAL_SYSTEM_CODES[st.spacegroup.crystal_system_str()],
        atoms=atoms,
        after_atoms=bytes.fromhex(tpl["after_atoms"]),
        profile=[_trv(n, n, version=v) for n in tpl["profile_names"]],
        mid=bytes.fromhex(tpl["mid"]),
        corrections=[_trv(n, n, version=v) for n in tpl["correction_names"]],
        harmonics_head=bytes.fromhex(laue["head"]),
        harmonics=harmonics,
        trailer=trailer,
        fmt=fmt,
    )

