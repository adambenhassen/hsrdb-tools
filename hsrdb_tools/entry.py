"""Turn a COD CIF into the fields of one reference pattern in a HighScore .hsrdb database."""

import json
import math
import os
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import gemmi
import numpy as np

from . import pattern, setting, sites, superspace

# Corundum (COD 1000032) for I/Ic
_CORUNDUM_CIF = """data_corundum
_cell_length_a 4.7605
_cell_length_b 4.7605
_cell_length_c 12.9956
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_symmetry_space_group_name_H-M 'R -3 c'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
Al1 Al 0 0 0.35216 0.0030
O1 O 0.30624 0 0.25 0.0035
"""

CRYSTAL_SYSTEM_CODES = {
    "triclinic": "A", "monoclinic": "M", "orthorhombic": "O",
    "tetragonal": "T", "trigonal": "H", "hexagonal": "H", "cubic": "C",
}
SUBFILE_INORGANIC, SUBFILE_ORGANIC, SUBFILE_MINERAL, SUBFILE_METALLIC = 0, 1, 2, 3
# elements that keep a compound out of the metallic subfile; as in the official COD databases, Ge and Sb count as
# metals and deuterium as hydrogen
_NON_METALS = {"H", "D", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar", "As", "Se", "Br", "Kr",
               "Te", "I", "Xe", "At", "Rn"}


class EntryError(Exception):
    """The CIF cannot be turned into a reference pattern."""


@dataclass
class Entry:
    cod_id: int
    crystal_system: str
    cell: tuple[float, float, float, float, float, float]
    volume: float
    sg_number: int
    sg_hm: str
    density: float
    iic: float
    elements: list[str]  # in order of first appearance
    formula: str  # unit-cell contents, e.g. "Cl4.00 C92.00 O28.00 H108.00"
    compound_name: str
    mineral_name: str
    comment: str
    literature: dict
    subfiles: list[int]
    lines: list[pattern.Line] = field(repr=False)
    # CIF atom_site row of each structure site; None for atoms added by a correction (corrections.json)
    site_rows: list[int | None] = field(repr=False, default_factory=list)
    setting_changed: bool = False
    dropped_sites: int = 0  # atom_site rows without element or coordinates
    invalid_adps: int = 0  # displacement parameters replaced because they were not physically valid
    dummy_sites: int = 0  # atoms flagged dum: kept in the structure record, not scattering
    correction: str = ""  # correction from the publication applied to this CIF (corrections.json), if any
    symbol_conflict: bool = False  # the listed operations or Hall symbol name another group than the H-M symbol
    symmetry_note: str = ""  # why the H-M symbol was used instead of the listed operations
    stated_general: bool = False  # sites on special positions that the CIF states to have lower site symmetry
    separate_stated: bool = False  # the images of those sites count separately (merge_note says why)
    merge_note: str = ""
    shell_sites: tuple = ()  # structure sites that are the centres of shell scatterers (not atoms)
    modulation: str = ""  # modulated structure: how its pattern was calculated (superspace)
    # only in the process that built the entry; set to None before it crosses process boundaries
    structure: gemmi.SmallStructure | None = field(repr=False, default=None)
    block: gemmi.cif.Block | None = field(repr=False, default=None)
    modulated: superspace.Modulated | None = field(repr=False, default=None)  # its superspace model, CIF setting
    images: sites.Images | None = field(repr=False, default=None)
    scattering: sites.Images | None = field(repr=False, default=None)  # images without dummy atoms
    site_orders: list | None = field(repr=False, default=None)  # site-symmetry orders the CIF states, per site

    @property
    def reference_code(self):
        """HighScore numbers COD entries as 96-<COD id + 1>."""
        return int(f"96{self.cod_id + 1}")


@lru_cache(maxsize=None)
def _corrections():
    return json.loads(Path(__file__).with_name("corrections.json").read_text())


_ANISO = ("u11", "u22", "u33", "u12", "u13", "u23")
MIN_UEQ = 1e-4  # Å²; a smaller equivalent value of an invalid tensor is no displacement parameter
# Å²; a tensor is invalid if a mean-square displacement is more negative than this (rounded values such as a
# U33 of -0.0001 pass); the factor keeps values exactly at the margin on the valid side despite rounding
ADP_MARGIN = 1e-4 * (1 + 1e-6)
_FIX_KEYS = {"source", "note", "sites", "types", "add", "spacegroup", "operations", "scattering_only", "omit_types",
             "cell", "superspace"}
_SITE_KEYS = {"species", "group", "occupancy_factor", "aniso", "u_iso", "stand_in", "shell", "omit"}


def _check_fix(cod_id, fix):
    """ValueError unless the correction is well-formed: a misspelt key or element must not be ignored silently."""
    def elements(d, what):
        if not d or not all(gemmi.Element(el).atomic_number > 0 and el == gemmi.Element(el).name for el in d):
            raise ValueError(f"corrections.json: COD {cod_id}: {what} must map element symbols, got {d!r}")

    unknown = set(fix) - _FIX_KEYS
    if unknown or not (fix.get("source") and fix.get("note")):
        raise ValueError(f"corrections.json: COD {cod_id}: unknown keys {sorted(unknown)} or no source/note")
    for label, change in fix.get("sites", {}).items():
        if set(change) - _SITE_KEYS or not change or {"species", "group"} <= set(change) or \
                ("stand_in" in change and not ({"species", "group"} & set(change) and change["stand_in"] is True)) or \
                ("omit" in change and change != {"omit": True}):
            raise ValueError(f"corrections.json: COD {cod_id} {label}: keys {sorted(change)}")
        if set(change.get("aniso", {})) - set(_ANISO):
            raise ValueError(f"corrections.json: COD {cod_id} {label}: aniso keys {sorted(change['aniso'])}")
        if "shell" in change and (set(change) != {"shell"} or set(change["shell"]) != {"radius", "per_cell"}
                                  or not change["shell"]["radius"] > 0):
            raise ValueError(f"corrections.json: COD {cod_id} {label}: shell needs radius and per_cell only")
        if "shell" in change:
            elements(change["shell"]["per_cell"], f"{label} shell")
        if "occupancy_factor" in change and not change["occupancy_factor"] > 0:
            raise ValueError(f"corrections.json: COD {cod_id} {label}: occupancy_factor must be positive")
        if "u_iso" in change and not change["u_iso"] > 0:  # 0 would read as missing (_default_b)
            raise ValueError(f"corrections.json: COD {cod_id} {label}: u_iso must be positive")
        for key in ("species", "group"):
            if key in change:
                elements(change[key], f"{label} {key}")
    for symbol, group in fix.get("types", {}).items():
        elements(group, f"type {symbol}")
    omit = fix.get("omit_types", [])
    if not isinstance(omit, list) or not all(isinstance(t, str) and t for t in omit):
        raise ValueError(f"corrections.json: COD {cod_id}: omit_types must list type symbols")
    for extra in fix.get("add", []):
        if not {"label", "xyz", "species"} <= set(extra) <= {"label", "xyz", "species", "u_iso"} or \
                len(extra["xyz"]) != 3 or not extra.get("u_iso", 0) >= 0:
            raise ValueError(f"corrections.json: COD {cod_id}: add needs label, xyz, species (and u_iso >= 0): "
                             f"{extra!r}")
        elements(extra["species"], f"add {extra['label']}")
    if "superspace" in fix:
        _superspace_block(cod_id, fix["superspace"])
    cell = fix.get("cell")
    if cell is not None and not (len(cell) == 6 and all(v > 0 for v in cell) and all(v < 180 for v in cell[3:])):
        raise ValueError(f"corrections.json: COD {cod_id}: cell must be a, b, c, alpha, beta, gamma: {cell!r}")
    if "spacegroup" in fix and gemmi.find_spacegroup_by_name(fix["spacegroup"]) is None:
        raise ValueError(f"corrections.json: COD {cod_id}: unknown spacegroup {fix['spacegroup']!r}")
    if "operations" in fix:
        from .setting import _is_group
        if not _is_group([gemmi.Op(o) for o in fix["operations"]]):
            raise ValueError(f"corrections.json: COD {cod_id}: the operations do not form a group")


def _superspace_block(cod_id, text):
    """The msCIF items of a correction's superspace model as a CIF block; ValueError unless it parses and holds
    wave vectors or superspace operations."""
    try:
        blk = gemmi.cif.read_string("data_superspace\n" + text).sole_block()
    except (RuntimeError, ValueError) as e:
        raise ValueError(f"corrections.json: COD {cod_id}: superspace is not CIF: {e}") from e
    if not (blk.find_values("_cell_wave_vector_x") or blk.find_values("_space_group_symop_ssg_operation_algebraic")):
        raise ValueError(f"corrections.json: COD {cod_id}: superspace gives no wave vectors or superspace operations")
    return blk


def _merge_items(block, extra):
    """The items and loops of extra written into block, replacing items and loops with the same tags (in memory;
    the CIF file is not changed)."""
    for item in extra:
        if item.pair is not None:
            block.set_pair(*item.pair)
        elif item.loop is not None:
            tags = list(item.loop.tags)
            prefix = os.path.commonprefix(tags)
            prefix = prefix[:prefix.rindex("_") + 1]
            loop = block.init_loop(prefix, [t[len(prefix):] for t in tags])
            for r in range(item.loop.length()):
                loop.add_row([item.loop[r, c] for c in range(item.loop.width())])


def _apply_corrections(st, cod_id):
    """Corrections from the publication for this CIF (corrections.json), applied to the structure as read; the
    CIF file is not changed. Returns the atom_site row of each resulting site (None for added atoms), its weight
    in the composition (formula, density, I/Ic), the correction's description (empty if none), the shell
    scatterers by atom_site row and the indices of the sites whose element the correction gives.

    Each entry needs "source" and "note"; _check_fix rejects unknown keys. "sites", per atom_site label:
    "species" (elements with absolute occupancies at the site), "group" (the atoms of a group scatterer, each with
    the site's occupancy), "occupancy_factor", "aniso" (tensor components u11...u23), "u_iso" (the published
    isotropic U in A^2, for a CIF that lacks it), "stand_in" (true: the
    species or group only stand in for the scattering of something else, such as guest molecules refined as
    single atoms, and do not count in formula, density, I/Ic and element lists), "shell" (the site is the centre
    of a spherical shell scatterer, sites.Shell: radius in Å and per_cell, the atoms spread over its shells), "omit"
    (true: the row is no atom, e.g. a ring centroid without a type, and is left out).
    "types", per type symbol: a
    group as above. "omit_types": type symbols of sites that are no atoms (e.g. interatomic scatterers for bonding
    density), left out. "add": atoms the CIF lacks, each with label, xyz and species. "spacegroup": the setting the
    coordinates are really in (e.g. "F d -3 m:1"); "operations": the symmetry operations, for a setting only the
    article lists. "superspace" (applied in _from_cif, before this): msCIF items with the article's superspace model
    (operations, wave vectors, Fourier waves) for a CIF that gives only the basic structure. Either replaces every symmetry statement of the CIF. "scattering_only": the occupancy factors
    and added atoms model scattering (e.g. an intergrowth of shifted layers) but are not atoms of the compound:
    weight 1/factor and 0. A label or type the CIF does not have is a KeyError: the file is wrong."""
    fix = _corrections().get(str(cod_id))
    if not fix:
        return list(range(len(st.sites))), [1.0] * len(st.sites), "", {}, set()
    _check_fix(cod_id, fix)
    by_site, types = fix.get("sites", {}), fix.get("types", {})
    missing = set(by_site) - {s.label for s in st.sites}
    if missing:
        raise KeyError(f"corrections.json: COD {cod_id} has no atom_site {sorted(missing)}")
    missing = (set(types) | set(fix.get("omit_types", []))) - {s.type_symbol for s in st.sites}
    if missing:
        raise KeyError(f"corrections.json: COD {cod_id} has no atom type {sorted(missing)}")
    if "spacegroup" in fix or "operations" in fix:  # replaces every symmetry statement of the CIF
        st.spacegroup_hm, st.spacegroup_hall = fix.get("spacegroup", ""), ""
        st.symops = fix.get("operations", [])
        st.determine_and_set_spacegroup("S" if st.symops else "1")  # H-M symbol with its setting, e.g. :1
    new, origin, weights, shells, assigned = [], [], [], {}, set()
    only_scattering = fix.get("scattering_only", False)
    for row, site in enumerate(st.sites):
        if site.type_symbol in fix.get("omit_types", []) or by_site.get(site.label, {}).get("omit"):
            # not atoms (e.g. bonding-density scatterers, a ring centroid)
            continue
        site = site.clone()
        change = by_site.get(site.label, {})
        weight = 0.0 if change.get("stand_in") else 1.0
        if "shell" in change:  # its centre, kept as a placeholder O site through the setting change
            shells[row] = change["shell"]
            site.element, site.type_symbol, site.occ, weight = gemmi.Element("O"), "O", 1.0, 0.0
        if "occupancy_factor" in change:
            site.occ *= change["occupancy_factor"]
            if only_scattering:
                weight = 1 / change["occupancy_factor"]
        if "aniso" in change:
            u = {k: getattr(site.aniso, k) for k in _ANISO} | change["aniso"]
            site.aniso = gemmi.SMat33d(*(u[k] for k in _ANISO))
        if "u_iso" in change:
            site.u_iso = change["u_iso"]
        if "species" in change:
            parts = [(el, occ, f"{site.label}_{el}") for el, occ in change["species"].items()]
        elif "group" in change or site.type_symbol in types:
            group = change.get("group") or types[site.type_symbol]
            parts = [(el, site.occ, site.label if k == 0 and j == 0 else f"{site.label}_{el}{k + 1}")
                     for j, (el, n) in enumerate(group.items()) for k in range(n)]
        else:
            parts = None
        for el, occ, label in parts or []:
            atom = site.clone()
            atom.element, atom.type_symbol, atom.occ, atom.label = gemmi.Element(el), el, occ, label
            assigned.add(len(new))
            new.append(atom)
            origin.append(row)
            weights.append(weight)
        if parts is None:
            new.append(site)
            origin.append(row)
            weights.append(weight)
    for extra in fix.get("add", []):
        for el, occ in extra["species"].items():
            atom = gemmi.SmallStructure.Site()
            atom.label, atom.type_symbol = f"{extra['label']}_{el}", el
            atom.element, atom.occ = gemmi.Element(el), occ
            atom.fract = gemmi.Fractional(*extra["xyz"])
            # the published U_iso, else B = 0.5 A^2 as for any atom without a displacement parameter
            atom.u_iso = extra.get("u_iso", sites.DEFAULT_B / (8 * math.pi**2))
            assigned.add(len(new))
            new.append(atom)
            origin.append(None)
            weights.append(0.0 if only_scattering else 1.0)
    st.sites = new
    return origin, weights, f"{fix['note']} ({fix['source']})", shells, assigned


def _read_aniso(block, st):
    """Anisotropic displacement parameters from the atom_site_aniso loop, as U_ij, B_ij or beta_ij.

    gemmi reads only U_ij and gives a tensor only to the first site of a label. Labels repeat in two ways:
    disordered parts (as many aniso rows as sites with the label, paired in order, or as many as sites the CIF
    marks Uani/Bani) and mixed sites (one row shared by atoms at the same position). U_ij = B_ij / (8 pi^2); U_ij = beta_ij / (2 pi^2 a*_i a*_j).
    Rows that are not numbers leave the atom isotropic.
    """
    for kind in ("U", "B", "beta"):
        table = block.find("_atom_site_aniso_", ["label"] + [f"{kind}_{ij}" for ij in ("11", "22", "33", "12", "13", "23")])
        if table.width() and len(table):
            break
    else:
        return
    rc = st.cell.reciprocal()
    star = [rc.a, rc.b, rc.c]
    scale = ([1.0] * 6 if kind == "U" else [1 / (8 * math.pi**2)] * 6 if kind == "B" else
             [1 / (2 * math.pi**2 * star[i] * star[j]) for i, j in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))])
    rows = {}
    for row in table:
        rows.setdefault(row.str(0), []).append([gemmi.cif.as_number(row[k]) * scale[k - 1] for k in range(1, 7)])
    adp_types = [gemmi.cif.as_string(v).strip().lower() for v in
                 block.find_values("_atom_site_adp_type") or block.find_values("_atom_site_thermal_displace_type")]
    sites, anisotropic = {}, {}
    for row, site in enumerate(st.sites):  # one site per atom_site row, in row order
        sites.setdefault(site.label, []).append(site)
        if row < len(adp_types) and adp_types[row] in ("uani", "bani"):
            anisotropic.setdefault(site.label, []).append(site)
    for label, group in sites.items():
        tensors = rows.get(label)
        if not tensors:
            continue
        for site in group:  # gemmi gave the first U row to the first site of the label
            site.aniso = gemmi.SMat33d(0, 0, 0, 0, 0, 0)
        if len(tensors) == len(group):
            pairs = zip(group, tensors)
        elif len(tensors) == len(anisotropic.get(label, [])):  # e.g. one disorder part Uani, the other Uiso
            pairs = zip(anisotropic[label], tensors)
        else:  # one tensor for the atoms of a mixed site
            first = group[0].fract
            pairs = ((s, tensors[0]) for s in group
                     if max(abs(a - b) for a, b in zip((s.fract.x, s.fract.y, s.fract.z), (first.x, first.y, first.z))) < 1e-4)
        for site, u in pairs:
            site.aniso = gemmi.SMat33d(*u) if all(map(math.isfinite, u)) else gemmi.SMat33d(0, 0, 0, 0, 0, 0)


_CHARGE = re.compile(r"\d*[+-]+\d*$")
_GROUP = re.compile(r"(?:[A-Z][a-z]?\d*)+")


def _element_of_type(type_symbol, formula_elements):
    """The element of a type symbol that is not one, when certain: it starts with an element symbol that the
    CIF's formula contains and the rest is a lowercase suffix or an atom label (Ow, Hy, Im1, HC13A, OCl(1)). None
    for generic or pseudo-atoms (T, M, SASH), mixed sites (ON, OW/Cl1) and groups (OH, CH, OH2)."""
    m = re.fullmatch(r"([A-Z])([a-z]?)(.*)", _CHARGE.sub("", type_symbol), re.S)
    if not m:
        return None
    element, rest = m.group(1) + m.group(2), m.group(3)
    if gemmi.Element(element).atomic_number == 0:
        element, rest = m.group(1), m.group(2) + rest
    if element not in formula_elements or "/" in rest:
        return None
    if re.search(r"[A-Z]", rest) and not re.search(r"\d", rest):
        return None  # capitals without a number: a code (SASH, SpHS, OiI), not an atom label such as HC13A
    if _GROUP.fullmatch(rest) and not re.search(r"\d[A-Za-z]", rest) and \
            all(gemmi.Element(e).atomic_number > 0 for e in re.findall(r"[A-Z][a-z]?", rest)):
        return None  # element symbols with counts: a group such as OH or OH2, or a mixed site such as ON
    return element


def _formula_elements(text):
    """Elements of a chemical formula; some old CIFs write it in lower case ("si12 zr4")."""
    if not re.search(r"[A-Z]", text):
        text = re.sub(r"[a-z]{1,2}", lambda m: m.group(0).capitalize(), text)
    return {e for e in re.findall(r"[A-Z][a-z]?", text) if gemmi.Element(e).atomic_number > 0}


def _capitals_misread(site, formula, has_types):
    """For an element read from two capitals (type or label, e.g. OS1, PB, CO3): None if the reading stands (the
    formula contains it, or a label like Co1 confirms it), else the element meant, or EntryError.

    Without a type column the capitals are a label: element + site letter in the AMCSD style (PB = P on site B,
    OS1 = O on site S1), so the first letter's element is taken if the formula has it, unless the rest is another
    element of the formula (CN). In a type column they are a code or group and are read as other types that are
    not elements (CO3 = carbonate: not certain)."""
    token = site.type_symbol if has_types else site.label
    if not re.match(r"[A-Z]{2}", token or "") or not formula or site.element.name in formula:
        return None
    if site.label.startswith(site.element.name) and site.label[:2] != site.label[:2].upper():
        return None  # label Co1, Mg: the element is meant; the formula leaves it out
    if not has_types:
        element, rest = token[0], token[1:]
        if rest in formula and element in formula:
            return _no_element(token)  # CN: carbon and nitrogen, not C on a site N
        return element if element in formula else _no_element(token)
    return _element_of_type(token, formula) or _no_element(token)


def _element_of_label(label, formula):
    """Element of an atom label that gemmi cannot read, in a CIF without types (HC12 = H on C12): the leading
    element if the formula contains it and the rest is not another element of the formula (CN). None otherwise."""
    m = re.match(r"([A-Z])([a-z]?)", label)
    if not m:
        return None
    for element, rest in ((m.group(1) + m.group(2), label[len(m.group(0)):]), (m.group(1), label[1:])):
        if gemmi.Element(element).atomic_number > 0 and element in formula:
            if rest in formula and not re.search(r"\d", rest):
                return None  # CN: carbon and nitrogen, not C with a site name
            return element
    return None


def _no_element(token):
    raise EntryError(f"atom type {token!r} is not an element")


def _described_elements(block):
    """Elements that the CIF's atom_type loop names for its type symbols: an _atom_type_description that is exactly
    an element symbol, with or without a charge (type TL described as Tl)."""
    out = {}
    for row in block.find("_atom_type_", ["symbol", "description"]):
        desc = _CHARGE.sub("", row.str(1).strip())
        if desc and gemmi.Element(desc).atomic_number > 0 and gemmi.Element(desc).name == desc:
            out[row.str(0)] = desc
    return out


def _resolve_types(block, st, origin, assigned=frozenset()):
    """Elements for atoms whose type symbol gemmi cannot read, or EntryError when the element is not certain:
    guessing it (a water O written WO read as tungsten) or leaving the atom out gives a wrong pattern. Elements
    read from two capitals that the formula does not contain are checked by _capitals_misread.

    Atoms whose type is '?' or '.', or whose type (or label, without a type column) names no element, are left out
    only when they are dummies (flagged dum or named Dummy/DUM, e.g. ring centroids, bond midpoints); '?' and '.'
    also at zero occupancy; otherwise the CIF is skipped. (Without a type column, the type is the
    label and is read as above.) Where the element would otherwise not be certain, the CIF's own atom_type
    description decides when it names one (_described_elements). origin: the atom_site row of each site;
    assigned: sites whose element a correction gives, which are not read again."""
    formula = _formula_elements(_text(block, "_chemical_formula_sum"))
    described = _described_elements(block)
    has_types = bool(block.find_values("_atom_site_type_symbol"))
    flags = [gemmi.cif.as_string(v).strip().lower() for v in block.find_values("_atom_site_calc_flag")]
    for i, site in enumerate(st.sites):
        if i in assigned or not all(map(math.isfinite, (site.fract.x, site.fract.y, site.fract.z))):
            continue
        if site.element != gemmi.Element("X"):
            try:
                meant = _capitals_misread(site, formula, has_types)
            except EntryError:
                if site.type_symbol not in described:
                    raise
                meant = described[site.type_symbol]
            if meant:
                site.element = gemmi.Element(meant)
            elif formula and site.element.name not in formula and not any(
                    re.match(rf"{site.element.name}(?![a-z])", t or "") for t in (site.type_symbol, site.label)):
                # read from a longer name (a ring centroid 'Cen' as cerium): an element the CIF does not name
                raise EntryError(f"atom {site.label!r} read as {site.element.name}, which the CIF's formula does "
                                 "not contain")
            continue
        row = origin[i]
        dummy = (row is not None and row < len(flags) and flags[row] == "dum") or site.label.lower().startswith("dum")
        if not site.type_symbol:
            if not dummy and not site.occ == 0:
                raise EntryError(f"atom {site.label!r} has no element")
            continue
        element = (_element_of_type(site.type_symbol, formula) if has_types
                   else _element_of_label(site.label, formula)) or described.get(site.type_symbol)
        if element is None:
            if dummy:  # a dummy without an element (a bond midpoint 'M1' in a CIF without type column) is left out
                continue
            raise EntryError(f"atom type {site.type_symbol!r} is not an element")
        site.element = gemmi.Element(element)


def _default_b(st, block, rows, ueq=None):
    """Isotropic displacement of sites without a tensor: U_iso_or_equiv, else B_iso_or_equiv (gemmi reads B_iso only
    without a U_iso column, so one given where U_iso is '?' or '.' is read here), else the equivalent value of a
    replaced invalid tensor (ueq, by site index, if at least MIN_UEQ), else B = sites.DEFAULT_B, HighScore's
    default, which the official COD databases also use for their patterns."""
    ueq = ueq or {}
    u_col, b_col = (block.find_values(t) for t in ("_atom_site_U_iso_or_equiv", "_atom_site_B_iso_or_equiv"))

    def value(column, row):
        return gemmi.cif.as_number(column[row]) if column and row is not None and row < len(column) else math.nan

    for i, (site, row) in enumerate(zip(st.sites, rows)):
        if row is None:  # an atom added by a correction: its value comes from the correction
            continue
        u, b = value(u_col, row), value(b_col, row)
        if site.aniso.nonzero() or site.u_iso != 0 or math.isfinite(u):
            continue
        if math.isfinite(b):
            site.u_iso = max(b, 0.0) / (8 * math.pi**2)
        elif ueq.get(i, 0.0) >= MIN_UEQ:
            site.u_iso = ueq[i]
        else:
            site.u_iso = sites.DEFAULT_B / (8 * math.pi**2)


def _clean_sites(st, spacegroup, shift=(0.0, 0.0, 0.0)):
    """Structure with only usable sites in the given space group, the atom_site row of each kept site, the
    number of displacement parameters that had to be replaced, and the equivalent value (trace / 3) of each
    replaced tensor, by site index, for _default_b.

    Drops sites without element or coordinates (e.g. read from stray atom_site loops) and moves the others by
    -shift. Reads missing occupancies as 1 and missing displacement parameters as 0 (see _default_b).
    Anisotropic tensors that are not physically valid (a negative mean-square displacement) are replaced by the
    isotropic value (see _default_b), and a negative isotropic value by 0. Raises EntryError for elements without X-ray scattering factors and when
    no site is left.
    """
    out = gemmi.SmallStructure()
    out.name = st.name
    out.cell = st.cell
    out.spacegroup_hm = spacegroup.xhm()
    out.determine_and_set_spacegroup("2")
    kept, ueq = [], {}
    invalid = 0
    for i, site in enumerate(st.sites):
        if site.element == gemmi.Element("X") or not all(map(math.isfinite, (site.fract.x, site.fract.y, site.fract.z))):
            continue
        if site.element.it92 is None:
            raise EntryError(f"no X-ray scattering factors for element {site.element.name}")
        n = site.clone()
        n.fract = gemmi.Fractional(n.fract.x - shift[0], n.fract.y - shift[1], n.fract.z - shift[2])
        if not math.isfinite(n.occ):
            n.occ = 1.0
        if not math.isfinite(n.u_iso):
            n.u_iso = 0.0
        u = n.aniso
        if not all(map(math.isfinite, (u.u11, u.u22, u.u33, u.u12, u.u13, u.u23))):
            n.aniso = gemmi.SMat33d(0, 0, 0, 0, 0, 0)
        elif any((u.u11, u.u22, u.u33, u.u12, u.u13, u.u23)) and \
                min(np.linalg.eigvalsh(ucart := sites.ucif_to_cart(n.aniso, st.cell))) < -ADP_MARGIN:
            # (gemmi's nonzero() tests only the trace, which misses tensors with an empty diagonal)
            n.aniso = gemmi.SMat33d(0, 0, 0, 0, 0, 0)  # falls back to the isotropic value (_default_b)
            ueq[len(kept)] = float(np.trace(ucart)) / 3
            invalid += 1
        if n.u_iso < 0:
            n.u_iso = 0.0
            invalid += 1
        out.add_site(n)
        kept.append(i)
    if not kept:
        raise EntryError("no atom coordinates")
    return out, kept, invalid, ueq


def _number(block, tag):
    """Numeric CIF item without its uncertainty; None when absent or not a number."""
    v = block.find_value(tag)
    if v is None:
        return None
    try:
        return float(gemmi.cif.as_string(v).split("(")[0])
    except ValueError:
        return None


def _site_orders(block, rows, n_ops):
    """Site-symmetry orders the CIF allows for each given atom_site row (None if it states nothing).

    From _atom_site_site_symmetry_order, else _atom_site_symmetry_multiplicity, which most programs fill with
    the multiplicity (n_ops / order) but SHELXL with the site-symmetry order, so both readings are allowed.
    """
    def column(tag):
        values = block.find_values(tag)
        return [gemmi.cif.as_string(v) for v in values] if values else []

    def positive_int(values, row):
        v = values[row].strip() if row is not None and row < len(values) else ""
        return int(v) if v.isdigit() and int(v) > 0 else None

    order, mult = column("_atom_site_site_symmetry_order"), column("_atom_site_symmetry_multiplicity")
    out = []
    for row in rows:
        k, m = positive_int(order, row), positive_int(mult, row)
        if k:
            out.append({k})
        elif m:
            out.append({m} | ({n_ops // m} if n_ops % m == 0 else set()))
        else:
            out.append(None)
    return out


def _text(block, tag):
    v = block.find_value(tag)
    if v is None:
        return ""
    v = gemmi.cif.as_string(v).strip()
    return "" if v in ("?", ".") else v


def _corundum_k():
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(_CORUNDUM_CIF).sole_block())
    return _k_factor(st, sites.unit_cell_images(st))


def _cell_mass(st, images, weights=None):
    """Unit-cell mass; weights scale each site's share (corrections that only model scattering count less)."""
    w = weights if weights is not None else [1.0] * len(st.sites)
    return (sum(st.sites[i].occ * w[i] * st.sites[i].element.weight for i in images.site)
            + sum(n * gemmi.Element(el).weight for sh in images.shells for el, n in sh.per_cell.items()))


def _k_factor(st, images, imax=None, weights=None):
    if imax is None:
        _, imax = pattern.calculate(st, 90.0, images=images)
    density = _cell_mass(st, images, weights) * 1.66054 / st.cell.volume
    return imax / st.cell.volume**2 / density


def _has_ch_bond(st, images, max_dist=1.25):
    carbons = [st.cell.orthogonalize(gemmi.Fractional(*x))
               for i, x in zip(images.site, images.xyz) if st.sites[i].element.name == "C"]
    if not carbons:
        return False
    for k in set(images.site):
        h = st.sites[k]
        if h.element.name not in ("H", "D"):
            continue
        ph = st.cell.orthogonalize(h.fract)
        if any(st.cell.find_nearest_image(ph, pc, gemmi.Asu.Any).dist() < max_dist for pc in carbons):
            return True
    return False


_K_CORUNDUM = None


def _symbol_group(block):
    """Space group of the CIF's Hermann-Mauguin symbol alone, or None."""
    probe = gemmi.SmallStructure()
    probe.spacegroup_hm = setting.normalise_symbol(
        _text(block, "_space_group_name_H-M_alt") or _text(block, "_symmetry_space_group_name_H-M"))
    probe.determine_and_set_spacegroup("12")
    return probe.spacegroup


def _composition_matches(e):
    """Whether the entry's unit-cell content agrees with the CIF's own statements (see _content_matches)."""
    have = {el: float(t[len(el):]) for el, t in zip(e.elements, e.formula.split())}
    return _content_matches(e.block, have, e.density)


def _content_matches(block, have, density):
    """Whether a unit-cell content (element -> atoms per cell) agrees with the CIF's own statements:
    _chemical_formula_sum x Z (_formula_matches); without formula or Z, the density calculated by the authors
    (_density_matches). None if the CIF states neither."""
    by_formula = _formula_matches(block, have)
    return by_formula if by_formula is not None else _density_matches(block, density)


def _formula_matches(block, have):
    """_chemical_formula_sum x Z against a unit-cell content, for the elements other than H and D (often not
    located), within 3%; None without formula or Z."""
    formula, z = _text(block, "_chemical_formula_sum"), _number(block, "_cell_formula_units_Z")
    if not (formula and z):
        return None
    want = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula):
        want[el] = want.get(el, 0.0) + (float(n) if n not in ("", ".") else 1.0) * z
    return all(abs(have.get(el, 0.0) - want.get(el, 0.0)) <= 0.03 * max(have.get(el, 0.0), want.get(el, 0.0))
               for el in set(want) | set(have) if el not in ("H", "D"))


def _density_matches(block, density):
    """The density calculated by the authors (_exptl_crystal_density_diffrn) within 3%; None if not given."""
    measured = _number(block, "_exptl_crystal_density_diffrn")
    return abs(density / measured - 1) < 0.03 if measured else None


def _crenel_occupancy_inside(block, model, st, atoms, weights):
    """Whether the CIF gives its crenel atoms' occupancies within their intervals rather than as averages. JANA
    writes averages (validated against the structure factors JANA lists, e.g. COD 2310780), but some CIFs
    (e.g. 2104133-2104141) give the values within the intervals. The CIF's formula x Z decides if exactly one
    reading agrees with it, else its density likewise; otherwise the average."""
    widths = {a.site: a.crenel[1] for a in model.atoms if a.crenel is not None and not a.occ_inside}
    if not widths:
        return False
    confined = {a.site: a.crenel[1] for a in model.atoms if a.occ_inside}  # sawtooth atoms: average occ x w

    def content(scale):
        counts, mass = {}, 0.0
        for i in atoms.site:
            s = st.sites[i]
            n = s.occ * weights[i] * scale.get(i, 1.0)
            counts[s.element.name] = counts.get(s.element.name, 0.0) + n
            mass += n * s.element.weight
        return counts, mass * 1.66054 / st.cell.volume

    (avg, avg_density), (inside, inside_density) = content(confined), content(confined | widths)
    for check in (lambda c, d: _formula_matches(block, c), lambda c, d: _density_matches(block, d)):
        a, b = check(avg, avg_density), check(inside, inside_density)
        if a is True and b is not True:
            return False
        if b is True and a is not True:
            return True
    return False


def from_cif(path, cod_id=None):
    """The reference pattern of a CIF.

    The listed symmetry operations (else the Hall symbol) define the space group. When they name another group
    than the Hermann-Mauguin symbol and do not give the CIF's density or formula (e.g. only half of the
    operations listed), the symbol is used if it does; if neither does, the CIF is skipped. Likewise, images of
    an atom on a special position that the CIF states as general are one atom, unless that contradicts the CIF's
    density or formula and counting them separately, as stated, does not."""
    e = _from_cif(path, cod_id)
    if e.symbol_conflict and _composition_matches(e) is False:
        by_symbol = _from_cif(path, cod_id, symbol_only=True)
        if not _composition_matches(by_symbol):
            raise EntryError(f"the listed symmetry operations ({e.sg_hm}) and the space-group symbol "
                             f"({by_symbol.sg_hm}) disagree, and neither gives the CIF's density or formula")
        by_symbol.symmetry_note = (f"symmetry from the space-group symbol {by_symbol.sg_hm}: the listed operations "
                                   f"({e.sg_hm}) do not give the CIF's density or formula")
        e = by_symbol
    if e.stated_general and _composition_matches(e) is False:
        separate = _from_cif(path, cod_id, symbol_only=bool(e.symmetry_note), separate_stated=True)
        if _composition_matches(separate):
            separate.symmetry_note = e.symmetry_note
            separate.merge_note = ("images of atoms on special positions that the CIF states as general (lower site "
                                   "symmetry) counted separately: merged, they do not give the CIF's density or formula")
            return separate
    return e


def _from_cif(path, cod_id=None, symbol_only=False, separate_stated=False):
    global _K_CORUNDUM
    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except RuntimeError as e:  # not CIF syntax, or not exactly one data block
        raise EntryError(f"unreadable CIF: {e}") from e
    cell_values = [_number(block, f"_cell_{t}") for t in ("length_a", "length_b", "length_c",
                                                          "angle_alpha", "angle_beta", "angle_gamma")]
    if not all(v is not None and v > 0 for v in cell_values) or not all(v < 180 for v in cell_values[3:]):
        raise EntryError("no valid unit cell")
    if cod_id is None:
        # the data block name when it is an id: COD 1010542.cif (data_1010542) carries another entry's
        # _cod_database_code (9016579)
        code = block.name if block.name.isdigit() else _text(block, "_cod_database_code") or block.name
        if not code.isdigit():
            raise EntryError(f"no COD id (data block {block.name!r})")
        cod_id = int(code)
    if not 1_000_000 <= cod_id <= 9_999_998:  # reference code 96 + (id + 1) must stay 9 digits
        raise EntryError(f"COD id {cod_id} outside the 7-digit range")
    fix = _corrections().get(str(cod_id), {})
    if "cell" in fix:  # the article's cell where the CIF's contradicts it
        _check_fix(cod_id, fix)
        cell_values = [float(v) for v in fix["cell"]]
    if "superspace" in fix:  # the article's superspace model, for a CIF that gives only the basic structure
        _check_fix(cod_id, fix)
        _merge_items(block, _superspace_block(cod_id, fix["superspace"]))
    try:
        st = gemmi.make_small_structure_from_block(block)
    except RuntimeError as e:
        raise EntryError(f"unreadable structure: {e}") from e
    c = st.cell
    if max(abs(x - y) for x, y in zip((c.a, c.b, c.c, c.alpha, c.beta, c.gamma), cell_values)) > 1e-6:
        st.cell = gemmi.UnitCell(*cell_values)  # gemmi misses values written on the line after their tag
    if not st.cell.volume > 0:
        raise EntryError("no valid unit cell")
    _read_aniso(block, st)
    origin, weights, correction, shell_specs, assigned = _apply_corrections(st, cod_id)
    _resolve_types(block, st, origin, assigned)
    shift = (0.0, 0.0, 0.0)
    if not setting.listed_operations(st, block) and not {"spacegroup", "operations"} & set(
            _corrections().get(str(cod_id), {})):
        # a modulated structure may list only superspace operations: their 3D parts are its basic space group
        st.symops = superspace.operations_3d(block)
        if st.symops:
            st.determine_and_set_spacegroup("S")
    stated = bool(st.symops or setting.listed_operations(st, block))  # listed operations or a Hall symbol
    if stated and not st.symops and st.spacegroup is not None and not symbol_only:
        # gemmi falls back to the H-M symbol when it does not know the Hall symbol's setting: the Hall symbol decides
        den = gemmi.Op.DEN
        key = lambda ops: {(tuple(map(tuple, o.rot)), tuple(t % den for t in o.tran)) for o in ops}  # noqa: E731
        if key(setting.listed_operations(st, block)) != key(st.spacegroup.operations()):
            st.determine_and_set_spacegroup("H")
    if symbol_only:
        st.symops, st.spacegroup_hall = [], ""
        st.spacegroup_hm = setting.normalise_symbol(st.spacegroup_hm)
        st.determine_and_set_spacegroup("12")
        if st.spacegroup is None:
            raise EntryError("unknown space group")
    elif st.spacegroup is None and not stated:
        # only a Hermann-Mauguin symbol, in a notation gemmi does not read directly. Listed operations or a Hall
        # symbol, when given, define the group even if the symbol is readable: they can have a shifted origin or
        # another setting than the symbol names
        st.spacegroup_hm = setting.normalise_symbol(st.spacegroup_hm)
        st.determine_and_set_spacegroup("12")
    spacegroup = st.spacegroup
    if spacegroup is None and st.symops:
        found = setting.identify_shifted(list(st.symops))  # operations of a known group with a shifted origin
        if found:
            spacegroup, shift = found
    ops = None
    if spacegroup is None:  # a setting gemmi does not know: standardise from the operations with spglib
        ops = setting.listed_operations(st, block)
        if not ops:
            raise EntryError("unknown space group")
        spacegroup = gemmi.SpaceGroup("P 1")  # until from_operations finds the group
    n_ops = len(ops) if ops else len(spacegroup.operations())  # in the CIF's own setting
    ops3d = ops or setting.listed_operations(st, block) or list(spacegroup.operations())  # in the CIF's setting
    dropped = len(st.sites)
    st, kept, invalid_adps, ueq = _clean_sites(st, spacegroup, shift)
    weights = [weights[i] for i in kept]
    dropped -= len(kept)
    site_rows = [origin[i] for i in kept]
    _default_b(st, block, site_rows, ueq)
    # a modulated structure whose superspace model cannot be evaluated (a composite, a function specific to JANA,
    # an incomplete CIF) is stored as its basic structure, as in COD24, and noted (unmodelled)
    unmodelled = ""
    try:
        modulated = superspace.read(block, st, site_rows, shift)
        if modulated:
            superspace.check_operations(modulated, ops3d)
    except superspace.ModulationError as e:
        modulated, unmodelled = None, str(e)
    basis = np.eye(3)  # x (final setting) = basis x (CIF setting) + t
    if ops:
        found = setting.from_operations(st, ops)
        if found is None:
            raise EntryError("unknown space group: the listed operations do not match the structure")
        st, basis = found
    to_ref = st.spacegroup.basisop.inverse()
    basis = np.array(to_ref.rot, dtype=float) / gemmi.Op.DEN @ basis
    st, setting_changed = setting.to_reference(st)
    st, reduced, cell_changed = setting.reduce_cell(st)
    basis = reduced @ basis
    setting_changed |= ops is not None or cell_changed
    corrected_symmetry = {"spacegroup", "operations"} & set(_corrections().get(str(cod_id), {}))
    symbol = _symbol_group(block) if stated and not symbol_only and not corrected_symmetry else None
    symbol_conflict = symbol is not None and symbol.number != st.spacegroup.number
    setting.constrain_cell(st)
    site_orders = _site_orders(block, site_rows, n_ops)
    images = sites.unit_cell_images(st, site_orders=site_orders, separate_stated=separate_stated)
    # dummy atoms (e.g. ring centroids) stay in the structure record, flagged, but do not scatter or count
    flags = [gemmi.cif.as_string(v).strip().lower() for v in block.find_values("_atom_site_calc_flag")]
    dummy = np.array([r is not None and r < len(flags) and flags[r] == "dum" for r in site_rows], dtype=bool)
    shell_sites = {k: shell_specs[r] for k, r in enumerate(site_rows) if r in shell_specs}
    shells = tuple(
        sites.Shell(images.xyz[images.site == k], spec["radius"], st.sites[k].u_iso, spec["per_cell"])
        for k, spec in shell_sites.items())
    real = ~(dummy | np.isin(np.arange(len(site_rows)), list(shell_sites)))[images.site]
    atoms = sites.Images(images.site[real], images.xyz[real], images.ucart[real], shells=shells)
    if not len(atoms.site):
        raise EntryError("no atoms other than dummy atoms")
    modulation = ""
    inside = False
    if modulated and _crenel_occupancy_inside(block, modulated, st, atoms, weights):
        inside = True
        for a in modulated.atoms:
            a.occ_inside |= a.crenel is not None
    if modulated:
        try:
            lines, imax, order, left = _modulated_pattern(modulated, atoms, dummy, shell_sites, basis, st.cell)
            check = superspace.agreement(modulated, block)  # the refinement's own structure factors, if listed
            if check and not check["R"] <= superspace.AGREEMENT_LIMIT:
                raise superspace.ModulationError(f"the structure factors the CIF lists are not reproduced (R "
                                                 f"{check['R']:.3f})")
        except superspace.ModulationError as e:
            modulated, unmodelled = None, str(e)
        else:
            for a in modulated.atoms:  # the basic structure holds the average occupancy
                if a.occ_inside:
                    st.sites[a.site].occ = a.occ * a.crenel[1]
            kind = ("commensurately modulated structure: pattern of its supercell" if modulated.t0 is not None
                    else f"modulated structure: pattern with satellites up to order {order}")
            modulation = (kind + (f" (superspace group {modulated.ssg})" if modulated.ssg else "")
                          + (f"; reproduces the structure factors the CIF lists (R {check['R']:.4f})" if check else "")
                          + ("; crenel occupancies read as the values within the intervals (the CIF's formula or "
                             "density)" if inside else "")
                          + (f"; higher orders left out (their strongest reflection {left:.3g} of 1000)" if left else "")
                          + "; the structure record holds the basic structure")
    if not modulated:
        try:
            lines, imax = pattern.stick_pattern(st, images=atoms)
        except ValueError as e:
            raise EntryError(str(e)) from e
        if unmodelled:
            modulation = (f"modulated structure, modulation not included ({unmodelled}): pattern of the basic "
                          "structure, as in COD24")
        elif superspace.states_modulation(block):
            modulation = ("modulated structure whose CIF gives only the average structure: pattern of the basic "
                          "structure, satellites missing")
    if not lines:
        raise EntryError("no reflections")
    if _K_CORUNDUM is None:
        _K_CORUNDUM = _corundum_k()
    iic = _k_factor(st, atoms, imax, weights) / _K_CORUNDUM
    density = _cell_mass(st, atoms, weights) * 1.66054 / st.cell.volume
    if not (math.isfinite(iic) and math.isfinite(density) and density > 0):
        raise EntryError("non-finite density or I/Ic")

    counts = {}
    for i in atoms.site:
        el = st.sites[i].element.name
        counts[el] = counts.get(el, 0.0) + st.sites[i].occ * weights[i]
    for sh in shells:
        for el, n in sh.per_cell.items():
            counts[el] = counts.get(el, 0.0) + n
    counts = {el: n for el, n in counts.items() if n > 0}  # zero occupancy, or atoms that only model scattering
    order = []
    for el in [site.element.name for site in st.sites] + [el for sh in shells for el in sh.per_cell]:
        if el in counts and el not in order:
            order.append(el)
    formula = " ".join(f"{el}{counts[el]:.2f}" for el in order)

    mineral = _text(block, "_chemical_name_mineral")
    compound = _text(block, "_chemical_name_systematic") or mineral or _text(block, "_chemical_name_common")
    compound = " ".join(compound.split()) or str(cod_id)

    comment = []
    colour = _text(block, "_exptl_crystal_colour")
    if colour:
        comment.append(f"Crystal color: {colour}")
    description = _text(block, "_exptl_crystal_description")
    if description:
        comment.append(f"Crystal description: {description}")
    title = _text(block, "_publ_section_title").replace("\n", " ")
    if title:
        comment.append(f"Publication title: {title}")
    if correction:
        comment.append(f"Corrected from the publication: {correction}")
    if modulation:
        comment.append(f"Pattern: {modulation}" + ("; satellite lines carry hkl 0 0 0" if modulated else ""))
    comment.append(f"COD database code: {cod_id}")

    first, last = _text(block, "_journal_page_first"), _text(block, "_journal_page_last")
    pages = (f"{first} - {last}" if last else first) if first.isdigit() else "0"
    authors = ", ".join(gemmi.cif.as_string(a) for a in block.find_values("_publ_author_name"))
    literature = {
        "volume": _text(block, "_journal_volume") or "0",
        "pages": pages,
        "year": _text(block, "_journal_year"),
        "authors": authors,
        "journal": _text(block, "_journal_name_full"),
    }

    subfiles = [SUBFILE_ORGANIC] if _has_ch_bond(st, atoms) else [SUBFILE_INORGANIC]
    if mineral:
        subfiles.append(SUBFILE_MINERAL)
    if all(el not in _NON_METALS for el in counts):
        subfiles.append(SUBFILE_METALLIC)

    c = st.cell
    return Entry(
        cod_id=cod_id,
        crystal_system=CRYSTAL_SYSTEM_CODES[st.spacegroup.crystal_system_str()],
        cell=(c.a, c.b, c.c, c.alpha, c.beta, c.gamma),
        volume=c.volume,
        sg_number=st.spacegroup.number,
        sg_hm=st.spacegroup.hm,
        density=density,
        iic=iic,
        elements=order,
        formula=formula,
        compound_name=compound,
        mineral_name=mineral,
        comment="\r".join(comment),
        literature=literature,
        subfiles=subfiles,
        lines=lines,
        site_rows=site_rows,
        setting_changed=setting_changed,
        symbol_conflict=symbol_conflict,
        shell_sites=tuple(shell_sites),
        stated_general="stated_general" in images.reasons,
        separate_stated=separate_stated,
        dropped_sites=dropped,
        invalid_adps=invalid_adps,
        dummy_sites=int(dummy.sum()),
        correction=correction,
        modulation=modulation,
        structure=st,
        block=block,
        modulated=modulated,
        images=images,
        scattering=atoms,
        site_orders=site_orders,
    )


def _modulated_pattern(model, atoms, dummy, shell_sites, basis, cell):
    """Stick pattern of a modulated structure (superspace.stick_pattern) in the stored setting. Its atoms are the
    scattering sites of the basic structure, each counted once per distinct image as the builder merges them.
    The stored cell may be smaller than the CIF's basic cell (a centred cell stored primitive): the images are
    counted per basic cell and the intensities scaled to the stored cell."""
    ratio = abs(np.linalg.det(basis))  # basic cell volume / stored cell volume
    per_basic = np.bincount(atoms.site, minlength=len(dummy)) * ratio
    model.atoms = [a for a in model.atoms if not dummy[a.site] and a.site not in shell_sites]
    for a in model.atoms:
        n = per_basic[a.site]
        if not n > 0 or abs(n - round(n)) > 1e-6 or len(model.ops) % round(n):
            raise superspace.ModulationError(f"{n:g} images of a site under {len(model.ops)} superspace operations")
        a.stab = len(model.ops) // round(n)
    lines, imax, order, left = superspace.stick_pattern(model, basis, cell)
    return lines, imax / ratio**2, order, left


def element_masks(elements):
    """Bit Z (Z < 64) in e1, bit Z-64 in e2."""
    e1 = e2 = 0
    for el in elements:
        z = gemmi.Element(el).atomic_number
        if z < 64:
            e1 |= 1 << z
        else:
            e2 |= 1 << (z - 64)
    # stored as signed 64-bit integers (bit 63 of e1, europium, is the sign bit)
    return tuple(v - (1 << 64) if v >= 1 << 63 else v for v in (e1, e2))

