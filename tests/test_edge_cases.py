"""Regression tests for problems found in full COD runs and code review, and edge cases of the formats."""

import math
import re
import sqlite3
import subprocess
import sys

import gemmi
import numpy as np
import pytest

from hsrdb_tools import entry, hsrdb, pattern, phase, setting, sites, structure

HEAD = """data_{cod}
_cod_database_code {cod}
_journal_name_full 'Test Journal'
_journal_year 2020
{extra}
_cell_length_a {a}
_cell_length_b {b}
_cell_length_c {c}
_cell_angle_alpha {al}
_cell_angle_beta {be}
_cell_angle_gamma {ga}
_symmetry_space_group_name_H-M '{sg}'
"""
ATOMS = """loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
"""


def cif(atoms, sg="F m -3 m", a=5.6402, b=None, c=None, al=90, be=90, ga=90, cod=1234567, extra="", tail=""):
    head = HEAD.format(cod=cod, extra=extra, a=a, b=b or a, c=c or a, al=al, be=be, ga=ga, sg=sg)
    return head + ATOMS + atoms + tail


NACL_ATOMS = "Na1 Na+ 0 0 0 1 0.012\nCl1 Cl- 0.5 0.5 0.5 1 0.011\n"


def build(tmp_path, text, name="1234567.cif"):
    p = tmp_path / name
    p.write_text(text)
    return entry.from_cif(p)


# --- bugs found by full runs --------------------------------------------------------------------------

def test_stray_site_without_element_is_dropped():
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(cif(NACL_ATOMS)).sole_block())
    bad = gemmi.SmallStructure.Site()
    bad.fract = gemmi.Fractional(math.nan, math.nan, math.nan)
    st.add_site(bad)
    clean, rows, _, _ = entry._clean_sites(st, st.spacegroup)
    assert [s.label for s in clean.sites] == ["Na1", "Cl1"] and rows == [0, 1]
    lines, _ = pattern.calculate(clean, 90)
    assert lines and all(math.isfinite(x.intensity) for x in lines)


def test_no_usable_sites():
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(cif(NACL_ATOMS)).sole_block())
    empty = gemmi.SmallStructure()
    empty.cell, empty.spacegroup_hm = st.cell, st.spacegroup.xhm()
    empty.determine_and_set_spacegroup("2")
    empty.add_site(gemmi.SmallStructure.Site())
    with pytest.raises(entry.EntryError, match="no atom coordinates"):
        entry._clean_sites(empty, empty.spacegroup)


def test_huge_page_number_and_text_volume(tmp_path):
    e = build(tmp_path, cif(NACL_ATOMS, extra="_journal_page_first 141007164223007\n_journal_volume 12a\n"
                                                "_journal_issue 99999999999"))
    p = structure.build(e, phase.HS4)
    assert p.trailer["page_first"] == 0 and p.trailer["volume"] == 0
    assert phase.encode_payload(p)  # fits the int32 fields


@pytest.mark.parametrize("text, expected", [
    ("0.1234(5)", (0.1234, 0.0005)), (".0210(10)", (0.021, 0.001)), ("0.012(3", (0.012, 0.003)),
    ("-1.5e-3(2)", (-0.0015, 0.0002)), ("1.", (1.0, 0.0)), ("-2", (-2.0, 0.0)),
    (".", None), ("?", None), ("abc", None),
])
def test_cif_numbers(text, expected):
    got = structure._num_str(text)
    assert got == expected if expected is None else got == pytest.approx(expected)


def test_missing_displacement_and_occupancy(tmp_path):
    e = build(tmp_path, cif("Na1 Na 0 0 0 . .\nCl1 Cl 0.5 0.5 0.5 1 0.011\n"))
    p = structure.build(e, phase.HS4)
    na = p.atoms[0]
    assert na.biso.value == pytest.approx(sites.DEFAULT_B) and na.sof.value == 1.0
    assert e.lines and all(math.isfinite(x.intensity) for x in e.lines)


def test_stated_zero_displacement_stays_zero_in_structure_record(tmp_path):
    # the pattern uses U = 0 as stated (and 0 for a negative U); the record must not write the 0.5 default
    e = build(tmp_path, cif("Na1 Na 0 0 0 1 0.0\nCl1 Cl 0.5 0.5 0.5 1 -0.01\n"))
    assert [s.u_iso for s in e.structure.sites] == [0.0, 0.0]
    p = structure.build(e, phase.HS4)
    lo = structure._templates()["trv"]["ATOM Biso"][3][3]  # HighScore's lower limit for Biso
    assert [a.biso.value for a in p.atoms] == [lo, lo]


def test_pseudo_element_is_rejected(tmp_path):
    with pytest.raises(entry.EntryError, match="scattering factors"):
        build(tmp_path, cif("Hg1 Hg 0 0 0 1 0.01\nCN1 CN 0.387 0.387 0.387 1 0.02\n"))


@pytest.mark.parametrize("cell_line, message", [
    ("_cell_length_a 5.6402\n", None),
    ("_cell_length_a ?\n", "no valid unit cell"),
])
def test_unit_cell_required(tmp_path, cell_line, message):
    text = cif(NACL_ATOMS).replace("_cell_length_a 5.6402\n", cell_line)
    if message:
        with pytest.raises(entry.EntryError, match=message):
            build(tmp_path, text)
    else:
        assert build(tmp_path, text).lines


def test_space_group_from_symmetry_operations(tmp_path):
    text = cif(NACL_ATOMS, sg="Q 9 9", tail="loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\n-x,-y,-z\n")
    assert build(tmp_path, text).sg_number == 2


def test_unknown_space_group(tmp_path):
    with pytest.raises(entry.EntryError, match="unknown space group"):
        build(tmp_path, cif(NACL_ATOMS, sg="Q 9 9"))


def test_cod_id_range(tmp_path):
    with pytest.raises(entry.EntryError, match="7-digit"):
        build(tmp_path, cif(NACL_ATOMS, cod=123))


# --- bugs found in review -----------------------------------------------------------------------------

ANISO_P422 = cif("C1 C 0.11 0.23 0.31 1 0.02\n", sg="P 4 2 2", a=5, c=7, tail="""loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
C1 0.01 0.03 0.02 0.005 0 0
""")


def test_anisotropic_tensors_rotate_with_symmetry():
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(ANISO_P422).sole_block())
    im = sites.unit_cell_images(st)
    orth = np.array(st.cell.orth.mat.tolist())
    for op, x, u in zip(st.spacegroup.operations(), im.xyz, im.ucart):
        r = np.array(op.rot) / gemmi.Op.DEN
        c = orth @ r @ np.linalg.inv(orth)
        u0 = sites.ucif_to_cart(st.sites[0].aniso, st.cell)
        assert np.allclose(u, c @ u0 @ c.T)


def test_anisotropic_intensities_match_direct_sum():
    """F(hkl) from the unit-cell images equals the sum over symmetry operations applied to the site."""
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(ANISO_P422).sole_block())
    im = sites.unit_cell_images(st)
    hkl = np.array([[1, 2, 3], [2, 1, 3], [3, 0, 1]])
    frac = np.array(st.cell.frac.mat.tolist())
    orth = np.array(st.cell.orth.mat.tolist())
    u0 = sites.ucif_to_cart(st.sites[0].aniso, st.cell)
    x0 = np.array(st.sites[0].fract.tolist()) if hasattr(st.sites[0].fract, "tolist") else \
        np.array([st.sites[0].fract.x, st.sites[0].fract.y, st.sites[0].fract.z])
    for h in hkl:
        s = h @ frac
        direct = 0
        for op in st.spacegroup.operations():
            r = np.array(op.rot) / gemmi.Op.DEN
            t = np.array(op.tran) / gemmi.Op.DEN
            c = orth @ r @ np.linalg.inv(orth)
            direct += np.exp(-2 * math.pi**2 * s @ c @ u0 @ c.T @ s) * np.exp(2j * math.pi * h @ (r @ x0 + t))
        ours = sum(np.exp(-2 * math.pi**2 * s @ u @ s) * np.exp(2j * math.pi * h @ x)
                   for x, u in zip(im.xyz, im.ucart))
        assert abs(ours - direct) < 1e-9


def test_special_position_given_to_three_decimals():
    text = cif("C1 C 0.333 0.667 0.25 1 0.02\n", sg="P 63/m m c", a=3, c=5, ga=120)
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(text).sole_block())
    assert list(sites.unit_cell_images(st).multiplicities(1)) == [2]
    assert structure.site_info(st, sites.unit_cell_images(st))[0] == (2, "2c")


def test_occupancy_and_limits_per_version(tmp_path):
    text = cif("C1 C 0.1 0.2 0.3 5 0.02\n", sg="P -1", a=5, b=6, c=7, al=80, be=85, ga=95, tail="""loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
C1 0.02 0.02 0.02 -0.005 0 0
""")
    e = build(tmp_path, text)
    for fmt in (phase.HS3, phase.HS4):
        a = structure.build(e, fmt).atoms[0]
        assert a.sof.value == 1.0  # HighScore's limit, as in the official databases
        lo, hi = a.aniso[3].values[3:5]
        assert lo <= a.aniso[3].value <= hi  # B12 within this version's limits


def test_duplicate_labels_use_their_own_rows(tmp_path):
    text = cif("X1 Na+ 0 0 0 1 0.012\nX1 Cl- 0.5 0.5 0.5 0.5 0.030\n")
    e = build(tmp_path, text)
    p = structure.build(e, phase.HS4)
    assert [a.element for a in p.atoms] == ["Na", "Cl"]
    assert p.atoms[1].sof.value == 0.5
    assert p.atoms[1].biso.value == pytest.approx(0.030 * structure.EIGHT_PI2)


def test_hs3_text_is_transliterated(tmp_path):
    e = build(tmp_path, cif(NACL_ATOMS, extra="_chemical_name_mineral 'Hématite α'"))
    payload = phase.encode_payload(structure.build(e, phase.HS3))
    p = phase.decode_payload(payload)
    assert p.trailer["mineral"] == "Hematite ?"
    p4 = phase.decode_payload(phase.encode_payload(structure.build(e, phase.HS4)))
    assert p4.trailer["mineral"] == "Hématite α"


def test_encoder_rejects_misaligned_fields(tmp_path):
    e = build(tmp_path, cif(NACL_ATOMS))
    p = structure.build(e, phase.HS4)
    p.atoms[0].tail = bytes(25)
    with pytest.raises(ValueError):
        phase.encode_payload(p)


def test_invalid_lines_are_rejected():
    with pytest.raises(ValueError):
        hsrdb.encode_lines([pattern.Line(2.0, math.nan, (1, 0, 0))])


# --- format details -----------------------------------------------------------------------------------

def test_element_masks():
    assert entry.element_masks(["Eu"]) == (-(2**63), 0)  # Z = 63: sign bit of the signed 64-bit column
    assert entry.element_masks(["Gd"]) == (0, 1)
    assert entry.element_masks(["H", "U"]) == (2, 1 << 28)
    db = sqlite3.connect(":memory:")
    db.execute("create table t(e1 integer)")
    db.execute("insert into t values (?)", (entry.element_masks(["Eu"])[0],))
    assert db.execute("select e1 from t").fetchone()[0] == -(2**63)


def test_reference_code(tmp_path):
    e = build(tmp_path, cif(NACL_ATOMS, cod=1000032))
    assert e.reference_code == 961000033
    db = hsrdb.create(tmp_path / "x.hsrdb")
    hsrdb.add(db, e, entry.element_masks(e.elements), hsrdb.encode_lines(e.lines), hsrdb.phase_blob(e), "hs4")
    hsrdb.update_prefixes(db)
    assert db.execute("select PhaseName from PatternPhase").fetchone()[0] == "PDF:96-100-0033"
    assert db.execute("select Value from Properties where Key='ReferenceCodePrefixes'").fetchone()[0] == "96-100\n"


def test_rhombohedral_axes_become_hexagonal(tmp_path):
    # corundum on rhombohedral axes
    a_r, alpha = 5.1284, 55.28
    text = cif("Al1 Al 0.35216 0.35216 0.35216 1 0.003\nO1 O 0.55624 0.94376 0.25 1 0.0035\n",
               sg="R -3 c :R", a=a_r, al=alpha, be=alpha, ga=alpha)
    e = build(tmp_path, text)
    assert e.setting_changed and e.sg_hm == "R -3 c"
    ca = math.cos(math.radians(alpha))
    a_h, c_h = 2 * a_r * math.sin(math.radians(alpha) / 2), a_r * math.sqrt(3 * (1 + 2 * ca))
    assert e.cell == pytest.approx((a_h, a_h, c_h, 90, 90, 120), abs=1e-6)
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(text).sole_block())
    before, _ = pattern.calculate(st, 90)
    after, _ = pattern.calculate(e.structure, 90)
    assert [round(x.d, 6) for x in before] == [round(x.d, 6) for x in after]
    assert [x.intensity for x in before] == pytest.approx([x.intensity for x in after], abs=1e-6)


def test_origin_choice_2(tmp_path):
    text = cif("Si1 Si 0 0 0 1 0.005\n", sg="F d -3 m :1", a=5.431)
    e = build(tmp_path, text)
    assert e.setting_changed and e.sg_hm == "F d -3 m"
    x = e.structure.sites[0].fract
    assert abs(((x.x - 0.125) * 8) % 2) < 1e-9 or abs(((x.x - 0.125) * 8) % 2 - 2) < 1e-9  # on the 8a/8b set
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(text).sole_block())
    before, _ = pattern.calculate(st, 90)
    after, _ = pattern.calculate(e.structure, 90)
    assert [(round(a.d, 6), round(a.intensity, 6)) for a in before] == \
        [(round(b.d, 6), round(b.intensity, 6)) for b in after]
    assert structure.build(e, phase.HS4).atoms[0].multiplicity == 8


def test_line_cap_and_order():
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(
        cif("C1 C 0.1 0.2 0.3 1 0.02\n", sg="P 1", a=20, b=21, c=22)).sole_block())
    lines, _ = pattern.stick_pattern(st)
    assert len(lines) == pattern.MAX_LINES
    assert all(x.d > y.d for x, y in zip(lines, lines[1:]))
    assert max(x.intensity for x in lines) == pytest.approx(1000)
    full, _ = pattern.calculate(st, 90)
    assert [x.d for x in lines] == [x.d for x in full[:pattern.MAX_LINES]]


def test_few_lines_extend_to_140_degrees():
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(
        cif("Cu1 Cu 0 0 0 1 0.005\n", a=3.615)).sole_block())
    lines, _ = pattern.stick_pattern(st)
    assert len(lines) > 5 and min(x.d for x in lines) < pattern.d_for_two_theta(90)


def test_coincident_reflections_merge():
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(
        cif("Na1 Na 0 0 0 1 0.01\n", sg="P m -3 m", a=4)).sole_block())
    lines, _ = pattern.calculate(st, 90)
    at = [x for x in lines if abs(x.d - 4 / 3) < 1e-6]  # (300) and (221) coincide
    assert len(at) == 1


def test_inverted_structure_gives_the_same_pattern():
    """Friedel mates share a powder line, so a structure and its inverted image (x -> -x) give one pattern even
    when anomalous scattering makes |F(h)| and |F(-h)| differ."""
    atoms = "U1 U 0.11 0.23 0.37 1 0.01\nO1 O 0.31 0.42 0.05 1 0.02\nN1 N 0.47 0.13 0.29 1 0.02\n"
    inverted = "".join(f"{lab} {el} {-float(x):.2f} {-float(y):.2f} {-float(z):.2f} {o} {u}\n"
                       for lab, el, x, y, z, o, u in (r.split() for r in atoms.splitlines()))
    patterns = []
    for rows in (atoms, inverted):
        st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(
            cif(rows, sg="P 21 21 21", a=5.1, b=6.3, c=7.2)).sole_block())
        lines, _ = pattern.calculate(st, 90)
        patterns.append([(round(x.d, 6), x.intensity) for x in lines])
    assert [d for d, _ in patterns[0]] == [d for d, _ in patterns[1]]
    assert np.allclose([i for _, i in patterns[0]], [i for _, i in patterns[1]], rtol=1e-9, atol=1e-9)


def test_many_atoms_round_trip(tmp_path):
    rows = "".join(f"C{i} C {(i * 0.0137) % 1:.4f} {(i * 0.0291) % 1:.4f} {(i * 0.0453) % 1:.4f} 1 0.02\n"
                   for i in range(300))
    e = build(tmp_path, cif(rows, sg="P 1", a=30, b=31, c=32))
    for fmt in (phase.HS3, phase.HS4):
        payload = phase.encode_payload(structure.build(e, fmt))
        assert phase.encode_payload(phase.decode_payload(payload)) == payload
        assert len(phase.decode_payload(payload).atoms) == 300


def test_no_reflections(tmp_path):
    with pytest.raises(entry.EntryError, match="no reflections"):
        build(tmp_path, cif("Na1 Na 0 0 0 0 0.01\n"))


# --- command line -------------------------------------------------------------------------------------

def run_cli(*args):
    return subprocess.run([sys.executable, "-m", "hsrdb_tools", "add", *map(str, args)],
                          capture_output=True, text=True)


def test_cli_skips_bad_cifs_and_reports(tmp_path):
    src = tmp_path / "cifs"
    src.mkdir()
    (src / "1234567.cif").write_text(cif(NACL_ATOMS))
    (src / "1234568.cif").write_text(cif(NACL_ATOMS, cod=1234568).replace("_cell_length_a 5.6402\n", ""))
    out = tmp_path / "db.hsrdb"
    r = run_cli(src, "--hs4", out, "--hs3", tmp_path / "db3.hsrdb", "-j", "1")
    assert r.returncode == 0, r.stderr
    assert "skipped 1" in r.stdout
    report = (tmp_path / "db.skipped.tsv").read_text()
    assert "1234568.cif" in report and "no valid unit cell" in report
    for path in (out, tmp_path / "db3.hsrdb"):
        db = sqlite3.connect(path)
        assert db.execute("select count(*) from general_indexed").fetchone()[0] == 1
        assert db.execute("select Value from Properties where Key='NoOfPatterns'").fetchone()[0] == "1"
    assert not list(tmp_path.glob("*.partial"))


def test_cli_append_duplicates_and_guards(tmp_path):
    src = tmp_path / "cifs"
    src.mkdir()
    (src / "a.cif").write_text(cif(NACL_ATOMS))
    (src / "b.cif").write_text(cif(NACL_ATOMS))  # same COD id
    out = tmp_path / "db.hsrdb"
    r = run_cli(src, "--hs4", out, "-j", "1")
    assert r.returncode == 0 and "already present 1" in r.stdout
    assert "b.cif\tnote: not added: same COD id as" in (tmp_path / "db.skipped.tsv").read_text()
    r = run_cli(src, "--hs4", out, "--append", "-j", "1")
    assert r.returncode == 0 and "added {}" in r.stdout
    assert run_cli(src, "--hs4", out, "-j", "1").returncode != 0  # exists
    r = run_cli(src, "--hs3", out, "--append", "-j", "1")
    assert r.returncode == 1 and "schema is for hs4, not hs3" in r.stderr  # wrong schema, clean exit
    assert run_cli(src, "--hs4", tmp_path / "c.hsrdb", "--hs3", tmp_path / "d.hsrdb", "--base", out).returncode != 0
    assert run_cli(src, "--hs4", tmp_path / "e.hsrdb", "--base", out, "--append").returncode != 0


def test_form_factors_match_international_tables():
    # International Tables Vol. C, Table 6.1.1.1, oxygen at sin(theta)/lambda = 0 ... 0.5; the IT92 fit is within 0.015
    stol = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    assert pattern._f0("O", stol**2) == pytest.approx([8.000, 7.250, 5.634, 4.094, 3.010, 2.338], abs=0.015)


def test_nacl_intensities_by_hand(tmp_path):
    """mult * |F|^2 * LP with F = 4 (f_Na T_Na +- f_Cl T_Cl), f = f0 + f' + i f'', T = exp(-8 pi^2 U s^2)."""
    e = build(tmp_path, cif(NACL_ATOMS))
    a, lam = 5.6402, pattern.CU_KA1
    expected = {}
    for hkl, mult in (((1, 1, 1), 8), ((2, 0, 0), 6), ((2, 2, 0), 12), ((3, 1, 1), 24), ((2, 2, 2), 8)):
        d = a / math.sqrt(sum(h * h for h in hkl))
        s2 = 1 / (4 * d * d)
        f = {el: gemmi.Element(el).it92.calculate_sf(s2) + complex(*pattern._anomalous(el, lam)) for el in ("Na", "Cl")}
        t = {"Na": math.exp(-8 * math.pi**2 * 0.012 * s2), "Cl": math.exp(-8 * math.pi**2 * 0.011 * s2)}
        sign = 1 if hkl[0] % 2 == 0 else -1
        F = 4 * (f["Na"] * t["Na"] + sign * f["Cl"] * t["Cl"])
        th = math.asin(lam / (2 * d))
        expected[round(d, 4)] = mult * abs(F) ** 2 * (1 + math.cos(2 * th) ** 2) / (math.sin(th) ** 2 * math.cos(th))
    top = max(expected.values())
    got = {round(x.d, 4): x.intensity for x in e.lines}
    for d, i in expected.items():
        assert got[d] == pytest.approx(1000 * i / top, rel=1e-6)  # gemmi evaluates the IT92 sum in float32


@pytest.mark.parametrize("type_symbol, label, formula, element", [
    ("CO3", "C1", "C Ca O3", None),   # carbonate group in the type column: not cobalt, not certain
    ("CO", "Co1", "Bi Fe O3", "Co"),  # the label confirms cobalt, a dopant the formula leaves out
])
def test_two_capitals_in_type_column(tmp_path, type_symbol, label, formula, element):
    text = cif(f"Ca1 Ca 0.1 0.2 0.3 1 0.01\n{label} {type_symbol} 0.3 0.1 0.2 1 0.01\n", sg="P 1", a=5, b=6, c=7,
               extra=f"_chemical_formula_sum '{formula}'")
    if element is None:
        with pytest.raises(entry.EntryError, match="not an element"):
            build(tmp_path, text)
    else:
        assert build(tmp_path, text).structure.sites[1].element.name == element


def test_zero_occupancy_element_is_not_in_the_compound(tmp_path):
    e = build(tmp_path, cif("Na1 Na 0 0 0 1 0.01\nK1 K 0 0 0 0 0.01\nCl1 Cl 0.5 0.5 0.5 1 0.01\n"))
    assert e.elements == ["Na", "Cl"] and e.formula == "Na4.00 Cl4.00"


def test_b_iso_only(tmp_path):
    text = cif(NACL_ATOMS).replace("_atom_site_U_iso_or_equiv", "_atom_site_B_iso_or_equiv")
    text = text.replace("Na1 Na+ 0 0 0 1 0.012", "Na1 Na+ 0 0 0 1 1.25(3)")
    e = build(tmp_path, text)
    na = structure.build(e, phase.HS4).atoms[0]
    assert na.biso.value == pytest.approx(1.25) and na.biso.values[1] == pytest.approx(0.03)


def test_invalid_anisotropic_tensor_falls_back_to_isotropic(tmp_path):
    # U12 far larger than U11 and U22: a negative mean-square displacement (a typo in some COD CIFs)
    text = cif("C1 C 0.1 0.2 0.3 1 0.02\nO1 O 0.3 0.1 0.2 1 0.03\n", sg="P -1", a=6.5, b=7.4, c=8.7,
               al=70.2, be=83.3, ga=66.3, tail="""loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
C1 0.02 0.02 0.02 0.9 0 0
""")
    e = build(tmp_path, text)
    assert e.invalid_adps == 1 and not e.structure.sites[0].aniso.nonzero()
    assert len(e.lines) > 10 and all(math.isfinite(x.intensity) for x in e.lines)
    c1 = structure.build(e, phase.HS4).atoms[0]
    assert c1.biso.value == pytest.approx(0.02 * structure.EIGHT_PI2)
    # without U_iso_or_equiv: the tensor's own equivalent value (trace / 3), not the 0.5 default
    no_iso = text.replace("_atom_site_U_iso_or_equiv\n", "").replace(" 1 0.02\n", " 1\n").replace(" 1 0.03\n", " 1\n")
    e = build(tmp_path, no_iso, name="no_iso.cif")
    cell = gemmi.UnitCell(6.5, 7.4, 8.7, 70.2, 83.3, 66.3)
    ueq = np.trace(sites.ucif_to_cart(gemmi.SMat33d(0.02, 0.02, 0.02, 0.9, 0, 0), cell)) / 3
    assert e.invalid_adps == 1 and e.structure.sites[0].u_iso == pytest.approx(ueq)
    assert e.structure.sites[1].u_iso == pytest.approx(sites.DEFAULT_B / structure.EIGHT_PI2)


def test_b_iso_where_u_iso_is_unknown(tmp_path):
    """U_iso '?' and B_iso 3.1 (COD 2005981): gemmi reads only the U column; B = 3.1 is the displacement."""
    rows = "Na1 Na 0 0 0 1 ? 3.1\nCl1 Cl 0.5 0.5 0.5 1 0.011 0.87\n"
    text = cif(rows).replace("_atom_site_U_iso_or_equiv\n", "_atom_site_U_iso_or_equiv\n_atom_site_B_iso_or_equiv\n")
    e = build(tmp_path, text)
    assert [s.u_iso for s in e.structure.sites] == pytest.approx([3.1 / structure.EIGHT_PI2, 0.011])


@pytest.mark.parametrize("u_iso, b_iso, expected", [
    ("0.03", "?", 0.03),                          # U_iso_or_equiv first
    ("?", "2.0", 2.0 / (8 * math.pi**2)),         # then B_iso_or_equiv, before the tensor's own value
    ("?", "?", None),                             # then trace / 3 (below)
])
def test_invalid_tensor_with_an_empty_diagonal(tmp_path, u_iso, b_iso, expected):
    """U11 = U22 = U33 = 0 with off-diagonal terms (COD 1100978, 2104388): invalid, although gemmi's nonzero()
    looks only at the trace."""
    rows = f"C1 C 0.1 0.2 0.3 1 {u_iso} {b_iso}\nO1 O 0.3 0.1 0.2 1 0.03 ?\n"
    tail = ANISO_LOOP + "C1 0.004 0 0 0.010 0 0\n"
    text = cif(rows, sg="P -1", a=6, b=7, c=8, al=80, be=85, ga=95, tail=tail)
    text = text.replace("_atom_site_U_iso_or_equiv\n", "_atom_site_U_iso_or_equiv\n_atom_site_B_iso_or_equiv\n")
    e = build(tmp_path, text)
    c1 = e.structure.sites[0]
    assert e.invalid_adps == 1 and not c1.aniso.nonzero()
    if expected is None:
        expected = np.trace(sites.ucif_to_cart(gemmi.SMat33d(0.004, 0, 0, 0.010, 0, 0), gemmi.UnitCell(6, 7, 8, 80, 85, 95))) / 3
    assert c1.u_iso == pytest.approx(expected)


def test_tensor_at_the_rounding_margin_is_valid(tmp_path):
    """U33 = -0.0001 (COD 9016426): within the rounding margin, the tensor is kept."""
    text = cif("Ce1 Ce 0 0.75 0.125 1 0.0025\n", sg="I 41/a m d :2", a=7.3308, c=6.4356,
               tail=ANISO_LOOP + "Ce1 0.0037 0.0037 -0.0001 0 0 0\n")
    e = build(tmp_path, text)
    assert e.invalid_adps == 0 and e.structure.sites[0].aniso.u33 == pytest.approx(-0.0001)


def test_tensor_with_only_off_diagonal_terms_is_invalid(tmp_path):
    rows = "C1 C 0.1 0.2 0.3 1 ?\nO1 O 0.3 0.1 0.2 1 0.03\n"
    text = cif(rows, sg="P 1", a=6, b=7, c=8, tail=ANISO_LOOP + "C1 0 0 0 0.010 0 0\n")
    e = build(tmp_path, text)
    assert e.invalid_adps == 1 and e.structure.sites[0].u_iso == pytest.approx(sites.DEFAULT_B / structure.EIGHT_PI2)


def test_cell_value_on_next_line(tmp_path):
    text = cif(NACL_ATOMS).replace("_cell_length_a 5.6402\n", "_cell_length_a\n5.6402\n")
    e = build(tmp_path, text)
    assert e.cell[0] == pytest.approx(5.6402)


def test_space_group_with_shifted_origin(tmp_path):
    # P 21/c operations with the origin moved by (1/4, 0, 0) and no usable symbol
    ops = ["x,y,z", "1/2-x,1/2+y,1/2-z", "1/2-x,-y,-z", "x,1/2-y,1/2+z"]
    tail = "loop_\n_symmetry_equiv_pos_as_xyz\n" + "".join(f"'{o}'\n" for o in ops)
    text = cif("C1 C 0.35 0.2 0.3 1 0.02\n", sg="P 21/c (x+1/4,y,z)", a=5, b=6, c=7, be=100, tail=tail)
    e = build(tmp_path, text)
    assert e.sg_number == 14
    unshifted = cif("C1 C 0.1 0.2 0.3 1 0.02\n", sg="P 1 21/c 1", a=5, b=6, c=7, be=100)
    ref = build(tmp_path, unshifted, name="ref.cif")
    assert [(round(x.d, 5), round(x.intensity, 3)) for x in e.lines] == \
        [(round(x.d, 5), round(x.intensity, 3)) for x in ref.lines]


# --- settings gemmi does not tabulate: operations standardised with spglib ------------------------------

def ops_loop(ops):
    return "loop_\n_symmetry_equiv_pos_as_xyz\n" + "".join(f"'{o.triplet()}'\n" for o in ops)


def listed_pattern(text):
    """Pattern of a CIF expanded with its own listed operations in its own cell (no standardisation)."""
    block = gemmi.cif.read_string(text).sole_block()
    st = gemmi.make_small_structure_from_block(block)
    ops = setting.listed_operations(st, block)
    st, _, _, _ = entry._clean_sites(st, gemmi.SpaceGroup("P 1"))
    images = sites.unit_cell_images(st, ops)
    p1 = gemmi.SmallStructure()
    p1.cell = st.cell
    p1.spacegroup_hm = "P 1"
    p1.determine_and_set_spacegroup("2")
    for i, x in zip(images.site, images.xyz):
        s = st.sites[i].clone()
        s.fract = gemmi.Fractional(*x)
        p1.add_site(s)
    return pattern.calculate(p1, 90)[0]


def same_lines(a, b):
    """Same lines: equal intensities, positions within the merge width (the reference computed in P 1 splits
    symmetry-equivalent reflections, so its strongest reflection of a merged line can differ)."""
    assert len(a) == len(b)
    two_theta = lambda lines: np.degrees(2 * np.arcsin(pattern.CU_KA1 / (2 * np.array([x.d for x in lines]))))  # noqa: E731
    assert np.abs(two_theta(a) - two_theta(b)).max() <= pattern.MERGE_TWO_THETA
    assert np.allclose([x.intensity for x in a], [x.intensity for x in b], atol=1e-6)


TWO_ATOMS = "Fe1 Fe 0.113 0.271 0.349 1 0.011\nS1 S 0.362 0.083 0.187 1 0.013\n"


@pytest.mark.parametrize("symbol, hall, cell, number", [
    # B-centred cell with a d glide: P 21/c on (a+c, b, c)-type axes
    ("B 1 21/d 1", "-P 2ybc (x+1/2*z,y,1/2*z)", dict(a=5.74, b=5.65, c=5.74, be=111.9), 14),
    # F-centred tetragonal cell of I -4 2 d
    ("F -4 d 2", "I -4 2bw (1/2*x-1/2*y,1/2*x+1/2*y,z)", dict(a=10.48, c=6.9), 122),
    # C 1 2 1 with the 2-fold along c on an F-centred cell (COD: 'F 1 1 2', Hall 'C 2y (1/2*z,x-1/2*z,y)')
    ("F 1 1 2", "C 2y (1/2*z,x-1/2*z,y)", dict(a=16.28, b=6.43, c=11.78, ga=133.6), 5),
    # supercell description of P 21/c, with its extra translations listed as operations
    ("P 1 21/c 1 (2*c,2*a+c,b)", "-P 2ybc (-1/4*x+1/2*z,1/2*x,y)", dict(a=7.83, b=5.22, c=5.54), 14),
    # trailing junk and a rhombohedral cell in the reverse setting
    ("R -3 m HR", '-R 3 2" (-x,-y,z)', dict(a=3.97, c=28.59, ga=120), 166),
])
def test_operations_in_unknown_settings(tmp_path, symbol, hall, cell, number):
    ops = list(gemmi.symops_from_hall(hall))
    assert gemmi.find_spacegroup_by_ops(gemmi.symops_from_hall(hall)) is None  # not in gemmi's table
    text = cif(TWO_ATOMS, sg=symbol, **cell, tail=ops_loop(ops))
    e = build(tmp_path, text)
    assert e.sg_number == number and e.setting_changed
    same_lines(pattern.calculate(e.structure, 90)[0], listed_pattern(text))


def test_origin_shift_by_a_third(tmp_path):
    # COD 1011159 'P 32 2 1 S': quartz with the origin moved by c/3 (identify_shifted tries eighths only)
    ops = ["x,y,z", "y,x,2/3-z", "-y,x-y,2/3+z", "-x,y-x,1/3-z", "y-x,-x,1/3+z", "x-y,-y,-z"]
    tail = "loop_\n_symmetry_equiv_pos_as_xyz\n" + "".join(f"'{o}'\n" for o in ops)
    text = cif("Si1 Si 0.469 0 0 1 0.01\nO1 O 0.403 0.253 0.122 1 0.015\n", sg="P 32 2 1 S", a=4.91, c=5.4, ga=120,
               tail=tail)
    e = build(tmp_path, text)
    assert e.sg_number == 154 and e.formula == "Si3.00 O6.00"
    same_lines(pattern.calculate(e.structure, 90)[0], listed_pattern(text))


def test_hall_symbol_with_change_of_basis_only(tmp_path):
    text = cif(TWO_ATOMS, sg="?", a=5.74, b=5.65, c=5.74, be=111.9).replace(
        "_symmetry_space_group_name_H-M '?'", "_symmetry_space_group_name_Hall '-P 2ybc (x+1/2*z,y,1/2*z)'")
    e = build(tmp_path, text)
    assert e.sg_number == 14
    same_lines(pattern.calculate(e.structure, 90)[0], listed_pattern(text))


@pytest.mark.parametrize("symbol, reference", [("C2:b1", "C 1 2 1"),
                                               ("F d d d {origin @ -1 @ d d d}", "F d d d :2")])
def test_symbol_notations(tmp_path, symbol, reference):
    kw = dict(a=14.16, b=14.44, c=10.06) if symbol.startswith("F") else dict(a=5.41, b=9.0, c=10.25, be=100.3)
    e = build(tmp_path, cif(TWO_ATOMS, sg=symbol, **kw))
    ref = build(tmp_path, cif(TWO_ATOMS, sg=reference, **kw), name="ref.cif")
    assert e.sg_number == ref.sg_number
    same_lines(e.lines, ref.lines)


def test_operations_that_do_not_fit_the_cell_are_rejected(tmp_path):
    # 2-fold screw along c with d glide, but the cell's oblique angle is beta (COD 5000046)
    ops = list(gemmi.symops_from_hall("-P 2cb (x+1/2*y,1/2*y,z)"))
    text = cif(TWO_ATOMS, sg="C 1 1 21/d", a=5.74, b=5.65, c=5.76, be=110.6, tail=ops_loop(ops))
    with pytest.raises(entry.EntryError, match="do not match"):
        build(tmp_path, text)


# --- atoms near special positions ---------------------------------------------------------------------

def near_mirror(tmp_path, rows, tail="", columns=""):
    """P 1 2/m 1 (mirror at y = 0, b = 10 A): a site at y = 0.015 has its mirror image 0.3 A away."""
    text = cif(rows, sg="P 1 2/m 1", a=6, b=10, c=7, be=100, tail=tail)
    if columns:
        text = text.replace("_atom_site_U_iso_or_equiv\n", "_atom_site_U_iso_or_equiv\n" + columns)
    e = build(tmp_path, text)
    return e.images.multiplicities(len(e.structure.sites)), e


def test_split_atom_near_a_mirror_keeps_both_positions(tmp_path):
    mult, e = near_mirror(tmp_path, "C1 C 0.1 0.015 0.2 0.5 0.02\n")
    assert list(mult) == [4] and e.formula == "C2.00"  # two half-occupied positions per mirror pair


def test_full_atom_near_a_mirror_sits_on_it(tmp_path):
    mult, e = near_mirror(tmp_path, "C1 C 0.1 0.015 0.2 1 0.02\n")  # two full atoms 0.3 A apart are impossible
    assert list(mult) == [2] and e.formula == "C2.00"
    assert np.allclose(e.images.xyz[:, 1] % 0.5, 0)  # moved onto the mirror


def test_mixed_site_near_a_mirror_counts_both_elements(tmp_path):
    mult, _ = near_mirror(tmp_path, "Fe1 Fe 0.1 0.015 0.2 0.6 0.02\nMg1 Mg 0.1 0.015 0.2 0.4 0.02\n")
    assert list(mult) == [2, 2]  # 0.6 + 0.4 on each of two images exceeds 1: one site on the mirror


@pytest.mark.parametrize("column, value, expected", [
    ("_atom_site_site_symmetry_order", 2, 2),  # stated on the mirror
    ("_atom_site_site_symmetry_order", 1, 4),  # stated general: a split atom
    ("_atom_site_symmetry_multiplicity", 2, 2),  # multiplicity of the mirror position
    ("_atom_site_symmetry_multiplicity", 1, 4),  # SHELXL: site-symmetry order 1
])
def test_stated_site_symmetry_decides(tmp_path, column, value, expected):
    mult, _ = near_mirror(tmp_path, f"C1 C 0.1 0.015 0.2 0.5 0.02 {value}\n", columns=column + "\n")
    assert list(mult) == [expected]


@pytest.mark.parametrize("occ, positions, formula", [
    (0.25, 8, "Tl2.00 O2.00"),  # four quarter-occupied positions around each 2a site: a split atom (COD 1521425)
    (0.5, 2, "Tl1.00 O2.00"),   # four half-occupied positions 0.36 A apart are impossible: one atom on 2a
])
def test_split_atom_around_a_special_position(tmp_path, occ, positions, formula):
    """Positions x,x,0 on a square of side 0.36 A and diagonal 0.51 A: merged or split as a whole, never in part.
    The occupancy rule counts distinct positions (4), not the operations that map onto them (16)."""
    rows = f"Tl1 Tl 0.048 0.048 0 {occ} 0.01\nO1 O 0.5 0.5 0 1 0.01\n"
    e = build(tmp_path, cif(rows, sg="I 4/m m m", a=3.7667, c=29.39))
    assert list(e.images.multiplicities(2)) == [positions, 2] and e.formula == formula
    if positions == 2:
        tl = e.images.xyz[e.images.site == 0]
        assert np.allclose(((tl + 0.25) % 0.5) - 0.25, 0)  # exactly on 0,0,0 and 1/2,1/2,1/2


def test_merged_images_merge_again_while_closer_than_the_tolerance(tmp_path):
    """COD 9013334: four images of a 0.58-occupied atom around a -4 point of P 42/n, 0.39 A apart in pairs whose
    centres are 0.35 A apart: one atom on the -4 point, not two over-occupied ones."""
    e = build(tmp_path, cif("C1 C 0.266 0.248 0.726 0.58 0.01\nCa1 Ca 0.1 0.03 0.74 1 0.01\n", sg="P 42/n :2",
                            a=12.1338, c=7.5755))
    assert list(e.images.multiplicities(2)) == [2, 8]
    c = e.images.xyz[e.images.site == 0]
    assert np.allclose(np.sort(c[:, 2]), [0.25, 0.75]) and np.allclose(c[:, :2] % 0.5, 0.25)


@pytest.mark.parametrize("carbon, expected, separate", [("0.25", "Pt8.00 C2.00", False), ("0.5", "Pt8.00 C4.00", True)])
def test_images_of_a_site_stated_general_on_a_special_position(tmp_path, carbon, expected, separate):
    """C1 exactly on a 2-fold axis of C2/c, occupancy 0.5, but stated general (SHELXL multiplicity 1, as in a
    negative PART; COD 4123499, 7126502). One atom per axis position, unless the CIF's formula x Z says the
    images count separately."""
    rows = "Pt1 Pt 0.1 0.2 0.3 1 0.01 1\nC1 C 0 0.3 0.25 0.5 0.02 1\n"
    text = cif(rows, sg="C 1 2/c 1", a=10, b=11, c=12, be=100,
               extra=f"_chemical_formula_sum 'C{carbon} Pt'\n_cell_formula_units_Z 8")
    text = text.replace("_atom_site_U_iso_or_equiv\n", "_atom_site_U_iso_or_equiv\n_atom_site_symmetry_multiplicity\n")
    e = build(tmp_path, text)
    assert e.stated_general and e.formula == expected
    assert e.separate_stated == separate and bool(e.merge_note) == separate


def test_rounded_special_position_is_one_atom(tmp_path):
    mult, e = near_mirror(tmp_path, "C1 C 0.1 0.0004 0.2 0.5 0.02\n")  # 0.008 A from its image
    assert list(mult) == [2] and e.formula == "C1.00"


def test_zero_intensity_reflection_does_not_join_lines(monkeypatch):
    st = gemmi.make_small_structure_from_block(gemmi.cif.read_string(
        cif("C1 C 0.1 0.2 0.3 1 0.02\n", sg="P 1", a=17.3, b=19.1, c=23.7, al=81, be=77, ga=69)).sole_block())
    hkl = np.asarray(gemmi.make_miller_array(st.cell, st.spacegroup, pattern.d_for_two_theta(90)))
    d = 1 / np.linalg.norm(hkl @ np.array(st.cell.frac.mat.tolist()), axis=1)
    order = np.argsort(np.degrees(2 * np.arcsin(pattern.CU_KA1 / (2 * d))))
    t = np.degrees(2 * np.arcsin(pattern.CU_KA1 / (2 * d[order])))
    step = pattern.MERGE_TWO_THETA
    i = next(i for i in range(len(t) - 2) if t[i + 1] - t[i] <= step and t[i + 2] - t[i + 1] <= step < t[i + 2] - t[i])
    outer, middle = {tuple(hkl[order[i]]), tuple(hkl[order[i + 2]])}, tuple(hkl[order[i + 1]])

    def fake(st, h, d, wavelength, images):
        return np.array([10.0 if tuple(x) in outer else 0.0 if tuple(x) == middle else 1.0 for x in h])

    monkeypatch.setattr(pattern, "_intensities", fake)
    lines, _ = pattern.calculate(st, 90)
    assert {round(d[order[i]], 9), round(d[order[i + 2]], 9)} <= {round(x.d, 9) for x in lines}


@pytest.mark.parametrize("kind", ["B", "beta"])
def test_anisotropic_b_and_beta(tmp_path, kind):
    """gemmi reads only aniso U_ij; B_ij and beta_ij must give the same displacement."""
    a, b, c = 5.2, 6.1, 7.3
    u = {"11": 0.012, "22": 0.018, "33": 0.015, "12": 0.002, "13": 0.0, "23": 0.0}
    star = {"1": 1 / a, "2": 1 / b, "3": 1 / c}
    conv = {"U": lambda ij, v: v, "B": lambda ij, v: 8 * math.pi**2 * v,
            "beta": lambda ij, v: 2 * math.pi**2 * star[ij[0]] * star[ij[1]] * v}

    def text(k):
        tail = "loop_\n_atom_site_aniso_label\n" + "".join(f"_atom_site_aniso_{k}_{ij}\n" for ij in u) + \
            "Ba1 " + " ".join(f"{conv[k](ij, v):.8f}" for ij, v in u.items()) + "\n"
        return cif("Ba1 Ba 0.1 0.2 0.3 1 0.015\nO1 O 0.3 0.4 0.1 1 0.02\n", sg="P 21 21 21", a=a, b=b, c=c, tail=tail)

    ref = build(tmp_path, text("U"), name="ref.cif")
    e = build(tmp_path, text(kind))
    assert e.structure.sites[0].aniso.u11 == pytest.approx(0.012) and e.structure.sites[0].aniso.u12 == pytest.approx(0.002)
    same_lines(e.lines, ref.lines)


def test_dummy_atoms_do_not_scatter(tmp_path):
    rows = "Fe1 Fe 0.113 0.271 0.349 1 0.011\nS1 S 0.362 0.083 0.187 1 0.013\n"
    columns = "_atom_site_calc_flag\n"
    plain = cif(rows.replace("\n", " d\n"), sg="P 1 21/c 1", a=5.7, b=5.6, c=5.7, be=112)
    with_dummy = cif(rows.replace("\n", " d\n") + "Cg1 C 0.25 0.25 0.25 1 0 dum\n", sg="P 1 21/c 1",
                     a=5.7, b=5.6, c=5.7, be=112)
    plain, with_dummy = (t.replace("_atom_site_U_iso_or_equiv\n", "_atom_site_U_iso_or_equiv\n" + columns)
                         for t in (plain, with_dummy))
    ref = build(tmp_path, plain, name="ref.cif")
    e = build(tmp_path, with_dummy)
    assert e.formula == ref.formula and e.elements == ["Fe", "S"] and e.density == pytest.approx(ref.density)
    assert e.dummy_sites == 1 and ref.dummy_sites == 0
    same_lines(e.lines, ref.lines)
    atoms = structure.build(e, phase.HS4).atoms
    assert [a.label for a in atoms] == ["Fe1", "S1", "Cg1"] and atoms[2].tail[20] == 2  # kept, flagged dummy


def test_missing_displacement_uses_highscore_default(tmp_path):
    missing = build(tmp_path, cif("Na1 Na 0 0 0 1 .\nCl1 Cl 0.5 0.5 0.5 1 .\n"))
    b05 = cif("Na1 Na 0 0 0 1 0.5\nCl1 Cl 0.5 0.5 0.5 1 0.5\n").replace("_atom_site_U_iso_or_equiv",
                                                                     "_atom_site_B_iso_or_equiv")
    same_lines(missing.lines, build(tmp_path, b05, name="ref.cif").lines)


@pytest.mark.parametrize("sg, given, stored", [
    # (the axes of P b c a also permuted cyclically to the shortest a, test_cell_choice_as_cod24)
    ("P b c a", dict(a=38.89, b=20.452, c=8.701, al=89, be=89), (8.701, 38.89, 20.452, 90, 90, 90)),
    ("P 1 21/c 1", dict(a=5.546, b=7.854, c=23.336, al=89.99, be=93.18, ga=89.98), (5.546, 7.854, 23.336, 90, 93.18, 90)),
    ("P 4/m n c", dict(a=14.274, b=14.214, c=15.104), (14.214, 14.214, 15.104, 90, 90, 90)),
    ("P -3", dict(a=2.99, c=4.72), (2.99, 2.99, 4.72, 90, 90, 120)),
    ("F m -3 m", dict(a=5.86, b=5.865, c=5.861), (5.861, 5.861, 5.861, 90, 90, 90)),
])
def test_cell_constrained_to_crystal_system(tmp_path, sg, given, stored):
    """As in the official COD databases: exact angles, a = b (tetragonal, hexagonal) or a = b = c (cubic)."""
    e = build(tmp_path, cif(NACL_ATOMS, sg=sg, **given))
    assert e.cell == pytest.approx(stored)


@pytest.mark.parametrize("atoms, metallic", [
    ("Ce1 Ce 0 0 0 1 0.01\nGe1 Ge 0.5 0.5 0.5 1 0.01\n", True),  # Ge and Sb count as metals, as in COD24
    ("Ni1 Ni 0 0 0 1 0.01\nSb1 Sb 0.5 0.5 0.5 1 0.01\n", True),
    ("La1 La 0 0 0 1 0.01\nD1 D 0.5 0.5 0.5 1 0.02\n", False),  # a deuteride is not a metal, like a hydride
    ("Ni1 Ni 0 0 0 1 0.01\nAs1 As 0.5 0.5 0.5 1 0.01\n", False),
])
def test_metallic_subfile(tmp_path, atoms, metallic):
    e = build(tmp_path, cif(atoms, sg="P m -3 m", a=3.5))
    assert (entry.SUBFILE_METALLIC in e.subfiles) == metallic


ANISO_LOOP = "loop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\n_atom_site_aniso_U_22\n_atom_site_aniso_U_33\n" \
             "_atom_site_aniso_U_12\n_atom_site_aniso_U_13\n_atom_site_aniso_U_23\n"


def test_disorder_parts_with_one_label_get_their_own_tensors(tmp_path):
    rows = "C1 C 0.10 0.20 0.30 0.8 0.02\nC1 C 0.12 0.22 0.31 0.2 0.05\nO1 O 0.3 0.1 0.2 1 0.02\n"
    tail = ANISO_LOOP + "C1 0.011 0.012 0.013 0 0 0\nC1 0.051 0.052 0.053 0 0 0\n"
    e = build(tmp_path, cif(rows, sg="P -1", a=6, b=7, c=8, al=80, be=85, ga=95, tail=tail))
    c1 = [s.aniso.u11 for s in e.structure.sites if s.label == "C1"]
    assert c1 == pytest.approx([0.011, 0.051])


@pytest.mark.parametrize("types, expected", [(("Uani", "Uiso"), [0.011, 0.0]), (("Uiso", "Uani"), [0.0, 0.011])])
def test_one_tensor_for_a_repeated_label_goes_to_the_uani_atom(tmp_path, types, expected):
    """Two disorder parts share a label, one refined anisotropically (COD 1515303): the one aniso row is its."""
    rows = f"C1 C 0.10 0.20 0.30 0.7 0.02 {types[0]}\nC1 C 0.14 0.25 0.33 0.3 0.05 {types[1]}\nO1 O 0.3 0.1 0.2 1 0.02 Uiso\n"
    text = cif(rows, sg="P -1", a=6, b=7, c=8, al=80, be=85, ga=95, tail=ANISO_LOOP + "C1 0.011 0.012 0.013 0 0 0\n")
    text = text.replace("_atom_site_U_iso_or_equiv\n", "_atom_site_U_iso_or_equiv\n_atom_site_adp_type\n")
    e = build(tmp_path, text)
    assert [s.aniso.u11 for s in e.structure.sites if s.label == "C1"] == pytest.approx(expected)


def test_mixed_site_shares_its_tensor(tmp_path):
    rows = "M1 Fe 0.1 0.2 0.3 0.6 0.02\nM1 Mn 0.1 0.2 0.3 0.4 0.02\nO1 O 0.3 0.1 0.2 1 0.02\n"
    tail = ANISO_LOOP + "M1 0.011 0.012 0.013 0.001 0 0\n"
    e = build(tmp_path, cif(rows, sg="P -1", a=6, b=7, c=8, al=80, be=85, ga=95, tail=tail))
    assert [s.aniso.u11 for s in e.structure.sites if s.label == "M1"] == pytest.approx([0.011, 0.011])


@pytest.mark.parametrize("type_symbol, formula, element", [
    ("Ow", "C14 H33 Cl3 N4 O2 Zn", "O"),      # water oxygen (COD 1557919)
    ("HC13A", "C28 H38 N2 O5 S", "H"),        # an atom label in the type column (COD 2006347)
    ("HN1B", "B2 C28 H38 N2 O5 S", "H"),      # label N1B, not nitrogen and boron
    ("OH2", "C8 H20 N O Si", None),           # a water group
    ("OCl(1)", "C14 H42 Cl2 O7 Pd S6", "O"),  # COD 2004740
    ("Im1", "Cu H2 I3 O2 Pb2", "I"),          # COD 1561471
    ("OH", "C8 H20 N O Si", None),            # a hydroxyl group scatterer, not an atom
    ("OH-", "Bi2 Ca2 H2.5 O12 Si3", None),
    ("ON", "C11 H21 Cd N7 Ni O", None),       # mixed O/N site
    ("OW/Cl1", "C26 H33 Cl2 N4 Ni O5", None),
    ("WO", "Si36 O76 N2 C8 H40", None),       # water read by COD24 as tungsten (COD 7201135)
    ("ST1", "F4 K3 Nb3 O11 Ti2", None),       # Nb/Ti site read by COD24 as sulfur (COD 2004310)
    ("SASH", "Al56 Na56 O384 Si136 S1", None),  # undefined pseudo-atom, even in a sulfur compound
    ("T", "Al0.024 Na0.015 O4 Si1.978", None),  # generic tetrahedral site
])
def test_atom_types_that_are_not_elements(type_symbol, formula, element):
    known = set(re.findall(r"[A-Z][a-z]?", formula))
    assert entry._element_of_type(type_symbol, known) == element


def test_uncertain_atom_type_skips_the_cif(tmp_path):
    text = cif("Si1 Si 0.1 0.2 0.3 1 0.01\nT1 T 0.3 0.1 0.2 1 0.01\n", sg="P 1", a=5, b=6, c=7,
               extra="_chemical_formula_sum 'O4 Si2'")
    with pytest.raises(entry.EntryError, match="not an element"):
        build(tmp_path, text)
    ok = build(tmp_path, text.replace("T1 T ", "O1 Ow "), name="ok.cif")
    assert ok.formula == "Si1.00 O1.00"


def test_atom_without_element_skips_the_cif_unless_dummy(tmp_path):
    # COD 9002229: interlayer sites "I1" with type '?' are Mg or O, not iodine
    rows = "Si1 Si 0.1 0.2 0.3 1 0.01\nI1 ? 0.3 0.1 0.2 1 0.02\n"
    with pytest.raises(entry.EntryError, match="has no element"):
        build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7))
    dummy = build(tmp_path, cif(rows.replace("I1 ?", "Dummy1 ?"), sg="P 1", a=5, b=6, c=7), name="d.cif")
    assert dummy.formula == "Si1.00" and dummy.dropped_sites == 1


def test_corrections_from_publications(tmp_path, monkeypatch):
    fixes = {"1234567": {"source": "A. Author, J. Test 1, 1 (2000)", "note": "test fixes",
                         "sites": {"I1": {"species": {"Mg": 0.25, "O": 0.25}}, "O1": {"aniso": {"u33": 0.03}}},
                         "types": {"CE": {"C": 1, "H": 2}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    rows = "Si1 Si 0.1 0.2 0.3 1 0.01\nI1 ? 0.3 0.1 0.2 1 0.02\nC1 CE 0.4 0.4 0.1 0.5 0.03\nO1 O 0.2 0.3 0.4 1 0.02\n"
    tail = ANISO_LOOP + "O1 0.011 0.012 0.013 0 0 0\n"
    text = cif(rows, sg="P 1", a=5, b=6, c=7, tail=tail)
    e = build(tmp_path, text)
    assert e.formula == "Si1.00 Mg0.25 O1.25 C0.50 H1.00"
    o1 = next(s for s in e.structure.sites if s.label == "O1")
    assert o1.aniso.u33 == pytest.approx(0.03) and o1.aniso.u11 == pytest.approx(0.011)
    assert "Corrected from the publication: test fixes (A. Author" in e.comment and e.correction
    assert (tmp_path / "1234567.cif").read_text() == text  # the CIF itself is not changed
    fixes["1234567"]["sites"]["Missing"] = {"occupancy_factor": 0.5}
    with pytest.raises(KeyError):
        build(tmp_path, text, name="again.cif")


def test_corrections_file_is_well_formed():
    for cod, fix in entry._corrections().items():
        if not cod.startswith("_"):
            entry._check_fix(cod, fix)


@pytest.mark.parametrize("change, message", [
    ({"scatering_only": True}, "unknown keys"),                               # misspelt top-level key
    ({"sites": {"O1": {"occupancy": 0.5}}}, "keys"),                          # misspelt site key
    ({"sites": {"O1": {"aniso": {"U11": 0.5}}}}, "aniso keys"),              # CIF-style component name
    ({"sites": {"O1": {"species": {"MG": 0.5}}}}, "element symbols"),        # element in capitals
    ({"spacegroup": "F m -3 mm"}, "unknown spacegroup"),
    ({"operations": ["x,y,z", "-x,-y,z+1/2", "x+1/2,y,z"]}, "do not form a group"),
    ({"omit_types": "IS4"}, "omit_types"),
    ({"cell": [5, 6, 7, 90, 90]}, "cell"),
    ({"sites": {"O1": {"omit": True, "occupancy_factor": 0.5}}}, "keys"),   # omit stands alone
    ({"add": [{"label": "X", "xyz": [0, 0, 0], "species": {"O": 1}, "u_iso": -0.01}]}, "u_iso"),
])
def test_malformed_correction_is_an_error(tmp_path, monkeypatch, change, message):
    fixes = {"1234567": {"source": "s", "note": "n", **change}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    with pytest.raises(ValueError, match=message):
        build(tmp_path, cif("O1 O 0.1 0.2 0.3 1 0.01\n", sg="P 1", a=5, b=6, c=7))


def test_correction_gives_the_whole_structure(tmp_path, monkeypatch):
    """A CIF without atoms and with a cell that contradicts its space group (COD 6000286): the article's cell and
    atoms, each with its published U_iso (also 0)."""
    text = cif("? ? ? ? ? ? ?\n", sg="P 3 c 1", a=10.1405, b=10.1496, c=21.66)
    fixes = {"1234567": {"source": "s", "note": "n", "cell": [10.14, 10.14, 21.66, 90, 90, 120], "spacegroup": "P 3 c 1",
                         "add": [{"label": "Ba1", "xyz": [0.0, 0.354, 0.016], "species": {"Ba": 1}, "u_iso": 0.009},
                                 {"label": "Cu10", "xyz": [0.927, 0, 0.015], "species": {"Cu": 0.212}, "u_iso": 0.0},
                                 {"label": "O1", "xyz": [0.326, 0.502, 0.005], "species": {"O": 1}}]}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, text)
    assert e.cell == pytest.approx((10.14, 10.14, 21.66, 90, 90, 120)) and e.sg_hm == "P 3 c 1"
    u = {s.label: s.u_iso for s in e.structure.sites}
    assert u == pytest.approx({"Ba1_Ba": 0.009, "Cu10_Cu": 0.0, "O1_O": 0.5 / (8 * math.pi**2)})


def test_correction_gives_a_site_its_published_u_iso(tmp_path, monkeypatch):
    """A CIF without displacement parameters whose article lists U_iso per atom (COD 2300057): the correction's
    value, kept through the species split; the other atoms get B = 0.5 as usual."""
    rows = "C1 C 0.1 0.2 0.3 1.17 ?\nO8 C 0.3 0.1 0.2 1 ?\nO10 O 0.2 0.3 0.1 1 ?\n"
    fixes = {"1234567": {"source": "s", "note": "n",
                         "sites": {"C1": {"species": {"C": 1, "H": 1}, "u_iso": 0.073}, "O8": {"species": {"O": 1}}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7))
    u = {s.label: s.u_iso for s in e.structure.sites}
    assert u == pytest.approx({"C1_C": 0.073, "C1_H": 0.073, "O8_O": 0.5 / (8 * math.pi**2),
                               "O10": 0.5 / (8 * math.pi**2)})
    assert e.formula == "C1.00 H1.00 O2.00"
    fixes["1234567"]["sites"]["C1"]["u_iso"] = 0.0  # would be replaced by the default B silently
    with pytest.raises(ValueError, match="u_iso"):
        build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7), name="bad.cif")


def test_dummy_without_element_in_cif_without_type_column(tmp_path):
    """Double-bond midpoints M1, M2 flagged dum in a CIF without type column (COD 2010123) are no atoms."""
    head = cif("", sg="P 1", a=5, b=6, c=7).split("loop_\n_atom_site_label")[0]
    loop = ("loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
            "_atom_site_U_iso_or_equiv\n_atom_site_calc_flag\n")
    rows = "Rh1 0.1 0.2 0.3 0.02 ?\nC1 0.3 0.1 0.2 0.03 ?\n"
    e = build(tmp_path, head + loop + rows + "M1 0.2 0.15 0.25 ? dum\nM2 0.25 0.1 0.2 ? dum\n")
    plain = build(tmp_path, head.replace("1234567", "1234568") + loop + rows, name="p.cif")
    assert [s.label for s in e.structure.sites] == ["Rh1", "C1"] and e.formula == plain.formula
    with pytest.raises(entry.EntryError, match="not an element"):  # the same rows without the flag: not certain
        build(tmp_path, head + loop + rows + "M1 0.2 0.15 0.25 ? ?\n", name="n.cif")


def test_element_given_by_a_correction_is_not_read_again(tmp_path, monkeypatch):
    """Rhenium labelled RE1 in a CIF without type column whose formula lacks Re (COD 7027327): the correction's
    element stands, although its label spells an element missing from the formula."""
    head = cif("", sg="P 1", a=5, b=6, c=7, extra="_chemical_formula_sum 'C N O'").split("loop_\n_atom_site_label")[0]
    loop = "loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
    text = head + loop + "RE1 0.1 0.2 0.3\nN1 0.3 0.1 0.2\nO1 0.2 0.3 0.1\n"
    with pytest.raises(entry.EntryError, match="RE1"):
        build(tmp_path, text, name="plain.cif")
    fixes = {"1234567": {"source": "s", "note": "n", "sites": {"RE1": {"species": {"Re": 1}}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    assert build(tmp_path, text).formula == "Re1.00 N1.00 O1.00"


def test_omitted_site_is_left_out(tmp_path, monkeypatch):
    """A ring centroid listed as an atom without a type column (COD 2003591, 'Cen', read as cerium)."""
    rows = "C1 C 0.1 0.2 0.3 1 0.02\nO1 O 0.3 0.1 0.2 1 0.02\n"
    fixes = {"1234567": {"source": "s", "note": "n", "sites": {"Cen": {"omit": True}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, cif(rows + "Cen Ce 0.2 0.15 0.25 1 0.03\n", sg="P 1", a=5, b=6, c=7))
    plain = build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7, cod=1234568), name="p.cif")
    assert [s.label for s in e.structure.sites] == ["C1", "O1"] and e.formula == plain.formula


def test_omitted_types_are_left_out(tmp_path, monkeypatch):
    """Interatomic scatterers for bonding density (COD 2017919, types IS4, IS5) are no atoms."""
    rows = "C1 C 0.1 0.2 0.3 1 0.02\nO1 O 0.3 0.1 0.2 1 0.02\nX1 IS4 0.2 0.15 0.25 0.8 0.03\n"
    fixes = {"1234567": {"source": "s", "note": "n", "omit_types": ["IS4"]}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7))
    plain = build(tmp_path, cif(rows.split("X1")[0], sg="P 1", a=5, b=6, c=7, cod=1234568), name="p.cif")
    assert [s.label for s in e.structure.sites] == ["C1", "O1"] and e.formula == plain.formula
    assert [x.intensity for x in e.lines] == pytest.approx([x.intensity for x in plain.lines], rel=1e-12)
    fixes["1234567"]["omit_types"] = ["IS5"]
    with pytest.raises(KeyError, match="no atom type"):
        build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7), name="bad.cif")


def test_correction_type_missing_from_the_cif_is_an_error(tmp_path, monkeypatch):
    fixes = {"1234567": {"source": "s", "note": "n", "types": {"CE": {"C": 1, "H": 2}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    with pytest.raises(KeyError, match="no atom type"):
        build(tmp_path, cif("O1 O 0.1 0.2 0.3 1 0.01\n", sg="P 1", a=5, b=6, c=7))


def _with_ops(text, ops):
    return text.replace("loop_\n_atom_site_label", "loop_\n_symmetry_equiv_pos_as_xyz\n" + "\n".join(ops)
                        + "\nloop_\n_atom_site_label", 1)


def test_listed_operations_with_a_shifted_origin_beat_a_readable_symbol(tmp_path):
    """P 1 21/m 1 with the origin moved by b/4 (COD 2014900 and ~90 others): the atoms lie on the listed mirror
    at y = 0; the symbol's own operations (mirror at y = 1/4) would double them."""
    ops = ["x,y,z", "-x,y+1/2,-z", "-x,-y+1/2,-z", "x,-y,z"]
    text = _with_ops(cif("Pb1 Pb 0.1 0 0.3 1 0.01\nS1 S 0.35 0 0.8 1 0.01\n", sg="P 1 21/m 1", a=5, b=6, c=7, be=100),
                     ops)
    e = build(tmp_path, text)
    assert e.formula == "Pb2.00 S2.00" and e.structure.spacegroup.number == 11 and not e.symmetry_note


def test_hall_symbol_gemmi_does_not_know_beats_the_symbol(tmp_path):
    """Hall 'P 2n -2a' (no table setting) next to 'P m c n' (COD 2001219): gemmi falls back to the symbol; the
    Hall symbol's four operations decide, here confirmed by the formula x Z."""
    text = cif("Co1 Co 0.13 0.21 0.34 1 0.01\nBr1 Br 0.37 0.08 0.71 1 0.01\n", sg="P m c n", a=5, b=6, c=7,
               extra="_symmetry_space_group_name_Hall 'P 2n -2a'\n_chemical_formula_sum 'Br Co'\n"
                     "_cell_formula_units_Z 4")
    e = build(tmp_path, text)
    assert e.structure.spacegroup.number == 26 and e.formula == "Co4.00 Br4.00"


@pytest.mark.parametrize("z, expected", [(4, "Pb4.00 S4.00"), (3, None)])
def test_incomplete_operations_fall_back_to_the_symbol(tmp_path, monkeypatch, z, expected):
    """Only the two P21 operations of P 1 21/c 1 listed: the symbol gives the CIF's formula x Z, the operations
    half of it. When neither agrees (Z = 3), the CIF is skipped."""
    text = _with_ops(cif("Pb1 Pb 0.1 0.2 0.3 1 0.01\nS1 S 0.35 0.1 0.8 1 0.01\n", sg="P 1 21/c 1", a=5, b=6, c=7,
                         be=100, extra=f"_chemical_formula_sum 'Pb S'\n_cell_formula_units_Z {z}"),
                     ["x,y,z", "-x,y+1/2,-z+1/2"])
    if expected is None:
        with pytest.raises(entry.EntryError, match="disagree"):
            build(tmp_path, text)
        # a correction that does not touch the symmetry keeps the check
        fixes = {"1234567": {"source": "s", "note": "n", "sites": {"Pb1": {"occupancy_factor": 1.0}}}}
        monkeypatch.setattr(entry, "_corrections", lambda: fixes)
        with pytest.raises(entry.EntryError, match="disagree"):
            build(tmp_path, text, name="c.cif")
    else:
        e = build(tmp_path, text)
        assert e.formula == expected and e.structure.spacegroup.number == 14 and "symbol" in e.symmetry_note


def test_correction_symmetry_replaces_the_cif_hall_symbol(tmp_path, monkeypatch):
    """Operations from the article replace every symmetry statement of the CIF, also a Hall symbol gemmi knows."""
    ops = ["x,y,z", "-x+1/2,y+1/2,-z+1/2", "-x+1/2,-y,-z", "x,-y+1/2,z+1/2"]  # P 1 21/c 1, origin at 1/4,0,0
    rows = "Pb1 Pb 0.1 0.2 0.3 1 0.01\nS1 S 0.3 0.1 0.4 1 0.02\n"
    text = cif(rows, sg="P 1 21/m 1", a=5, b=6, c=7, be=100,
               extra="_symmetry_space_group_name_Hall '-P 2yb'")
    fixes = {"1234567": {"source": "s", "note": "n", "operations": ops}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, text)
    assert e.structure.spacegroup.hm == "P 1 21/c 1"
    fixes["1234567"] = {"source": "s", "note": "n", "spacegroup": "P 1 21/c 1"}
    assert build(tmp_path, text, name="sg.cif").structure.spacegroup.hm == "P 1 21/c 1"


def test_scattering_only_corrections_do_not_count_as_atoms(tmp_path, monkeypatch):
    rows = "Si1 Si 0.1 0.2 0.3 1 0.01\nO1 O 0.3 0.1 0.2 1 0.02\n"
    plain = build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7))
    fixes = {"1234567": {"source": "s", "note": "n", "scattering_only": True,
                         "sites": {"O1": {"occupancy_factor": 1.15}},
                         "add": [{"label": "T3", "xyz": [0.6, 0.5, 0.4], "species": {"Si": 0.15}}]}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7), name="c.cif")
    assert e.formula == plain.formula and e.density == pytest.approx(plain.density)
    # the pattern is that of the refined model: O1 at 1.15, T3 added with B = 0.5
    u = sites.DEFAULT_B / (8 * math.pi**2)
    model = build(tmp_path, cif(f"Si1 Si 0.1 0.2 0.3 1 0.01\nO1 O 0.3 0.1 0.2 1.15 0.02\nT3 Si 0.6 0.5 0.4 0.15 {u!r}\n",
                                sg="P 1", a=5, b=6, c=7, cod=1234568), name="m.cif")  # no correction
    assert [(x.d, x.hkl) for x in e.lines] == [(x.d, x.hkl) for x in model.lines]
    assert [x.intensity for x in e.lines] == pytest.approx([x.intensity for x in model.lines], rel=1e-12)
    assert [x.intensity for x in e.lines] != pytest.approx([x.intensity for x in plain.lines])


def test_stand_in_atoms_scatter_but_are_not_in_the_compound(tmp_path, monkeypatch):
    """Guest molecules refined as spherical stand-in atoms (melanophlogite, COD 1544611: F and Mg in the cages)."""
    rows = "Si1 Si 0.1 0.2 0.3 1 0.01\nO1 O 0.3 0.1 0.2 1 0.02\nM1 M 0.6 0.5 0.4 1 0.05\n"
    fixes = {"1234567": {"source": "s", "note": "n", "sites": {"M1": {"species": {"Mg": 0.9}, "stand_in": True}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7))
    model = build(tmp_path, cif(rows.replace("M1 M 0.6 0.5 0.4 1", "M1 Mg 0.6 0.5 0.4 0.9"), sg="P 1", a=5, b=6, c=7,
                                cod=1234568), name="m.cif")
    assert e.formula == "Si1.00 O1.00" and e.elements == ["Si", "O"]
    assert [x.intensity for x in e.lines] == pytest.approx([x.intensity for x in model.lines], rel=1e-12)
    assert e.density < model.density
    fixes["1234567"]["sites"]["M1"] = {"stand_in": True}  # a stand-in needs the species it stands in with
    with pytest.raises(ValueError):
        build(tmp_path, cif(rows, sg="P 1", a=5, b=6, c=7), name="bad.cif")


def test_spherical_shell_scatterer(tmp_path, monkeypatch):
    """A water layer refined as a hollow sphere (faujasite, COD 1545186; FullProf SASH): its atoms times
    sin(qR)/(qR) at the centre, counted as water in the formula, absent from the structure record."""
    rows = "Na1 Na 0 0 0 1 0.01\nSPS SASH 0.5 0.5 0.5 0.8 0.02\n"
    fixes = {"1234567": {"source": "s", "note": "n",
                         "sites": {"SPS": {"shell": {"radius": 3.0, "per_cell": {"O": 1.5, "H": 3.0}}}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, cif(rows, sg="P 1", a=10, b=11, c=12))
    assert e.formula == "Na1.00 O1.50 H3.00" and e.shell_sites == (1,)
    assert [a.label for a in structure.build(e, phase.HS4).atoms] == ["Na1"]
    hkl = np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1], [2, 1, 0]])
    d = 1 / np.sqrt((hkl / [10, 11, 12]) ** 2 @ [1, 1, 1])
    got = pattern._intensities(e.structure, hkl, d, pattern.CU_KA1, e.scattering)
    s2 = 1 / (4 * d**2)
    na = (np.array([gemmi.Element("Na").it92.calculate_sf(x) for x in s2]) + complex(*pattern._anomalous("Na", pattern.CU_KA1)))
    na *= np.exp(-8 * math.pi**2 * 0.01 * s2)
    qr = 2 * math.pi * 3.0 / d
    atoms = sum(n * (np.array([gemmi.Element(el).it92.calculate_sf(x) for x in s2])
                     + complex(*pattern._anomalous(el, pattern.CU_KA1))) for el, n in (("O", 1.5), ("H", 3.0)))
    shell = atoms * np.sin(qr) / qr * np.exp(-2 * math.pi**2 * 0.02 / d**2) * np.exp(1j * math.pi * hkl.sum(axis=1))
    th = np.arcsin(pattern.CU_KA1 / (2 * d))
    want = 2 * np.abs(na + shell) ** 2 * (1 + np.cos(2 * th) ** 2) / (np.sin(th) ** 2 * np.cos(th))
    assert got == pytest.approx(want, rel=1e-6)  # gemmi's IT92 sum is float32


def test_group_correction_replaces_a_site_by_its_atoms(tmp_path, monkeypatch):
    fixes = {"1234567": {"source": "s", "note": "n", "sites": {"C1": {"group": {"C": 1, "H": 3}}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, cif("O1 O 0.1 0.2 0.3 1 0.01\nC1 C 0.3 0.1 0.2 0.5 0.02\n", sg="P 1", a=5, b=6, c=7))
    assert [s.label for s in e.structure.sites] == ["O1", "C1", "C1_H1", "C1_H2", "C1_H3"]
    assert e.formula == "O1.00 C0.50 H1.50" and e.site_rows == [0, 1, 1, 1, 1]


@pytest.mark.parametrize("label, formula, element", [
    ("PB", "Ca5 F O12 P3", "P"),        # AMCSD: P on site B, not lead (COD 9016317)
    ("YB", "Ba C5 Ca Na3 O15 Sr2 Y0.1", "Y"),  # Y on site B, not ytterbium (COD 9013840)
    ("OS1", "Al3 Ca Na3 O14 S Si3", "O"),  # O on site S1, not osmium (COD 9012653)
    ("HO2a", "C38 H51 Mn2 N6 O14", "H"),   # H on O2, not holmium (COD 7125504)
    ("CN", "C K N", None),             # cyanide: C and N, ambiguous
    ("TM1", "Al12 Cr Fe2", None),      # a mixed Cr/Fe site, not thulium (COD 2101903)
    ("Co1", "Bi Fe O3", "Co"),         # a dopant the formula leaves out
])
def test_two_capitals_read_against_the_formula(tmp_path, label, formula, element):
    text = HEAD.format(cod=1234567, extra=f"_chemical_formula_sum '{formula}'", a=5, b=6, c=7, al=90, be=90, ga=90,
                       sg="P 1")
    text += ("loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
             f"Si9 0.1 0.2 0.3\n{label} 0.3 0.1 0.2\n")
    text = text.replace("Si9", "Ca9") if "Ca" in formula else text
    if element is None:
        with pytest.raises(entry.EntryError, match="not an element"):
            build(tmp_path, text)
    else:
        assert build(tmp_path, text).structure.sites[1].element.name == element


def test_correction_of_the_setting(tmp_path, monkeypatch):
    """A correction can state the setting the coordinates are really in (COD 5910081: origin choice 1)."""
    rows = "Hg1 Hg 0 0 0 1 0.01\nK1 K 0.625 0.625 0.625 1 0.02\nCN CN 0.387 0.387 0.387 1 0.03\n"
    text = cif(rows, sg="F d -3 m :2", a=12.76, extra="_chemical_formula_sum 'C Hg K2 N4'")
    fixes = {"1234567": {"source": "s", "note": "n", "spacegroup": "F d -3 m:1", "types": {"CN": {"C": 1, "N": 1}}}}
    monkeypatch.setattr(entry, "_corrections", lambda: fixes)
    e = build(tmp_path, text)
    assert e.formula.split() == ["Hg8.00", "K16.00", "C32.00", "N32.00"]


def test_atom_type_description_names_the_element(tmp_path):
    """Type TL (a 3% thallium substitution the formula leaves out, COD 5000363): the CIF's atom_type loop
    describes it as Tl. Without that description the CIF is skipped."""
    rows = "Zn1 ZN 0.1 0.2 0.3 0.97 0.02\nTL1 TL 0.11 0.19 0.29 0.03 0.02\nN1 N 0.4 0.1 0.2 1 0.02\n"
    types = "loop_\n_atom_type_symbol\n_atom_type_description\nZN Zn\nTL Tl\nN N\n"
    text = cif(rows, sg="P 1", a=5, b=6, c=7, extra="_chemical_formula_sum 'N Zn'", tail=types)
    assert [s.element.name for s in build(tmp_path, text).structure.sites] == ["Zn", "Tl", "N"]
    with pytest.raises(entry.EntryError):
        build(tmp_path, text.replace("TL Tl\n", "TL ?\n"), name="nodesc.cif")


def test_label_that_only_starts_like_an_element_is_not_that_element(tmp_path):
    """No type column: 'Cen' is a phenyl ring centroid (COD 2003591), which gemmi reads as cerium; the formula
    has no cerium, and nothing says which atom it is: the CIF is skipped. 'Ce1' in the same CIF would stand."""
    text = HEAD.format(cod=1234567, extra="_chemical_formula_sum 'C6 H5 B F3 K'", a=5, b=6, c=7, al=90, be=90,
                       ga=90, sg="P 1")
    text += "loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
    text += "K 0.1 0.2 0.3\nC(1) 0.3 0.1 0.2\nCen 0.35 0.15 0.25\n"
    with pytest.raises(entry.EntryError, match="read as Ce"):
        build(tmp_path, text)
    assert build(tmp_path, text.replace("Cen ", "Ce1 "), name="ce1.cif").elements == ["K", "C", "Ce"]


def test_label_without_type_names_its_element(tmp_path):
    """No type column: HC12 is a hydrogen on C12 (COD 2101850), not H with a C12 group."""
    text = HEAD.format(cod=1234567, extra="_chemical_formula_sum 'C12 H10 Pb S2'", a=5, b=6, c=7, al=90, be=90,
                       ga=90, sg="P 1")
    text += "loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
    text += "Pb 0.1 0.2 0.3\nC12 0.3 0.1 0.2\nHC12 0.35 0.15 0.25\n"
    assert [s.element.name for s in build(tmp_path, text).structure.sites] == ["Pb", "C", "H"]


MONO_ATOMS = "C1 C 0.123 0.231 0.345 1 0.02\nO1 O 0.412 0.118 0.271 1 0.03\nN1 N 0.301 0.402 0.067 1 0.025\n"


@pytest.mark.parametrize("sg, a, b, c, be, want", [
    # COD 2000028: P 1 21/n 1 -> P 1 21/c 1 with c' = a + c, as COD24 (not c - a: 21.03 A, 140.0 deg)
    ("P 1 21/n 1", 12.942, 9.264, 13.877, 103.25, (12.942, 9.264, 16.666, 125.85)),
    # COD 2001062: P 1 21/a 1 -> P 1 21/c 1 with beta obtuse (not 69.4 deg)
    ("P 1 21/a 1", 14.47, 19.749, 17.50, 110.6, (17.50, 19.749, 14.47, 110.6)),
    # COD 1517594: already P 1 21/c 1, c' = c + 2a
    ("P 1 21/c 1", 12.2756, 14.1032, 25.8603, 129.437, (12.2756, 14.1032, 21.562, 112.14)),
    # COD 1516255: acute beta made obtuse
    ("P 1 21/c 1", 4.6011, 12.7317, 15.9122, 78.949, (4.6011, 12.7317, 15.9122, 101.051)),
    # COD 1000420: C 1 2/c 1, c' = c + a
    ("C 1 2/c 1", 7.7526, 7.5228, 7.4477, 124.081, (7.7526, 7.5228, 7.132, 120.12)),
    # COD 1000140: P 1 21 1, a < c
    ("P 1 21 1", 9.0652, 7.1468, 5.4700, 98.782, (5.4700, 7.1468, 9.0652, 98.782)),
    # orthorhombic, axes permuted where the symbol allows: COD 1000001 all, 1000047 cyclic, 1000073 a and b only
    ("P 21 21 21", 48.48, 21.72, 10.74, 90, (10.74, 21.72, 48.48, 90)),
    ("P b c a", 18.251, 8.814, 5.181, 90, (5.181, 18.251, 8.814, 90)),
    ("C c c a :2", 8.789, 8.768, 12.958, 90, (8.768, 8.789, 12.958, 90)),
    ("P n m a", 5.44, 7.639, 5.38, 90, (5.44, 7.639, 5.38, 90)),  # no other cell keeps P n m a
])
def test_cell_choice_as_cod24(tmp_path, monkeypatch, sg, a, b, c, be, want):
    """Monoclinic and orthorhombic entries are stored in the cell of the reference setting that the official COD
    databases choose (monoclinic: beta >= 90, the shortest a + c, then a; orthorhombic: the shortest a, then b);
    the pattern does not change."""
    text = cif(MONO_ATOMS, sg=sg, a=a, b=b, c=c, be=be)
    e = build(tmp_path, text)
    assert e.cell[:3] == pytest.approx(want[:3], abs=2e-3) and e.cell[4] == pytest.approx(want[3], abs=0.02)
    assert e.volume == pytest.approx(gemmi.UnitCell(a, b, c, 90, be, 90).volume, rel=1e-9)
    monkeypatch.setattr(setting, "reduce_cell", lambda st: (st, np.eye(3), False))
    fixed = build(tmp_path, text, name="fixed.cif")
    assert [(round(x.d, 8), round(x.intensity, 6)) for x in e.lines] == \
        [(round(x.d, 8), round(x.intensity, 6)) for x in fixed.lines]
    assert e.density == pytest.approx(fixed.density, rel=1e-12) and e.iic == pytest.approx(fixed.iic, rel=1e-9)


def test_triclinic_cell_is_niggli_reduced_as_cod24(tmp_path, monkeypatch):
    """COD 1000034: P -1 stored in the Niggli-reduced cell, as in the official COD databases."""
    text = cif(MONO_ATOMS, sg="P -1", a=8.173, b=12.869, c=14.165, al=93.113, be=115.913, ga=91.261)
    e = build(tmp_path, text)
    assert e.cell == pytest.approx((8.173, 12.869, 12.894, 85.779, 81.154, 88.739), abs=2e-3)
    monkeypatch.setattr(setting, "reduce_cell", lambda st: (st, np.eye(3), False))
    fixed = build(tmp_path, text, name="fixed.cif")
    assert [(round(x.d, 8), round(x.intensity, 6)) for x in e.lines] == \
        [(round(x.d, 8), round(x.intensity, 6)) for x in fixed.lines]


def test_data_block_name_gives_the_cod_id(tmp_path):
    """COD 1010542.cif (data_1010542) carries _cod_database_code 9016579, another entry's id: the data block name
    decides, else that entry would replace the real 9016579."""
    text = cif(NACL_ATOMS).replace("_cod_database_code 1234567", "_cod_database_code 7654321")
    assert build(tmp_path, text).cod_id == 1234567
    named = text.replace("data_1234567", "data_nacl")  # a block name that is no id: the tag gives it
    assert build(tmp_path, named, name="n.cif").cod_id == 7654321
