import sqlite3

import gemmi
import pytest

from hsrdb_tools import entry, hsrdb, pattern, phase, setting, structure

NACL = """data_1234567
_cod_database_code 1234567
_chemical_name_mineral Halite
_chemical_formula_sum 'Cl Na'
_journal_name_full 'Test Journal'
_journal_year 2020
_journal_volume 1
_journal_page_first 10
_journal_page_last 12
_cell_length_a 5.6402
_cell_length_b 5.6402
_cell_length_c 5.6402
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F m -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
Na1 Na+ 0 0 0 0.012
Cl1 Cl- 0.5 0.5 0.5 0.011
"""

# a general-position structure in a non-reference setting (P 1 21/n 1)
MONO = """data_test
_cell_length_a 7.1
_cell_length_b 9.3
_cell_length_c 8.2
_cell_angle_alpha 90
_cell_angle_beta 103.5
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1 21/n 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
C1 C 0.123 0.231 0.345 0.02
O1 O 0.412 0.118 0.271 0.03
N1 N 0.301 0.402 0.067 0.025
"""


@pytest.fixture
def nacl(tmp_path):
    p = tmp_path / "1234567.cif"
    p.write_text(NACL)
    return entry.from_cif(p)


def test_entry_fields(nacl):
    assert nacl.reference_code == 961234568
    assert nacl.sg_hm == "F m -3 m"
    assert nacl.crystal_system == "C"
    assert nacl.formula == "Na4.00 Cl4.00"
    assert nacl.subfiles == [entry.SUBFILE_INORGANIC, entry.SUBFILE_MINERAL]
    strongest = max(nacl.lines, key=lambda x: x.intensity)
    assert strongest.intensity == pytest.approx(1000)
    assert strongest.d == pytest.approx(5.6402 / 2, abs=1e-4)  # (200)
    assert [x.d for x in nacl.lines] == sorted((x.d for x in nacl.lines), reverse=True)


@pytest.mark.parametrize("fmt", [phase.HS3, phase.HS4])
def test_phase_record_round_trip(nacl, fmt):
    p = structure.build(nacl, fmt)
    payload = phase.encode_payload(p)
    decoded = phase.decode_payload(phase.decode_blob(phase.encode_blob(payload)))
    assert decoded.fmt == fmt
    assert phase.encode_payload(decoded) == payload
    assert [a.label for a in decoded.atoms] == ["Na1", "Cl1"]
    assert [a.wyckoff for a in decoded.atoms] == ["4a", "4b"]
    assert [a.charge for a in decoded.atoms] == [1, -1]


@pytest.mark.parametrize("target", ["hs3", "hs4"])
def test_database(tmp_path, nacl, target):
    path = tmp_path / f"t_{target}.hsrdb"
    db = hsrdb.create(path, target)
    hsrdb.add(db, nacl, entry.element_masks(nacl.elements), hsrdb.encode_lines(nacl.lines),
              hsrdb.phase_blob(nacl, target), target)
    hsrdb.update_prefixes(db)
    db.commit()
    db.close()
    db = sqlite3.connect(path)
    assert db.execute("pragma integrity_check").fetchone()[0] == "ok"
    assert db.execute("select Value from Properties where Key='NoOfPatterns'").fetchone()[0] == "1"
    assert db.execute("select PhaseName from PatternPhase").fetchone()[0] == "PDF:96-123-4568"
    cols = [r[1] for r in db.execute("pragma table_info(general_stored)")]
    assert ("LinesI" in cols) == (target == "hs4")


def test_reference_setting_keeps_pattern(tmp_path):
    p = tmp_path / "mono.cif"
    p.write_text(MONO)
    st = gemmi.read_small_structure(str(p))
    ref, changed = setting.to_reference(st)
    assert changed and ref.spacegroup.hm == "P 1 21/c 1"
    a, _ = pattern.calculate(st, 60)
    b, _ = pattern.calculate(ref, 60)
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert x.d == pytest.approx(y.d, abs=1e-6)
        assert x.intensity == pytest.approx(y.intensity, abs=1e-3)
