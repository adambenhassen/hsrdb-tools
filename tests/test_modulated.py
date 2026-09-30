"""Modulated structures (superspace descriptions): analytic structure factors and equivalences."""

import math

import numpy as np
import pytest

from hsrdb_tools import entry, pattern, superspace

HEAD = """data_{cod}
_cod_database_code {cod}
_journal_name_full 'Test Journal'
_journal_year 2020
_cell_length_a 4.0
_cell_length_b 5.0
_cell_length_c 6.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_cell_modulation_dimension 1
loop_
_space_group_symop_operation_xyz
{ops3}
loop_
_space_group_symop_ssg_operation_algebraic
{ops4}
loop_
_cell_wave_vector_seq_id
_cell_wave_vector_x
_cell_wave_vector_y
_cell_wave_vector_z
1 {q} 0 0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
{atoms}
"""
P1 = ("x,y,z", "x1,x2,x3,x4")
P_1 = ("x,y,z\n-x,-y,-z", "x1,x2,x3,x4\n-x1,-x2,-x3,-x4")
U = 0.01


def text(atoms, tail="", group=P1, q=0.3, cod=1234567):
    return HEAD.format(cod=cod, ops3=group[0], ops4=group[1], q=q, atoms=atoms) + tail


def displace(rows):
    return ("loop_\n_atom_site_displace_Fourier_atom_site_label\n_atom_site_displace_Fourier_axis\n"
            "_atom_site_displace_Fourier_wave_vector_seq_id\n_atom_site_displace_Fourier_param_cos\n"
            "_atom_site_displace_Fourier_param_sin\n" + rows)


def assert_basic(e, reason):
    """The entry is the basic structure (as in COD24), noted with why the modulation is not included."""
    assert e.modulated is None and all(x.hkl != (0, 0, 0) for x in e.lines)
    assert reason in e.modulation and "as in COD24" in e.modulation and "modulation not included" in e.comment


def build(tmp_path, cif, name="1234567.cif"):
    p = tmp_path / name
    p.write_text(cif)
    return entry.from_cif(p)


def f_t(e, h, m):
    """|f| T of the single Cu atom at the reflection (h, m) of entry e's model."""
    model = e.modulated
    S = np.array(h, float) + m * model.q[0]
    s2 = S @ np.array(model.cell.reciprocal_metric_tensor().as_mat33().tolist()) @ S
    f = pattern._f0("Cu", np.array([s2 / 4]))[0] + complex(*pattern._anomalous("Cu", pattern.CU_KA1))
    return abs(f) * math.exp(-2 * math.pi**2 * U * s2)


def bessel(m, z):
    t = (np.arange(20000) + 0.5) / 20000 * math.pi
    return float(np.mean(np.cos(m * t - z * np.sin(t))))


def test_displacement_wave_gives_bessel_satellites(tmp_path):
    A = 0.02  # fractional amplitude along a
    e = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace(f"Cu1 x 1 0 {A}\nCu1 y 1 0 0\nCu1 z 1 0 0\n")))
    H = [(2, 1, 0, 0), (2, 1, 0, 1), (2, 1, 0, -2), (1, 0, 1, 3)]
    F = np.abs(superspace.structure_factors(e.modulated, H))
    for (h, k, l, m), got in zip(H, F):
        z = 2 * math.pi * (h + m * 0.3) * A
        assert got == pytest.approx(f_t(e, (h, k, l), m) * abs(bessel(m, z)), rel=1e-9)


def test_crenel_occupancy_gives_sinc_satellites(tmp_path):
    w, occ = 0.4, 0.4  # JANA writes the average occupancy: the site is full within the interval
    tail = ("loop_\n_atom_site_occ_special_func_atom_site_label\n_atom_site_occ_special_func_crenel_c\n"
            f"_atom_site_occ_special_func_crenel_w\nCu1 0.5 {w}\n")
    e = build(tmp_path, text(f"Cu1 Cu 0 0 0 {occ} {U}", tail))
    H = [(1, 1, 1, 0), (1, 1, 1, 1), (0, 2, 1, -2), (2, 0, 0, 3)]
    F = np.abs(superspace.structure_factors(e.modulated, H))
    for (h, k, l, m), got in zip(H, F):
        sinc = 1.0 if m == 0 else abs(math.sin(math.pi * m * w) / (math.pi * m * w))
        assert got == pytest.approx(f_t(e, (h, k, l), m) * occ * sinc, rel=1e-9)


def test_occupancy_wave(tmp_path):
    tail = ("loop_\n_atom_site_occ_Fourier_atom_site_label\n_atom_site_occ_Fourier_wave_vector_seq_id\n"
            "_atom_site_occ_Fourier_param_cos\n_atom_site_occ_Fourier_param_sin\nCu1 1 0.3 0\n")
    e = build(tmp_path, text(f"Cu1 Cu 0 0 0 0.6 {U}", tail))
    F = np.abs(superspace.structure_factors(e.modulated, [(1, 0, 0, 0), (1, 0, 0, 1), (1, 0, 0, 2)]))
    assert F[0] == pytest.approx(f_t(e, (1, 0, 0), 0) * 0.6, rel=1e-9)
    assert F[1] == pytest.approx(f_t(e, (1, 0, 0), 1) * 0.15, rel=1e-9)
    assert F[2] == pytest.approx(0, abs=1e-9)


def test_legendre_order_one_is_the_sawtooth(tmp_path):
    crenel = ("loop_\n_atom_site_occ_special_func_atom_site_label\n_atom_site_occ_special_func_crenel_c\n"
              "_atom_site_occ_special_func_crenel_w\nCu1 0.25 0.5\n")
    legendre = ("loop_\n_jana_atom_site_displace_legendre_atom_site_label\n_jana_atom_site_displace_legendre_axis\n"
                "_jana_atom_site_displace_legendre_param_order\n_jana_atom_site_displace_legendre_param_coeff\n"
                "Cu1 x 1 0.03\nCu1 y 1 -0.01\n")
    saw = ("loop_\n_atom_site_displace_special_func_atom_site_label\n_atom_site_displace_special_func_sawtooth_ax\n"
           "_atom_site_displace_special_func_sawtooth_ay\n_atom_site_displace_special_func_sawtooth_az\n"
           "_atom_site_displace_special_func_sawtooth_c\n_atom_site_displace_special_func_sawtooth_w\n"
           "Cu1 0.03 -0.01 0 0.25 0.5\n")
    a = build(tmp_path, text(f"Cu1 Cu 0 0 0 0.5 {U}", crenel + legendre), "a.cif")  # average occupancy 0.5
    b = build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", saw), "b.cif")  # sawtooth alone: occupancy 1 within the interval
    H = [(1, 0, 0, m) for m in range(-4, 5)] + [(2, 3, 1, 1)]
    np.testing.assert_allclose(superspace.structure_factors(a.modulated, H),
                               superspace.structure_factors(b.modulated, H), rtol=1e-12)


def test_symmetry_images_equal_explicit_atoms(tmp_path):
    """P-1 image of an atom with u(v) = B cos 2pi v + A sin 2pi v: at -xbar with u'(v) = -u(-v) = -B cos + A sin."""
    wave = "Cu1 x 1 0.01 0.02\nCu1 y 1 -0.015 0.005\nCu1 z 1 0 0.01\n"
    image = "Cu2 x 1 -0.01 0.02\nCu2 y 1 0.015 0.005\nCu2 z 1 0 0.01\n"
    sym = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace(wave), group=P_1), "s.cif")
    explicit = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}\nCu2 Cu -0.1 -0.2 -0.3 1 {U}",
                                    displace(wave + image)), "e.cif")
    H = [(h, k, 1, m) for h in range(-2, 3) for k in (0, 1) for m in (-2, -1, 0, 1, 2)]
    np.testing.assert_allclose(superspace.structure_factors(sym.modulated, H),
                               superspace.structure_factors(explicit.modulated, H), rtol=1e-9, atol=1e-9)


def test_atom_on_special_position_counts_once(tmp_path):
    wave = "Cu1 x 1 0 0.02\nCu1 y 1 0 0\nCu1 z 1 0 0\n"  # odd in v: compatible with the inversion centre
    sym = build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", displace(wave), group=P_1), "s.cif")
    single = build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", displace(wave)), "p.cif")
    H = [(1, 0, 0, m) for m in (-2, -1, 0, 1, 2)] + [(2, 1, 1, 1)]
    np.testing.assert_allclose(superspace.structure_factors(sym.modulated, H),
                               superspace.structure_factors(single.modulated, H), rtol=1e-9, atol=1e-9)


def test_commensurate_equals_explicit_supercell(tmp_path):
    """q = a*/2 with section t0 = 0.1: the supercell 2a holds the atom at v = 0.1 and v = 0.6."""
    A, t0 = 0.03, 0.1
    tail = displace(f"Cu1 x 1 {A} 0\nCu1 y 1 0 0\nCu1 z 1 0 0\n") + "".join(
        f"_jana_cell_commen_supercell_matrix_{i}_{j} {2 if i == j == 1 else int(i == j)}\n"
        for i in (1, 2, 3) for j in (1, 2, 3)) + f"_jana_cell_commen_t_section_1 {t0}\n"
    mod = build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", tail, q=0.5), "m.cif")
    x0 = A * math.cos(2 * math.pi * t0) / 2
    x1 = (1 + A * math.cos(2 * math.pi * (t0 + 0.5))) / 2
    supercell = ("data_1234568\n_cod_database_code 1234568\n_journal_name_full 'Test Journal'\n_journal_year 2020\n"
                 "_cell_length_a 8.0\n_cell_length_b 5.0\n_cell_length_c 6.0\n_cell_angle_alpha 90\n"
                 "_cell_angle_beta 90\n_cell_angle_gamma 90\n_symmetry_space_group_name_H-M 'P 1'\nloop_\n"
                 "_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n"
                 f"_atom_site_fract_z\n_atom_site_occupancy\n_atom_site_U_iso_or_equiv\nCu1 Cu {x0} 0 0 1 {U}\n"
                 f"Cu2 Cu {x1} 0 0 1 {U}\n")
    ref = build(tmp_path, supercell, "1234568.cif")
    assert [round(x.d, 9) for x in mod.lines] == [round(x.d, 9) for x in ref.lines]
    np.testing.assert_allclose([x.intensity for x in mod.lines], [x.intensity for x in ref.lines], rtol=1e-9)
    # per basic cell against per supercell: |F|^2 and the cell volume both double, so I/Ic agrees
    assert mod.iic == pytest.approx(ref.iic, rel=1e-9)


def test_entry_of_modulated_structure(tmp_path):
    e = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace("Cu1 x 1 0 0.05\nCu1 y 1 0 0\nCu1 z 1 0 0\n")))
    satellites = [x for x in e.lines if x.hkl == (0, 0, 0)]
    assert satellites and len(satellites) < len(e.lines)
    assert min(x.intensity for x in satellites) >= superspace.SATELLITE_LIMIT
    assert "modulated structure: pattern with satellites up to order" in e.modulation
    assert "satellite lines carry hkl 0 0 0" in e.comment
    # without the wave the same CIF gives only the average structure: main lines only, noted
    basic = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}"), "b.cif")
    assert basic.modulated is None and all(x.hkl != (0, 0, 0) for x in basic.lines)
    assert "gives only the average structure" in basic.modulation


def test_unsupported_modulation_function_falls_back_to_basic_structure(tmp_path):
    """JANA's crenel-orthogonalised functions are not reproduced; their definition alone (written by JANA also when
    unused) and zero coefficients do not count."""
    defined = ("loop_\n_jana_atom_site_crenel_ortho_func_id\n_jana_atom_site_crenel_ortho_func_c\n"
               "_jana_atom_site_crenel_ortho_func_w\n1 0.5 0.5\n")
    coefficients = ("loop_\n_jana_atom_site_displace_crenel_ortho_atom_site_label\n"
                    "_jana_atom_site_displace_crenel_ortho_axis\n_jana_atom_site_displace_crenel_ortho_id\n"
                    "_jana_atom_site_displace_crenel_ortho_param_order\n"
                    "_jana_atom_site_displace_crenel_ortho_param_coeff\nCu1 x 1 1 {}\n")
    wave = displace("Cu1 x 1 0 0.01\n")
    for i, tail in enumerate((defined + wave, defined + coefficients.format(0) + wave)):
        assert build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", tail), f"{i}.cif").modulated
    assert_basic(build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", defined + coefficients.format(0.01) + wave)),
                 "not supported: crenel_ortho")


def test_superspace_operations_must_match_space_group(tmp_path):
    cif = text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace("Cu1 x 1 0 0.01\n"), group=("x,y,z\n-x,-y,-z", P1[1]))
    assert_basic(build(tmp_path, cif), "do not match the space group")


def test_centred_basic_cell_stored_primitive(tmp_path):
    """C-1 with a displacement wave is stored in its primitive cell: same pattern and I/Ic as the four images
    listed in P1 (inversion image: -xbar, u'(v) = -u(-v); centring image: xbar + (1/2, 1/2, 0), same u)."""
    c_1 = ("x,y,z\nx+1/2,y+1/2,z\n-x,-y,-z\n-x+1/2,-y+1/2,-z",
           "x1,x2,x3,x4\nx1+1/2,x2+1/2,x3,x4\n-x1,-x2,-x3,-x4\n-x1+1/2,-x2+1/2,-x3,-x4")
    wave = "Cu1 x 1 0.01 0.02\nCu1 y 1 0 0\nCu1 z 1 -0.01 0.005\n"
    sym = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace(wave), group=c_1, q=0.31), "s.cif")
    atoms = [("Cu1", (0.1, 0.2, 0.3), 1), ("Cu2", (0.6, 0.7, 0.3), 1), ("Cu3", (-0.1, -0.2, -0.3), -1),
             ("Cu4", (0.4, 0.3, -0.3), -1)]
    rows = "".join(f"{lab} x 1 {0.01 * s} 0.02\n{lab} y 1 0 0\n{lab} z 1 {-0.01 * s} 0.005\n" for lab, _, s in atoms)
    explicit = build(tmp_path, text("\n".join(f"{lab} Cu {x} {y} {z} 1 {U}" for lab, (x, y, z), _ in atoms),
                                    displace(rows), q=0.31), "e.cif")
    assert sym.setting_changed and sym.sg_hm == "P -1"
    np.testing.assert_allclose(sorted(x.d for x in sym.lines), sorted(x.d for x in explicit.lines), rtol=1e-9)
    by_d = lambda e: [x.intensity for x in sorted(e.lines, key=lambda x: x.d)]  # noqa: E731
    np.testing.assert_allclose(by_d(sym), by_d(explicit), rtol=1e-7, atol=1e-7)
    assert sym.iic == pytest.approx(explicit.iic, rel=1e-9)


def test_wave_vectors_listed_as_components(tmp_path):
    """JANA may list each Fourier term's wave vector by its components, leaving out a column (here y)."""
    waves = ("loop_\n_atom_site_Fourier_wave_vector_seq_id\n_atom_site_Fourier_wave_vector_x\n"
             "_atom_site_Fourier_wave_vector_z\n1 0.6 0\n2 0.3 0\n")
    listed = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", waves + displace("Cu1 x 1 0 0.02\nCu1 x 2 0.01 0\n")),
                   "a.cif")
    by_order = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace("Cu1 x 2 0 0.02\nCu1 x 1 0.01 0\n")),
                     "b.cif")
    H = [(1, 0, 0, m) for m in range(-3, 4)]
    np.testing.assert_allclose(superspace.structure_factors(listed.modulated, H),
                               superspace.structure_factors(by_order.modulated, H), rtol=1e-12)


def test_wave_vector_components_with_a_column_left_out(tmp_path):
    """Older JANA CIFs omit the y column although q has a y component: wave i <= d is q_i when the listed
    components agree."""
    two_q = HEAD.replace("_cell_modulation_dimension 1", "_cell_modulation_dimension 2").replace(
        "1 {q} 0 0\n", "1 0.3 0.3 0\n2 -0.3 0.3 0\n")
    ops = ("x,y,z", "x1,x2,x3,x4,x5")
    atoms = f"Cu1 Cu 0.1 0.2 0.3 1 {U}"
    rows = "Cu1 x 1 0 0.02\nCu1 y 2 0.01 0\n"
    waves = ("loop_\n_atom_site_Fourier_wave_vector_seq_id\n_atom_site_Fourier_wave_vector_x\n"
             "_atom_site_Fourier_wave_vector_z\n1 0.3 0\n2 -0.3 0\n")
    coeffs = ("loop_\n_jana_atom_site_fourier_wave_vector_seq_id\n_jana_atom_site_fourier_wave_vector_q1_coeff\n"
              "_jana_atom_site_fourier_wave_vector_q2_coeff\n1 1 0\n2 0 1\n")
    a = build(tmp_path, two_q.format(cod=1234567, ops3=ops[0], ops4=ops[1], q=0, atoms=atoms) + waves + displace(rows),
              "a.cif")
    b = build(tmp_path, two_q.format(cod=1234567, ops3=ops[0], ops4=ops[1], q=0, atoms=atoms) + coeffs + displace(rows),
              "b.cif")
    H = [(1, 0, 0, m1, m2) for m1 in (-1, 0, 1) for m2 in (-1, 0, 1)]
    np.testing.assert_allclose(superspace.structure_factors(a.modulated, H),
                               superspace.structure_factors(b.modulated, H), rtol=1e-12)


def test_fourier_terms_by_id(tmp_path):
    """msCIF's second layout: a loop of term ids (atom, axis, wave) and a loop of their parameters."""
    by_id = ("loop_\n_atom_site_displace_Fourier_id\n_atom_site_displace_Fourier_atom_site_label\n"
             "_atom_site_displace_Fourier_axis\n_atom_site_displace_Fourier_wave_vector_seq_id\n"
             "Cu1_x_1 Cu1 x 1\nCu1_z_1 Cu1 z 1\nloop_\n_atom_site_displace_Fourier_param_id\n"
             "_atom_site_displace_Fourier_param_cos\n_atom_site_displace_Fourier_param_sin\n"
             "Cu1_x_1 0.01 0.02\nCu1_z_1 0 -0.01\n")
    a = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", by_id), "a.cif")
    b = build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace("Cu1 x 1 0.01 0.02\nCu1 z 1 0 -0.01\n")), "b.cif")
    H = [(1, 0, 1, m) for m in range(-3, 4)]
    np.testing.assert_allclose(superspace.structure_factors(a.modulated, H),
                               superspace.structure_factors(b.modulated, H), rtol=1e-12)


def test_modulation_of_unknown_atom_is_not_evaluated(tmp_path):
    assert_basic(build(tmp_path, text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace("Cu9 x 1 0 0.02\n"))),
                 "not an atom_site label")


def test_composite_is_stored_as_basic_structure(tmp_path):
    cif = text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", "_cell_subsystems_number 2\n" + displace("Cu1 x 1 0 0.02\n"))
    assert_basic(build(tmp_path, cif), "composite structure")


def test_sawtooth_interval_confines_the_atom(tmp_path):
    """Without a crenel, a sawtooth's interval is where the atom exists, at the CIF's occupancy (JANA; COD 4002590):
    the basic structure holds the average, occupancy x width."""
    saw = ("loop_\n_atom_site_displace_special_func_atom_site_label\n_atom_site_displace_special_func_sawtooth_ax\n"
           "_atom_site_displace_special_func_sawtooth_ay\n_atom_site_displace_special_func_sawtooth_az\n"
           "_atom_site_displace_special_func_sawtooth_c\n_atom_site_displace_special_func_sawtooth_w\n"
           "Cu1 0.03 0 0 0.25 0.5\n")
    crenel = ("loop_\n_atom_site_occ_special_func_atom_site_label\n_atom_site_occ_special_func_crenel_c\n"
              "_atom_site_occ_special_func_crenel_w\nCu1 0.25 0.5\n")
    a = build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", saw), "a.cif")  # occupancy 1 within the interval
    b = build(tmp_path, text(f"Cu1 Cu 0 0 0 0.5 {U}", crenel + saw.replace("0.25 0.5\n", "0.25 0.5\n")), "b.cif")
    H = [(1, 0, 0, m) for m in range(-3, 4)]
    np.testing.assert_allclose(superspace.structure_factors(a.modulated, H),
                               superspace.structure_factors(b.modulated, H), rtol=1e-12)
    assert a.structure.sites[0].occ == pytest.approx(0.5) and a.density == pytest.approx(b.density)


def test_structure_factors_listed_in_the_cif_must_be_reproduced(tmp_path):
    """A modulated CIF whose own listed F^2_calc (e.g. JANA's) the calculation does not reproduce is stored as its
    basic structure."""
    cif = text(f"Cu1 Cu 0.1 0.2 0.3 1 {U}", displace("Cu1 x 1 0 0.03\nCu1 z 1 0.01 0\n"))
    model = build(tmp_path, cif, "m.cif").modulated
    H = [(h, k, l, m) for h in range(3) for k in range(2) for l in range(2) for m in (-1, 0, 1)][1:]
    f2 = np.abs(superspace.structure_factors(model, H, 0.71073)) ** 2

    def refln(values):
        return ("_diffrn_radiation_wavelength 0.71073\nloop_\n_refln_index_h\n_refln_index_k\n_refln_index_l\n"
                "_refln_index_m_1\n_refln_F_squared_calc\n"
                + "".join(f"{h} {k} {l} {m} {v:.6g}\n" for (h, k, l, m), v in zip(H, values)))

    good = build(tmp_path, cif + refln(3 * f2), "good.cif")
    assert "reproduces the structure factors the CIF lists (R 0.0000)" in good.modulation
    wrong = f2.copy()
    wrong[:: 3] *= 2
    assert_basic(build(tmp_path, cif + refln(wrong), "bad.cif"), "structure factors the CIF lists are not reproduced")


def test_high_order_satellites_near_the_origin_are_left_out(tmp_path):
    """With q = 0.3001 a*, the order-10 satellite (-3 0 0 10) lies at d = a / 0.001; an occupancy crenel keeps it
    strong, and the Lorentz factor would make it the strongest line."""
    crenel = ("loop_\n_atom_site_occ_special_func_atom_site_label\n_atom_site_occ_special_func_crenel_c\n"
              "_atom_site_occ_special_func_crenel_w\nCu1 0.5 0.5\n")
    e = build(tmp_path, text(f"Cu1 Cu 0 0 0 0.5 {U}", crenel, q=0.3001))
    first = 1 / 0.3001 * 4.0  # d of (0 0 0 1), the largest of the main and first-order reflections
    assert max(x.d for x in e.lines) <= first * (1 + 1e-9)


def test_crenel_occupancy_reading_follows_the_cif_formula(tmp_path):
    """A crenel atom of occupancy 1 and width 0.5 holds 1 atom per cell if 1 is the average (JANA's convention) and
    0.5 if it is the value within the interval; the CIF's formula x Z decides."""
    crenel = ("loop_\n_atom_site_occ_special_func_atom_site_label\n_atom_site_occ_special_func_crenel_c\n"
              "_atom_site_occ_special_func_crenel_w\nCu1 0.5 0.5\nO1 0 0.5\n")
    atoms = f"Cu1 Cu 0 0 0 1 {U}\nO1 O 0 0 0.5 1 {U}"
    for z, occ in (("1", 1.0), ("0.5", 0.5)):
        e = build(tmp_path, text(atoms, f"_chemical_formula_sum 'Cu O'\n_cell_formula_units_Z {z}\n" + crenel),
                  f"z{z}.cif")
        assert e.structure.sites[0].occ == pytest.approx(occ)
        assert ("within the intervals" in e.modulation) == (occ == 0.5)


def test_crenel_occupancy_reading_follows_the_density_if_the_formula_does_not_decide(tmp_path):
    """Formula x Z matching neither reading (another Z convention) leaves the choice to the CIF's density."""
    crenel = ("loop_\n_atom_site_occ_special_func_atom_site_label\n_atom_site_occ_special_func_crenel_c\n"
              "_atom_site_occ_special_func_crenel_w\nCu1 0.5 0.5\n")
    inside = build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", crenel), "a.cif")
    assert "within the intervals" not in inside.modulation  # nothing decides: the average
    density = 0.5 * 63.546 * 1.66054 / 120.0  # half a Cu atom in the 4 x 5 x 6 A cell
    e = build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", "_chemical_formula_sum 'Cu'\n_cell_formula_units_Z 7\n"
                             f"_exptl_crystal_density_diffrn {density:.5f}\n" + crenel), "b.cif")
    assert "within the intervals" in e.modulation and e.structure.sites[0].occ == pytest.approx(0.5)


def test_rounded_commensurate_wave_vector_is_read_as_its_fraction(tmp_path):
    """A commensurate q written rounded (1/3 as 0.33333; COD 2310038: 13/47 as 0.276596) is the fraction: otherwise
    the satellite m q - h of highest order misses the origin by the rounding and becomes a 'line' at d ~ 10^5 A."""
    def commensurate(q, cod):
        tail = displace("Cu1 x 1 0.03 0\nCu1 y 1 0 0\nCu1 z 1 0 0\n") + "".join(
            f"_jana_cell_commen_supercell_matrix_{i}_{j} {3 if i == j == 1 else int(i == j)}\n"
            for i in (1, 2, 3) for j in (1, 2, 3)) + "_jana_cell_commen_t_section_1 0.1\n"
        return build(tmp_path, text(f"Cu1 Cu 0 0 0 1 {U}", tail, q=q, cod=cod), f"{cod}.cif")
    rounded, exact = commensurate("0.33333", 1234567), commensurate("0.3333333333", 1234568)
    assert max(x.d for x in rounded.lines) <= 12.0 + 1e-9  # the supercell 3a
    assert [round(x.d, 6) for x in rounded.lines] == [round(x.d, 6) for x in exact.lines]
    np.testing.assert_allclose([x.intensity for x in rounded.lines], [x.intensity for x in exact.lines], rtol=1e-6)
