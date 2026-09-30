"""Calculated powder stick patterns from crystal structures, in the form HighScore stores them."""

import math
from dataclasses import dataclass
from functools import lru_cache

import gemmi
import numpy as np

from . import sites

CU_KA1 = 1.540598
MAX_LINES = 203


@dataclass
class Line:
    d: float
    intensity: float  # relative, strongest line = 1000
    hkl: tuple[int, int, int]


def d_for_two_theta(two_theta_deg, wavelength=CU_KA1):
    return wavelength / (2 * math.sin(math.radians(two_theta_deg) / 2))


@lru_cache(maxsize=None)
def _anomalous(element, wavelength):
    z = gemmi.Element(element).atomic_number
    return gemmi.cromer_liberman(z=z, energy=12398.4193 / wavelength)


def _equivalents(hkl, spacegroup):
    """Integer keys of all hkl equivalent under the group's point operations plus Friedel pairs, shape (n, 2*ops)."""
    rots = np.array([op.rot for op in spacegroup.operations().sym_ops]) // gemmi.Op.DEN
    eq = np.einsum("nj,rjk->nrk", hkl, rots)  # h' = h R
    eq = np.concatenate([eq, -eq], axis=1)
    keys = (eq[..., 0] + 512) * 1048576 + (eq[..., 1] + 512) * 1024 + (eq[..., 2] + 512)
    keys.sort(axis=1)
    return keys


def _multiplicities(hkl, spacegroup):
    """Number of distinct reflections equivalent to each hkl under the Laue group."""
    return 1 + np.count_nonzero(np.diff(_equivalents(hkl, spacegroup), axis=1), axis=1)


def _intensities(st, hkl, d, wavelength, images, chunk=2048):
    """mult * |F|^2 * Lorentz-polarisation for each reflection (no monochromator).

    |F|^2 is the Friedel average (|F(h)|^2 + |F(-h)|^2) / 2: with anomalous scattering the two differ in
    non-centrosymmetric structures, and a powder line holds both (the multiplicity counts Friedel mates).
    """
    elements = sorted({s.element.name for s in st.sites})
    el_index = np.array([elements.index(s.element.name) for s in st.sites])[images.site]
    occ = np.array([s.occ for s in st.sites])[images.site]
    frac = np.array(st.cell.frac.mat.tolist())
    anomalous = np.array([complex(*_anomalous(el, wavelength)) for el in elements])  # f' + i f''
    F2 = np.empty(len(hkl))
    for lo in range(0, len(hkl), chunk):
        h = hkl[lo:lo + chunk]
        stol2 = 1 / (4 * d[lo:lo + chunk] ** 2)
        f = np.stack([_f0(el, stol2) for el in elements], axis=1) + anomalous  # (n, elements)
        s_cart = h @ frac  # reciprocal-space vectors, |s| = 1/d
        t = np.exp(-2 * math.pi**2 * np.einsum("ni,mij,nj->nm", s_cart, images.ucart, s_cart))
        phase = np.exp(2j * math.pi * (h @ images.xyz.T))
        w = occ * f[:, el_index] * t
        fp, fm = np.sum(w * phase, axis=1), np.sum(w * phase.conj(), axis=1)  # F(h), F(-h)
        for sh in images.shells:  # atoms on spherical shells: sum(n f) * sin(qR)/(qR), q = 2 pi / d
            dd = d[lo:lo + chunk]
            qr = 2 * math.pi * sh.radius / dd
            atoms = sum(n / len(sh.xyz) * (_f0(el, stol2) + complex(*_anomalous(el, wavelength)))
                        for el, n in sh.per_cell.items())
            fs = atoms * np.sin(qr) / qr * np.exp(-2 * math.pi**2 * sh.u_iso / dd**2)
            ph = np.exp(2j * math.pi * (h @ sh.xyz.T)).sum(axis=1)
            fp, fm = fp + fs * ph, fm + fs * ph.conj()
        F2[lo:lo + chunk] = (np.abs(fp) ** 2 + np.abs(fm) ** 2) / 2
    theta = np.arcsin(wavelength / (2 * d))
    lp = (1 + np.cos(2 * theta) ** 2) / (np.sin(theta) ** 2 * np.cos(theta))  # powder Lorentz-polarisation
    return _multiplicities(hkl, st.spacegroup) * F2 * lp


@lru_cache(maxsize=None)
def _it92(element):
    c = gemmi.Element(element).it92.get_coefs()
    return np.array(c[:4]), np.array(c[4:8]), c[8]


def _f0(element, stol2):
    """IT92 atomic form factor: sum of a_i exp(-b_i sin^2(theta)/lambda^2) + c."""
    a, b, c = _it92(element)
    return np.exp(-np.outer(stol2, b)) @ a + c


MERGE_TWO_THETA = 0.005  # lines closer than this (degrees 2-theta, Cu Ka1) are stored as one


def calculate(st, two_theta_max=90.0, wavelength=CU_KA1, images=None):
    """All lines up to two_theta_max, sorted by decreasing d, and the absolute intensity of the strongest line.

    st is a gemmi.SmallStructure with chemical (as-published) occupancies; images is
    sites.unit_cell_images(st), computed here if not given. Each reflection within MERGE_TWO_THETA of the
    previous one joins its line (so a line can span more than MERGE_TWO_THETA); a line carries the summed
    intensity and the d and hkl of its strongest reflection. Reflections of zero intensity are left out.
    """
    return _calculate(st, d_for_two_theta(two_theta_max, wavelength), wavelength,
                      images if images is not None else sites.unit_cell_images(st))


def _calculate(st, dmin, wavelength, images):
    hkl = np.asarray(gemmi.make_miller_array(st.cell, st.spacegroup, dmin), dtype=np.int64)
    if len(hkl) == 0:
        return [], 0.0
    if np.abs(hkl).max() >= 512:  # _equivalents packs indices into 10 bits
        raise ValueError("unit cell too large for the reflection range")
    frac = np.array(st.cell.frac.mat.tolist())
    d = 1 / np.linalg.norm(hkl @ frac, axis=1)
    return merge_lines(hkl, d, _intensities(st, hkl, d, wavelength, images), wavelength)


def merge_lines(hkl, d, inten, wavelength=CU_KA1):
    """Reflections (hkl, d, intensity) as lines by decreasing d, and the intensity of the strongest line (see
    calculate)."""
    if not len(inten) or not inten.max() > 0:
        return [], 0.0
    # reflections extinct by the atom positions (numerically zero) are no lines and must not join two lines
    present = inten > 1e-12 * inten.max()
    hkl, d, inten = hkl[present], d[present], inten[present]
    two_theta = np.degrees(2 * np.arcsin(wavelength / (2 * d)))
    merged = []  # [two_theta of last member, total intensity, strongest intensity, d, hkl]
    for i in np.argsort(two_theta, kind="stable"):
        t, inte = float(two_theta[i]), float(inten[i])
        if merged and t - merged[-1][0] <= MERGE_TWO_THETA:
            m = merged[-1]
            m[0] = t
            m[1] += inte
            if inte > m[2] * (1 + 1e-9):  # equal within rounding: the first (lower-angle) reflection stays
                m[2], m[3], m[4] = inte, float(d[i]), tuple(int(v) for v in hkl[i])
        else:
            merged.append([t, inte, inte, float(d[i]), tuple(int(v) for v in hkl[i])])
    imax = max(m[1] for m in merged)
    if imax <= 0:
        return [], 0.0
    lines = [Line(m[3], 1000 * m[1] / imax, m[4]) for m in merged if 1000 * m[1] / imax >= 1e-4]
    return lines, imax


def stick_pattern(st, wavelength=CU_KA1, images=None):
    """Lines to store for a structure: up to 90 deg 2-theta (140 deg if fewer than 10), at most MAX_LINES by decreasing d.

    Also returns the absolute intensity of the strongest line, for I/Ic.
    """
    if images is None:
        images = sites.unit_cell_images(st)
    d90 = d_for_two_theta(90.0, wavelength)
    # start from the d-spacing that should give a few times MAX_LINES reflections, widening as needed
    n_ops = len(st.spacegroup.operations())
    dmin = max(d90, (4 * math.pi * st.cell.volume / (3 * 4 * MAX_LINES * n_ops)) ** (1 / 3))
    while True:
        lines, imax = _calculate(st, dmin, wavelength, images)
        if len(lines) >= MAX_LINES or dmin <= d90:
            break
        dmin = max(d90, dmin * 0.8)
    if len(lines) < 10:
        lines, imax = calculate(st, 140.0, wavelength, images)
    lines = lines[:MAX_LINES]
    top = max((x.intensity for x in lines), default=0)
    if 0 < top < 1000:  # the strongest line fell outside the kept range: rescale to 1000
        imax *= top / 1000
        lines = [Line(x.d, 1000 * x.intensity / top, x.hkl) for x in lines]
    return lines, imax
