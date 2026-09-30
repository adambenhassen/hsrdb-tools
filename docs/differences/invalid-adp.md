# Physically invalid anisotropic displacement tensors

**COD24:** uses the tensor as given.
**hsrdb-tools:** replaces it by the atom's isotropic value (`_atom_site_U_iso_or_equiv` or `_atom_site_B_iso_or_equiv`; if the CIF gives neither, the tensor's own equivalent value, trace / 3, or B = 0.5 Å² if that is not positive). A negative isotropic value is set to 0; both count in the same note.

## Why

An anisotropic displacement tensor U describes the mean-square displacement of an atom in every direction, so it must be positive definite. About 2,800 COD CIFs contain a tensor with a negative mean-square displacement in some direction, usually from a typing error or an unstable refinement. Such a tensor turns the Debye–Waller factor exp(−2π² sᵀUs) into an amplification for reflections along that direction, which gives unphysical intensities.

Criterion: the smallest eigenvalue of the tensor in Cartesian coordinates is below −10⁻⁴ Å² (a margin for rounded values). The build report notes each CIF (`note: replaced N physically invalid displacement parameters`).

## Evidence that COD24 keeps them

Median difference of strong lines from COD24 (0–1000 scale), with the tensor replaced (hsrdb-tools) and as given (calculated with cctbx), `dev/invalid_adp_check.py`:

| COD | Invalid tensors | Replaced | As given |
|---|---|---|---|
| [1000217](https://www.crystallography.net/cod/1000217.cif) | 1 | 1.4 | 0.0 |
| [9009434](https://www.crystallography.net/cod/9009434.cif) | 2 | 3.2 | 0.5 |
| [2001473](https://www.crystallography.net/cod/2001473.cif) | 14 | 1.6 | 0.5 |

Test: `test_invalid_anisotropic_tensor_falls_back_to_isotropic`.
