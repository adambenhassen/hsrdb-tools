# Dummy atoms

**COD24:** atoms flagged as dummies (`_atom_site_calc_flag dum`) are treated as real atoms.
**hsrdb-tools:** they do not scatter and are left out of the formula, density and element search filters. They stay in the structure record, flagged as dummy atoms, as in the official records.

## Why

COD CIFs use the dummy flag for two kinds of site, neither of which is an atom at that position:

- **Geometric points** a refinement program uses, such as ring centroids or metal–ring midpoints (e.g. `C(100)`, `Cg1`). They hold no electrons.
- **Atoms whose position was never determined**, which COD lists with coordinates −1, −1, −1 so that the formula stays complete. The atoms exist, but placing them at −1, −1, −1 (the cell origin) puts them where they are not.

In both cases, scattering from the flagged site adds intensity the crystal does not produce there. hsrdb-tools calculates the pattern of the atoms whose positions were published, in the same way that hydrogen atoms not located in a structure are simply absent.

## Evidence

438 COD CIFs contain dummy atoms (`dev/cif_features.py`). Three examples:

| COD | CIF formula | Dummy sites | COD24 unit-cell formula | hsrdb-tools |
|---|---|---|---|---|
| [1010599](https://www.crystallography.net/cod/1010599.cif) | C6 K2 N6 Pt S6 | C1 and N1 at −1, −1, −1 (positions not determined) | Pt1.00 K2.00 S6.00 **C1.00 N1.00** | Pt1.00 K2.00 S6.00 |
| [1101029](https://www.crystallography.net/cod/1101029.cif) | Al3 H2 K O12 Si3 | H1 at −1, −1, −1 | K4.00 Si12.00 Al12.00 O48.00 **H1.00** | K4.00 Si12.00 Al12.00 O48.00 |
| [2005855](https://www.crystallography.net/cod/2005855.cif) | C20 H12 Ge O8 Os2 | `C(100)`, `C(200)`, marked as dummy atoms in the CIF | Os16.00 Ge8.00 O64.00 **C162.00** H96.00 | Os16.00 Ge8.00 O64.00 C160.00 H96.00 |

In 2005855, Z = 8 gives C160, as in the CIF's own formula; COD24's two extra carbons are the dummy points.

The build report notes each CIF with dummy atoms (`note: N dummy atoms do not scatter`). Test: `test_dummy_atoms_do_not_scatter`.
