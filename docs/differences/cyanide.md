# Cyanide pseudo-atom `CN`

**COD24:** reads `CN` as carbon, so the nitrogen is lost.
**hsrdb-tools:** C and N at the CN site, and for the three K₂M(CN)₄ the coordinates in origin choice 1 ([corrections.md](../corrections.md)).

## The four CIFs

| COD | Compound | Problem |
|---|---|---|
| [5910081](https://www.crystallography.net/cod/5910081.cif) | K₂Hg(CN)₄ | one `CN` site; coordinates in origin choice 1 while the CIF says F d −3 m origin choice 2, which puts 16 Hg in the cell instead of 8 |
| [5910114](https://www.crystallography.net/cod/5910114.cif) | K₂Zn(CN)₄ | same |
| [5910134](https://www.crystallography.net/cod/5910134.cif) | K₂Cd(CN)₄ | same |
| [5910129](https://www.crystallography.net/cod/5910129.cif) | KCN | one `CN` site |

COD24 unit-cell formulas: K8 Hg16 C32, K8 Zn16 C32, K8 Cd16 C32 and K4 C4, all without nitrogen.

## Background

All four come from Wyckoff's 1931 compilation *The Structure of Crystals*, which merged C and N into one site. The original determinations placed them separately: Dickinson, J. Am. Chem. Soc. 44, 774 (1922) for the three K₂M(CN)₄, and Bozorth, J. Am. Chem. Soc. 44, 317 (1922) for KCN. Room-temperature KCN is orientationally disordered, as later neutron studies showed.

COD has correct entries for the same compounds, and hsrdb-tools includes them: [1010079](https://www.crystallography.net/cod/1010079.cif), [1010080](https://www.crystallography.net/cod/1010080.cif) and [1010081](https://www.crystallography.net/cod/1010081.cif) (Dickinson's data, C and N separate, origin choice 1), and for KCN [1541430](https://www.crystallography.net/cod/1541430.cif) and [1541550](https://www.crystallography.net/cod/1541550.cif) (neutron diffraction).

## Options considered

| Option | Effect | Decision |
|---|---|---|
| Skip the four | superseded by correct COD entries, which are included | earlier behaviour |
| Read `CN` as carbon (COD24's way) | formula and pattern lack nitrogen; the three K₂M(CN)₄ keep the wrong origin | rejected |
| C and N on the same site, same occupancy (and origin choice 1 for the K₂M(CN)₄) | the compilation's own model: Wyckoff describes the cyanide group as a single spherical scatterer | **chosen** |
| Split C and N along the site's 3-fold axis with a C–N distance of 1.15 Å | chemically right bond, but invents coordinates the CIF does not give | rejected |

Decided from Wyckoff's supplement ([archive.org](https://archive.org/details/structureofcryst030914mbp), pp. 15 and 49): the CN group acts geometrically as a single atom, a sphere about the size of a bromide ion, and separate C and N positions could not be established. The corrected K₂Hg(CN)₄ has the same content and density as Dickinson's split model (COD 1010081), and its strongest lines agree with it.
