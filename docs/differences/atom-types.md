# Atom types that are not elements

**COD24:** guesses the element from the first letters of the type symbol.
**hsrdb-tools:** takes the element where it is certain (6 CIFs), builds 27 CIFs from their publication or the CIF's own statements ([corrections.md](../corrections.md)), leaves out the dummy atoms of 2 CIFs, and skips the rest (1 CIF, listed in [SKIPPED.md](../SKIPPED.md) with the type that could not be read). Two capitals that spell an element not in the formula are read again (13 CIFs corrected, 4 skipped; see below).

## The problem

`_atom_site_type_symbol` should hold a chemical element, optionally with a charge (`O2-`, `Fe3+`). 36 COD CIFs use other codes: labels copied into the type column (`HC13A`), water or hydroxyl markers (`Ow`, `OH-`), generic sites whose element is not fixed (`T` for a tetrahedral Si/Al site, `M` for a metal site), mixed sites (`ON`, `OW/Cl1`, `ST1` for Nb/Ti), and pseudo-atoms defined only in the authors' software (`SASH`, `SpHS`, `IS4`).

Reading the first letters as an element, as COD24 does, gives wrong elements: water oxygen `WO` becomes tungsten, the Nb/Ti site `ST1` sulfur, `IS4` iodine in a compound without iodine, `SASH` sulfur.

## Options considered

| Option | Effect | Decision |
|---|---|---|
| Resolve where certain, correct from the publication, skip the rest | 6 CIFs get the right element, 27 are corrected, 2 lose only dummy atoms; 1 CIF skipped and listed | **chosen**, implemented |
| Resolve where certain, leave the other atoms out | 6 CIFs fixed; 30 CIFs keep incomplete patterns | rejected: an incomplete pattern looks complete to a user |
| COD24's way: first letters | all 36 included, many with wrong elements | rejected: wrong elements give wrong patterns |
| Leave the atoms out (earlier behaviour) | 36 incomplete patterns, noted in the report | replaced |

## Rule

Take element E when the type starts with the symbol of a real element E that the CIF's own `_chemical_formula_sum` contains, and the rest is a lowercase suffix or an atom label (`Ow`, `Oc`, `Hy`, `Im1`, `HC13A`, `HN1B`, `OCl(1)`). Skip the CIF when, for any atom with coordinates, the element cannot be found this way:

- generic sites and pseudo-atoms: `T`, `M`, `X`, `Z`, `SASH`, `SpHS`, `IS4`; a rest in capitals without a number is a code, not a label;
- groups and mixed sites, where the rest is itself made of element symbols: `OH`, `OH-`, `CH`, `OH2` (group scatterers such as O+H, as the authors of COD 7201135 used for CH₂, CH₃ and water), `ON`, `OW/Cl1`;
- E not in the formula, or no formula (`ST1` in a compound without sulfur, `WO` in one without tungsten).

Code: `entry._element_of_type`, `entry._resolve_types`; tests `test_atom_types_that_are_not_elements`, `test_uncertain_atom_type_skips_the_cif`.

In CIFs without a type column the label is the type and is read the same way (`HC12` = H on C12); across COD, one such CIF lists a ring centroid as an atom: `Cen` in COD 2003591, which would read as cerium, missing from the CIF's formula, and is left out by a [correction](../corrections.md) (`dev/label_element_scan.py`). Atoms whose type is `?` or `.` are covered below.

## The article behind COD 7201135

[Ikeda, Kayamori and Mizukami, J. Mater. Chem. 19, 5518 (2009)](https://doi.org/10.1039/b905415d), p. 5524: the authors refined the CH₂ and CH₃ groups of the template molecules and the adsorbed water (or OH⁻) as virtual atoms, marked A–C in Table 4. In the CIF they appear as `CE` (CH₂; sites C1, C2, C5), `CH` (CH₃; C3, C4, C6), `WO` (water) and `OH`. The compound is Si₃₆O₆₈(OH)₈·2TEAOH·12.7H₂O, with neither cerium nor tungsten, yet gemmi reads `CE` as cerium and COD24 stores W12.7 and Ce8. hsrdb-tools builds it with these groups as their atoms, from the article ([corrections.md](../corrections.md)).

## All 36 CIFs

| COD | Types that are not elements | CIF formula | COD24 unit-cell formula | hsrdb-tools |
|---|---|---|---|---|
| [1544611](https://www.crystallography.net/cod/1544611.cif) | `M` | O2 Si | not in COD24 | corrected: stand-in atoms of the refinement ([corrections.md](../corrections.md)) |
| [1544612](https://www.crystallography.net/cod/1544612.cif) | `M` | O2 Si | not in COD24 | corrected: stand-in atoms of the refinement ([corrections.md](../corrections.md)) |
| [1544613](https://www.crystallography.net/cod/1544613.cif) | `M` | O2 Si | not in COD24 | corrected: stand-in atoms of the refinement ([corrections.md](../corrections.md)) |
| [1544625](https://www.crystallography.net/cod/1544625.cif) | `T`, `X`, `Z` | Al7 B3 Ca Fe3 H4 O31 Si5 | not in COD24 | corrected: cation assignment of the refinement ([corrections.md](../corrections.md)) |
| [1545186](https://www.crystallography.net/cod/1545186.cif) | `SASH` | Al56 H600 Na56 O684 Si136 | O505.51 Si136.01 Al55.99 Na55.59 S6.64 | corrected: spherical-shell scatterer ([corrections.md](../corrections.md)) |
| [1545415](https://www.crystallography.net/cod/1545415.cif) | `SASH` | Al56 Na56 O384 Si136 | O401.89 Si136.01 Al55.99 Na47.97 S0.68 | corrected: spherical-shell scatterer ([corrections.md](../corrections.md)) |
| [1545416](https://www.crystallography.net/cod/1545416.cif) | `SASH` | Al56 H600 Na56 O684 Si136 | O556.12 Si136.01 Al55.99 Na56.95 S4.40 | corrected: spherical-shell scatterer ([corrections.md](../corrections.md)) |
| [1545417](https://www.crystallography.net/cod/1545417.cif) | `SASH` | Al56 H300 Na56 O534 Si136 | O460.36 Si136.01 Al55.99 Na56.34 S2.01 | corrected: spherical-shell scatterer ([corrections.md](../corrections.md)) |
| [1546291](https://www.crystallography.net/cod/1546291.cif) | `T` | Al0.024 Na0.015 O4 Si1.978 | T68.00 O136.00 C18.88 | corrected: Si/Al from the article's composition ([corrections.md](../corrections.md)) |
| [1549188](https://www.crystallography.net/cod/1549188.cif) | `SpHS` | C62 H4 Br2 Cl2 | C3.00 H4.00 Cl2.00 Br2.00 | corrected: spherical-shell C₆₀ ([corrections.md](../corrections.md)) |
| [1557919](https://www.crystallography.net/cod/1557919.cif) | `Ow` | C14 H33 Cl3 N4 O2 Zn | Zn2.00 Cl6.00 N8.00 O4.00 C28.00 | resolve |
| [1561471](https://www.crystallography.net/cod/1561471.cif) | `Im1`, `Om1` | Cu H2 I3 O2 Pb2 | Pb64.00 I96.00 Cu32.00 O64.00 | resolve |
| [1563539](https://www.crystallography.net/cod/1563539.cif) | `OH-` | Bi2.15 Ca2 H2.5 O12.475 Si3 | Ca4.00 Bi4.30 Si6.00 O24.95 | corrected: hydroxyl oxygen, H not refined ([corrections.md](../corrections.md)) |
| [1569813](https://www.crystallography.net/cod/1569813.cif) | `OiI` | Al8.08 B3.51 Ca0.25 F0.07 H3.36 Li0.76 Mn0.16 Na0.42 Oi31 Si5.49 | Na1.37 Ca1.63 Li6.98 Mn2.02 Al18.00 B9.00 Si18.00 O91.19 F1.76 H9.00 | corrected: O²⁻ scattering-factor label ([corrections.md](../corrections.md)) |
| [1569815](https://www.crystallography.net/cod/1569815.cif) | `OiI` | Al24.24 B10.53 Ca0.74 F0.2 H10.08 Li2.27 Mn0.49 Na1.26 Oi93 Si16.47 | Na0.94 Ca2.06 Li7.07 Mn1.93 Al18.00 B9.00 Si18.00 O91.24 F1.76 H9.00 | corrected: O²⁻ scattering-factor label ([corrections.md](../corrections.md)) |
| [1569816](https://www.crystallography.net/cod/1569816.cif) | `OiI` | Al8.08 B3.51 Ca0.25 F0.07 H3.36 Li0.76 Mn0.16 Na0.42 Oi31 Si5.49 | Na1.50 Ca1.50 Li5.93 Mn3.07 Al18.00 B9.00 Si18.00 O91.25 F1.76 H9.00 | corrected: O²⁻ scattering-factor label ([corrections.md](../corrections.md)) |
| [2004310](https://www.crystallography.net/cod/2004310.cif) | `OF1`, `OF2`, `OF3`, `ST1`, `ST2` | F4 K3 Nb3 O11 Ti2 | S10.00 K6.00 F2.00 O28.00 | corrected: Nb/Ti and O/F mixed sites of the refinement ([corrections.md](../corrections.md)) |
| [2004555](https://www.crystallography.net/cod/2004555.cif) | `ON` | C11 H21 Cd N7 Ni O | Cd4.00 Ni4.00 O8.00 N24.00 C44.00 H84.00 | corrected: N/O site of a ligand on a twofold axis ([corrections.md](../corrections.md)) |
| [2004556](https://www.crystallography.net/cod/2004556.cif) | `ON` | C12 H25 Cd N7 Ni O2 | Cd4.00 Ni4.00 O12.00 N24.00 C48.00 H100.00 | corrected: N/O site of a ligand on a twofold axis ([corrections.md](../corrections.md)) |
| [2004740](https://www.crystallography.net/cod/2004740.cif) | `OCl(1)` | C14 H42 Cl2 O7 Pd S6 | Pd4.00 S24.00 C56.00 Cl12.00 O24.00 H176.00 | resolve |
| [2006347](https://www.crystallography.net/cod/2006347.cif) | `HC13A`, `HC13B`, `HC14A`, `HC14B`, `HC16A`, `HC16B`, `HC19A`, `HC19B`, `HC20A`, `HC20B`, `HC21A`, `HC21B`, `HC22A`, `HC22B`, `HC23A`, `HC23B`, `HC4A`, `HC4B`, `HN1A`, `HN1B` | C28.5 H38 N2 O5.5 S | S8.00 O22.00 N8.00 C110.00 H148.00 | resolve |
| [2007448](https://www.crystallography.net/cod/2007448.cif) | `Oc` | C14 H42 Cl2 O7 Pt S6 | Pt4.00 S24.00 C56.00 Cl12.00 O24.00 H176.00 | resolve |
| [2010123](https://www.crystallography.net/cod/2010123.cif) | `M1`, `M2`, `M3`, `M4` | C16 H24 Cl2 Rh2 | Rh8.00 Cl8.00 C64.00 H96.00 | left out: double-bond midpoints flagged `dum` |
| [2017919](https://www.crystallography.net/cod/2017919.cif) | `IS4`, `IS5`, `IS9` | C23 H26 F2 N2 O4 | C90.05 O15.57 H104.00 F7.82 N7.88 I105.55 | corrected: interatomic scatterers left out ([corrections.md](../corrections.md)) |
| [2019464](https://www.crystallography.net/cod/2019464.cif) | `M` | Cu1.04 Ni5.2 Sn8.7 Zn4.16 | not in COD24 | corrected: Ni/Cu/Zn sites with the analysed composition ([corrections.md](../corrections.md)) |
| [2100388](https://www.crystallography.net/cod/2100388.cif) | `M` | Hf0.006 O1.968 Zr0.994 | Zr2.00 O3.94 | corrected: the site label Zr0.994Hf0.006 and the formula ([corrections.md](../corrections.md)) |
| [2100389](https://www.crystallography.net/cod/2100389.cif) | `M` | (none) | Zr2.00 O3.94 | corrected: the site label Zr ([corrections.md](../corrections.md)) |
| [2209646](https://www.crystallography.net/cod/2209646.cif) | `X` | Cl X O4 | not in COD24 | skip: the corrigendum makes the cation ammonium, without a model ([corrections.md](../corrections.md)) |
| [4060015](https://www.crystallography.net/cod/4060015.cif) | `X` | C40 H58 Ce F5 | Ce2.00 F10.00 C80.00 H116.00 | left out: ring centroids flagged `dum` |
| [4518778](https://www.crystallography.net/cod/4518778.cif) | `M` | M4 Na2.97 Sb | not in COD24 | corrected: S/Cl anion site of the refinement ([corrections.md](../corrections.md)) |
| [5910139](https://www.crystallography.net/cod/5910139.cif) | `Cb` | Cb K O3 | K1.00 C1.00 O3.00 | corrected: columbium = niobium ([corrections.md](../corrections.md)) |
| [7012100](https://www.crystallography.net/cod/7012100.cif) | `OW/Cl1`, `OW/Cl2`, `OW/Cl3` | C26.5 H33 Cl2 N4 Ni O5 | Ni2.00 Cl9.00 O5.00 N8.00 C53.00 H52.00 | corrected: water/chloride sites with the compositions the CIF's own atom-type descriptions give (0.90 O + 0.10 Cl etc.; [corrections.md](../corrections.md)) |
| [7034643](https://www.crystallography.net/cod/7034643.cif) | `Hy` | C70 H122 N6 O Si4 Zn4 | Zn16.00 H488.00 N24.00 C280.00 Si16.00 O4.00 | resolve |
| [7201135](https://www.crystallography.net/cod/7201135.cif) | `CH`, `OH`, `WO` | (none) | Si36.00 O76.00 N2.00 Ce8.00 C8.00 W12.70 | corrected from the article |
| [7201136](https://www.crystallography.net/cod/7201136.cif) | `CH`, `OH` | (none) | Si36.00 O76.00 N3.25 C19.51 H3.25 | corrected from the article |
| [7232188](https://www.crystallography.net/cod/7232188.cif) | `SpHS` | C62 H2 Br2 Cl4 | C3.00 H2.00 Cl4.00 Br2.00 | corrected: spherical-shell C₆₀ ([corrections.md](../corrections.md)) |

## Two capitals read as the wrong element

Some types and labels are two capital letters that happen to spell an element: `OS1`, `PB`, `CO3`. gemmi, like the official database, reads them as that element (osmium, lead, cobalt). When the CIF's formula does not contain that element, hsrdb-tools reads the capitals again:

- **CIF without a type column** (element taken from the label): the label is element + site, in the style of the American Mineralogist database (`PB` = P on site B, `OS1` = O on site S1, `HO2a` = H on O2). The first letter's element is taken if the formula contains it, unless the rest is another element of the formula (`CN`: skipped).
- **Type column**: the capitals are a code or a group and are read like the other types above; if the element is not certain, the CIF is skipped.
- A label in normal case that confirms the element (`Co1`, `Mg`) keeps it: the formula leaves out a dopant.

Result across COD (`dev/element_formula_scan.py`):

| COD | Capitals | Read by gemmi and COD24 as | Now |
|---|---|---|---|
| 9004436, 9005290, 9007680, 9012653, 9014495, 9015706, 9016053 | `OS1`, `OS2`, … | osmium | O (oxygen of a sulfate or on site S) |
| [7125504](https://www.crystallography.net/cod/7125504.cif) | `HO2a`, … | holmium | H |
| [9016317](https://www.crystallography.net/cod/9016317.cif) | `PB`, `SB` | lead, antimony | P, S (fluorapatite) |
| [9013840](https://www.crystallography.net/cod/9013840.cif) | `YB` | ytterbium | Y (burbankite) |
| [9004480](https://www.crystallography.net/cod/9004480.cif) | `PT` | platinum | P |
| [9004790](https://www.crystallography.net/cod/9004790.cif) | `CO3` | cobalt | C (carbonate carbon) |
| [2004294](https://www.crystallography.net/cod/2004294.cif) | `CS` | caesium | C |
| [2101903](https://www.crystallography.net/cod/2101903.cif) | `TM(1)`–`TM(14)` | thulium | corrected: transition-metal sites, Fe/Cr in the alloy's ratio ([corrections.md](../corrections.md)) |
| [7027327](https://www.crystallography.net/cod/7027327.cif), [7012620](https://www.crystallography.net/cod/7012620.cif) | `RE1`, `CO`/`SI` | rhenium, cobalt/silicon | corrected: Re, and Co and Si, as the articles give them; the CIFs' formulas lack Re or split Co and Si into C, O, S and I ([corrections.md](../corrections.md)) |
| [5000363](https://www.crystallography.net/cod/5000363.cif) | `TL` | thallium | Tl: the CIF's atom_type loop describes type `TL` as Tl (a 3% substitution on the Zn site that the formula leaves out) |

Where a type would otherwise be skipped, the CIF's own `_atom_type_description` decides when it is exactly an element symbol (type `TL` described as `Tl`, COD 5000363). Code: `entry._capitals_misread`, `entry._described_elements`; tests `test_two_capitals_read_against_the_formula`, `test_atom_type_description_names_the_element`.

## Atoms with no type at all: COD 9002229 and 9002230

Atoms without a type symbol get the element gemmi reads from their label. When that fails, the atom is left out only if it is a dummy (flagged `dum` or named `Dummy`/`DUM`, e.g. ring centroids) or has zero occupancy; otherwise the CIF is skipped (`atom 'I1' has no element`). A dummy whose type (or label, without a type column) names no element is likewise left out: the double-bond midpoints `M1`–`M4` of COD 2010123 and the ring centroids of type `X` in 4060015. Across COD this concerns two CIFs, which are instead built from their article ([corrections.md](../corrections.md)):

[9002229](https://www.crystallography.net/cod/9002229.cif) and [9002230](https://www.crystallography.net/cod/9002230.cif) (AMCSD) are the triclinic and monoclinic "modified" chlorites of [Guggenheim and Zhan, Am. Mineral. 84, 1415 (1999)](https://www.degruyterbrill.com/document/doi/10.2138/am-1999-0920/html). Their interlayer sites `I1`–`I6` (monoclinic: `X1`–`X4`) have type `?`, and the formula counts them as iodine (`I12`) or `X12`. According to the article (pp. 1416–1418), these sites cannot be identified as cation or oxygen sites; each is on average about half occupied by Mg or O, and they were refined as split atoms. The CIFs give these sites occupancy 1 and omit the 15% intergrowth tetrahedral atoms the article refined, so the published model cannot be rebuilt from the CIFs alone; the correction rebuilds the interlayer from the article's composition.
