# hsrdb-tools

Build reference databases for Malvern Panalytical HighScore (`.hsrdb`) from [Crystallography Open Database](https://www.crystallography.net) (COD) CIF files, for HighScore 3.x and for HighScore 4.x and later.

Official COD databases for HighScore appear every year or two; the latest is from May 2024. COD itself grows every week. With hsrdb-tools you build a current database yourself, from all of COD or from any set of CIF files.

## What goes into a database

Each CIF becomes one reference pattern with the same content as the official COD databases:

- **Stick pattern** for search-match: d-spacing, intensity and hkl of up to 203 lines, plus I/Ic (reference intensity ratio against corundum)
- **Search filters**: cell, space group, crystal system, formula, element list and subfile (organic, inorganic, mineral, metallic)
- **Names and literature**: systematic, mineral and common name; authors, title, journal, volume, year and pages
- **Crystal structure** (atoms, occupancies, displacement parameters), which HighScore Plus uses for Rietveld refinement

Structures are stored in the space group's reference setting (standard axes, origin choice 2, hexagonal axes for rhombohedral groups), the setting HighScore assumes for a Hermann–Mauguin symbol. CIFs in other settings (unusual axes or centring, a shifted origin, supercells) are converted from their listed symmetry operations or Hall symbol with spglib, and only when spglib finds exactly the listed group. Within the reference setting, the cell is the one the official COD databases choose: triclinic cells Niggli-reduced; monoclinic cells with β ≥ 90° and the shortest a + c among those that keep the space-group symbol (e.g. P2₁/n and P2₁/a become P2₁/c with c′ = a + c, not β ≈ 140°); orthorhombic axes permuted, where the symbol allows it, to the shortest a, then b.

## Calculation model

Each line intensity is multiplicity × |F|² × Lorentz-polarisation, for the structure exactly as published:

| | Choice |
|---|---|
| Radiation | Cu Kα1, 1.540598 Å; lines up to 90° 2θ (140° if that gives fewer than 10 lines), at most 203 |
| Form factors | IT92 (International Tables Vol. C, Table 6.1.1.4), neutral atoms; charges in type symbols are ignored |
| Anomalous scattering | f′, f″ by Cromer–Liberman at 8047.8 eV (gemmi: Brennan–Cowan with the Kissel–Pratt correction). f″ agrees with other tables; f′ of heavy atoms differs between tables by up to about 0.9 electrons (e.g. Pb: −3.95 here, −4.31 Chantler/NIST, −4.82 International Tables), up to about 2% of \|F\|² for reflections dominated by atoms heavier than Ba |
| Friedel pairs | \|F(h)\|² and \|F(−h)\|² averaged: both fall on the same powder line |
| Displacement | published U or B, isotropic or anisotropic (U, B or β tensors), rotated with each symmetry image; B = 0.5 Å² for atoms without one (HighScore's default, also used by the official databases); physically invalid tensors are replaced by the isotropic value (the tensor's own equivalent value, trace / 3, if the CIF gives none), a negative isotropic value by 0; anharmonic (Gram–Charlier) terms are not used, as in COD24 |
| Occupancy | as published (chemical occupancy) |
| Atoms near special positions | symmetry images closer than 0.01 Å are one atom, placed exactly on the special position; images up to 0.4 Å apart are one atom when the CIF says so (site-symmetry order or multiplicity) or when they could not all be occupied (occupancy × images above 1.01, a margin for rounded occupancies); otherwise they are the separate positions of a split, disordered atom. Images of an atom on a special position that the CIF states as general count separately only if merging them contradicts the CIF's formula (or density) and counting them does not |
| Dummy atoms | atoms flagged `dum` (e.g. ring centroids) do not scatter; they stay in the structure record, flagged |
| Modulated structures | from the CIF's superspace description (msCIF, as JANA writes it): main and satellite reflections, each atom's scattering averaged over its internal coordinate with its displacement, occupancy and displacement-parameter functions; satellite reflections weaker than 0.5 on the 0–1000 scale are left out, and a satellite line carries hkl 0 0 0. The structure record holds the basic structure. Composites and functions specific to JANA are stored as the basic structure, as in COD24; see [Modulated structures](docs/differences/modulated.md) |
| Non-atomic scatterers | where an article refined matter as a spherical shell (a water layer in zeolite cages), a [correction](docs/corrections.md) adds its atoms spread over the sphere, form factor n·f(q)·sin(qR)/(qR) as in FullProf; it scatters and counts in formula and density, but the structure record, which holds atoms only, leaves it out |
| Lorentz-polarisation | (1 + cos²2θ)/(sin²θ cos θ): Bragg–Brentano geometry, no monochromator |
| Lines | a reflection within 0.005° 2θ of the previous one joins its line; reflections of zero intensity are left out |
| I/Ic | strongest line against corundum (COD 1000032), 50:50 by weight, same model |
| Density | IUPAC atomic weights (gemmi) |

## Install

```
pip install .
```

Requires Python 3.10+, [gemmi](https://gemmi.readthedocs.io), numpy and [spglib](https://spglib.readthedocs.io). A `Dockerfile` with all dependencies is included.

## Build a database

Mirror COD. The first run downloads about 110 GB (535,000 files); later runs fetch only changes:

```
rsync -a --delete rsync://www.crystallography.net/cif/ cod/
```

Build databases for both HighScore generations in one pass:

```
hsrdb-tools add cod/ --hs4 COD_HS4x.hsrdb --hs3 COD_HS3x.hsrdb -j 8
```

| Option | Meaning |
|---|---|
| `--hs4 PATH` | database for HighScore 4.x and later |
| `--hs3 PATH` | database for HighScore 3.x |
| `-j N` | worker processes (default: all cores) |
| `--base DB` | copy an existing database and add only CIFs it does not contain (one output only) |
| `--append` | add to the given databases in place |

All of COD takes about 6 hours on a 4-core ARM64 Ampere Altra machine and gives two databases of about 7 GB each. Other CIF files work too if each holds one data block and carries a 7-digit COD-style id (a numeric data block name, as in COD, else `_cod_database_code`; COD's own 1010542 carries another entry's `_cod_database_code`).

A build writes `<name>.partial` and renames it only on success, so an interrupted build leaves no half-written database (with `--append`, nothing is committed). CIFs that cannot be converted, for example without coordinates, unit cell or a recognisable space group, are skipped and listed with the reason in a report next to the first output (`--hs4` if given): `COD_HS4x.skipped.tsv` for `COD_HS4x.hsrdb`. Unexpected errors are printed with a traceback and make the command exit with status 1; the databases are still written, without those CIFs.

### In HighScore

Customize → Manage Databases → Add HighScore Database, select the file, then tick **Use**.

**Use it instead of the official COD database, not next to it.** Both use COD's `96-xxx-xxxx` reference codes, and HighScore refuses two active databases that declare the same code prefixes.

## Compatibility

| HighScore | Status |
|---|---|
| 4.x and later (`--hs4`) | opened and used for search-match in HighScore 4.7 |
| 3.x (`--hs3`) | patterns opened in HighScore 3.x; structure records follow the official 2014 format but have not yet been opened in 3.x |
| Plus (Rietveld, Convert Pattern to Phase) | structure records not yet tested |

Reports from HighScore Plus and 3.x users are welcome.

## How close is it to the official databases?

A full build from COD (September 2026) gives 533,455 entries: every entry of the official COD24 database and 21,872 newer ones. 1,669 CIFs are skipped, mostly because they contain no atom coordinates (1,451) or no usable symmetry (170); none of them is in COD24.

Every entry of both databases passes `dev/audit_database.py` (valid peak lists, consistent metadata, element filters, and structure records that decode and re-encode identically). All 511,583 entries of the official COD24 database were compared:

| | Agreement with COD24 |
|---|---|
| Database schema | identical for each HighScore version |
| Structure record format | byte-identical re-encoding of official records (2,000 per version checked) |
| Space group | 99.47% |
| Cell volume within 0.5% | 99.99% |
| Element set | 99.56% |
| I/Ic within 10% | 97.17% |
| All strong official lines present | 99.98% |
| Largest strong-line intensity difference | median 12, 90th percentile 111 (scale 0–1000) |

7% of entries have a strong line whose intensity differs by more than 150/1000. The main causes are the first two differences below.

## Independent verification

Every pattern is checked against [cctbx](https://cctbx.github.io), an independent crystallographic library, with `dev/verify_cctbx.py` (in the `dev` stage of the Dockerfile) over all of COD:

- **Same structure**: cctbx calculates |F|² for the reflections of each stored pattern from the structure as hsrdb-tools reads it (same f′, f″). This checks structure factors, form factors, displacement factors, symmetry expansion, multiplicities and the Lorentz-polarisation factor.
- **Own reading of the CIF**: cctbx reads the raw CIF with its own parser and symmetry handling, in the CIF's own setting, and calculates the lines. This checks parsing, settings and site handling.

`dev/verify_explain.py` assigns every difference a cause or lists it as unexplained. The full run over all 533,275 CIFs of the build then current gave no calculation differences and no unexplained reading differences; for this release, the 31,317 CIFs whose entries changed since (corrections, modulated structures, recovered CIFs, and a sample of the cells now stored as COD24 stores them) were verified again with the final code, and against the finished database, with the same result (the corrected and modulated entries differ from cctbx's reading of the raw CIF by design). The known causes are cases where cctbx reads a CIF differently by design or by limitation (for example, it drops a whole displacement column when one value is `?`, and reads dummy atoms as scatterers); each one is listed in the script.

## Differences from the official COD databases

Where hsrdb-tools deliberately does something different from the official COD24 database, and why. Each point was checked against COD24 entries; the linked pages give the evidence and the options that were considered. [docs/decisions.md](docs/decisions.md) lists every decision in the calculation, including those where hsrdb-tools matches COD24.

| | COD24 | hsrdb-tools | Why |
|---|---|---|---|
| [Coincident reflections](docs/differences/coincident-reflections.md) | reflections that are not symmetry-equivalent but fall at the same angle (Laue classes −3, 6/m, 4/m, m−3, pseudo-symmetric cells) appear to count once; e.g. COD 9009599 lines are up to 40% low | all are added | a measured line contains all of them |
| [Atom types that are not elements](docs/differences/atom-types.md) (36 CIFs, and 17 with two capitals) | guessed from the first letters: water oxygen `WO` becomes tungsten, a mixed Nb/Ti site `ST1` sulfur, `IS4` iodine in an iodine-free compound, sulfate oxygen `OS1` osmium, `PB` (P on site B) lead | taken where certain (6 CIFs), built from the publication or the CIF's own statements (27), dummies left out (2), skipped where neither decides (1); see [Atom types](#atom-types-that-are-not-elements) | a wrong element gives a wrong pattern |
| [Cyanide pseudo-atom `CN`](docs/differences/cyanide.md) (4 CIFs) | read as carbon; the nitrogen is lost, and three keep a wrong origin (16 metal atoms per cell instead of 8) | C and N at the site, origin corrected | the compilation's model of the cyanide group |
| [Dummy atoms](docs/differences/dummy-atoms.md) (`calc_flag dum`: ring centroids, atoms with unknown positions; 438 CIFs) | counted as atoms (extra H, C or N in the formula) | not scattering and not in formula, density or element filters; kept, flagged, in the structure record | a dummy atom is a geometric point, not an atom |
| [Physically invalid anisotropic tensors](docs/differences/invalid-adp.md) (about 2,800 CIFs) | used as given | replaced by the isotropic value | a tensor with a negative mean-square displacement gives unphysical intensities |
| [Modulated structures](docs/differences/modulated.md) (198 CIFs) | the basic structure: main lines only, as if every atom sat at its average position | main and satellite lines from the CIF's superspace model; the structure factors agree with those the refinement program lists in the CIF (30 CIFs, R at most 0.9%) | the modulation changes a main line by 10/1000 or more in 147 of them and gives satellites of 10/1000 or more in 138, which a measured pattern shows |
| [Corrections from publications](docs/corrections.md) (57 CIFs) | CIF as deposited: interlayer sites as iodine (COD 9002229), virtual atoms as cerium and tungsten (7201135) | the article's model, with the CIF file unchanged and the entry marked | the article states what the CIF leaves out |
| [Non-standard settings](docs/differences/non-standard-settings.md) (81 CIFs, e.g. `B 1 21/d 1`, supercells) | 29 of them included, 5 with a wrong space group | all standardised from the CIF's own symmetry operations with spglib | see [Skipped CIFs](#skipped-cifs) |
| COD content | May 2024 | the mirror you build from; about 30% of entries with large differences are CIFs COD revised after May 2024 | newer data |

The same as COD24, checked: cells are constrained to the crystal system (exact angles, a = b or a = b = c where required) when a CIF's cell does not quite fit its space group; atoms without displacement parameters get B = 0.5 Å²; anisotropic B and U tensors, split atoms next to special positions, radiation, Lorentz-polarisation, the line cap, the reference setting and its cell choice (checked on a sample of 16,600 entries: 99.97% of monoclinic, 99.9% of orthorhombic and 98.2% of triclinic cells as COD24; the rest are near-ties or COD24 keeping a CIF's axes under the reference symbol) and the subfiles (inorganic, organic, mineral, metallic).

Other known differences:

- Atom coordinates can differ from the official records by a symmetry-equivalent choice; both describe the same structure and give the same pattern.
- Wyckoff letters stay blank for a small number of disordered sites where spglib does not reproduce the published space group.

## Skipped CIFs

A CIF is skipped when no correct pattern can be calculated from it. Every skipped CIF is listed with its reason in the build report (`<name>.skipped.tsv`); [docs/SKIPPED.md](docs/SKIPPED.md) lists those of the COD build with a download link each (`dev/skipped_md.py` writes it from the report). For COD (September 2026):

| Reason | CIFs | Details |
|---|---|---|
| No atom coordinates | 1,451 | the CIF has a cell and symmetry but no atom positions (structure not determined, or coordinates only in the paper). Five CIFs of this kind whose official COD24 entries were built from unusable data (2300247, 2300248, 2300253, 2300257, 5900030) are included as COD24 has them, marked as not correct structures ([corrections](docs/corrections.md)) |
| No symmetry information | 155 | neither symmetry operations nor a Hall or Hermann–Mauguin symbol; nearly all are empty entries from papers of 1926–1962, the rest modulated or composite structures whose superspace model the CIF lacks |
| Invalid or ambiguous symbol, no operations | 10 | e.g. `Pbc2`, `?P?`, `unknown`, `P21 or P21/m`: cell-only entries (three are COD duplicates of built entries) and modulated structures without their superspace model, each checked against its article ([corrections](docs/corrections.md)) |
| Superspace symbol without operations | 4 | `X4bm`, `Cmca(00γ)s00`: modulated structures whose superspace operations the CIF does not list |
| Symmetry operations that do not fit the structure | 1 | COD 2100427 |
| No valid unit cell | 46 | 45 are empty entries that give only a formula and a reference (no atoms, cell or symmetry); one (2104629) is a 1,440-atom supercell model without its cell |
| Atom type that is not an element and cannot be resolved | 2 | 2209646 (the corrigendum makes its cation ammonium but gives no model) and 2101649 (metal sites `M1`–`M20` of an electron-microscopy model, oxygen not published) |

Skipped CIFs are checked against their original articles; those that an article can complete are built with a [correction](docs/corrections.md), and those checked but not correctable are listed there with the reason.

Some atoms are left out of a pattern without skipping the CIF, and noted in the same report: atom_site rows without coordinates or without any element (e.g. ring centroids `CNT1` whose type is `.`), and dummy atoms, which stay in the structure record, flagged.

### Atom types that are not elements

36 COD CIFs give some atoms a type symbol that is not a chemical element: water or label codes (`Ow`, `HC13A`), group scatterers (`OH`, `CH`), generic sites (`T`, `M`), mixed sites (`ON`, `OW/Cl1`) and undefined pseudo-atoms (`SASH`, `SpHS`, `IS4`). hsrdb-tools takes the element when the type starts with an element that the CIF's own formula contains and the rest is only a suffix or an atom label (6 CIFs), and otherwise builds the CIF from its publication or its own statements with a [correction](docs/corrections.md) (27 CIFs, e.g. COD 7201135, 7201136, the melanophlogites and faujasites), leaves out two sets of dummy atoms and skips one CIF (2209646) that neither decides; corrections also complete two CIFs whose interlayer atoms have no type at all (9002229, 9002230). Types or labels of two capitals that spell an element missing from the CIF's formula are read again: sulfate oxygens `OS1` read by gemmi and COD24 as osmium, `PB` (P on site B) as lead, `HO2a` as holmium; 13 CIFs are corrected this way and 4 skipped. The official COD24 database guesses from the first letters instead (water `WO` read as tungsten, a mixed Nb/Ti site `ST1` as sulfur). Details and the full list: [docs/differences/atom-types.md](docs/differences/atom-types.md).

### Cyanide pseudo-atom `CN`

Four CIFs from Wyckoff's 1931 compilation give cyanide as one pseudo-atom `CN`. hsrdb-tools puts C and N at that site, the compilation's own model of the cyanide group as a single spherical scatterer, and for the three K₂M(CN)₄ reads the coordinates in origin choice 1, which the file mislabels as origin choice 2. COD24 reads `CN` as carbon only. Details: [docs/differences/cyanide.md](docs/differences/cyanide.md).

## Format

[docs/FORMAT.md](docs/FORMAT.md) documents the `.hsrdb` layout, including the binary structure record in both versions.

## Development

```
pip install -e '.[test]'
pytest
```

The tests cover the formats, space-group settings, pattern calculation, and regressions for problems found in full COD runs.

`dev/` holds the scripts used to decode the format and to check output:

- `audit_database.py` checks every entry of a built database
- `compare_databases.py` compares a built database with an official one
- the other scripts compare individual fields; most expect the official `COD24_HS4x.hsrdb` in the working directory

## Data, credits and trademarks

COD data are open access; see [crystallography.net](https://www.crystallography.net) for the terms, including attribution of the original publications. Pattern and symmetry calculations use [gemmi](https://gemmi.readthedocs.io) and [spglib](https://spglib.readthedocs.io).

HighScore is a trademark of Malvern Panalytical. This project is not affiliated with or endorsed by Malvern Panalytical.

## Licence

MIT, see [LICENSE](LICENSE).
