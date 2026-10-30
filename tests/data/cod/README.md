# Real COD CIFs for the regression test

`tests/test_cod_regression.py` runs `entry.from_cif` and the structure-record encoders on each CIF here and compares the result with `expected.json`. `python dev/update_expected.py` regenerates `expected.json` from the current code and prints which ids changed and which fields moved.

## Source and terms

The CIFs are unmodified copies of files of the Crystallography Open Database, taken from a local mirror of https://www.crystallography.net/cod/ (file `https://www.crystallography.net/cod/<id>.cif`).

Terms, from the COD home page https://www.crystallography.net/cod/ (checked 2026-10-24):

> All data in the COD and the database itself are dedicated to the public domain and licensed under the CC0 License. Users of the data should acknowledge the original authors of the structural data.

CC0 text: https://creativecommons.org/publicdomain/zero/1.0/. The headers of most CIFs here repeat it ("All data on this site have been placed in the public domain by the contributors"); those from IUCr journals and from the American Mineralogist Crystal Structure Database say instead that the file "may be used within the scientific community so long as proper attribution is given to the journal article from which the data were obtained", and 3000357 has no header. Redistribution is allowed under the database's CC0 dedication; the acknowledgement is met by keeping the files unchanged, with their `_publ_author_name`, `_journal_*` and `_cod_database_code` items. How to cite COD: https://wiki.crystallography.net/citingcod/.

## Files

| COD id | Why it is here |
|---|---|
| 1000032 | corundum: plain inorganic, rhombohedral in hexagonal axes, mineral name and mineral subfile |
| 1000041 | NaCl: plain cubic inorganic |
| 1008063 | anisotropic displacement parameters given as B_ij |
| 4319054 | anisotropic displacement parameters given as beta_ij |
| 1524007 | Fe-Cr alloy: metallic subfile |
| 1524039 | triclinic P-1, Niggli reduction (beta 89.76 becomes 90.24) |
| 9001314 | triclinic P-1 mineral (werdingite), Niggli reduction; split O19A/O19B atoms, O19B near an inversion centre |
| 2010123 | organic (C-H bonds: organic subfile), P 1 21/n 1 stored as P 1 21/c 1 (beta 91.79 becomes 127.82); 4 atom_site rows without coordinates dropped |
| 1537595 | setting gemmi does not know: `P 1 21/c 1 (2*a+c,b,c)` supercell symbol, converted via spglib |
| 1524367 | orthorhombic axis permutation: `P n a m` stored as `P n m a` |
| 2003598 | symbol `P m n b` (axis permutation), and the listed operations name another group: the symbol is used (symmetry_note) |
| 2100788 | symmetry stated twice, disagreeing: H-M symbol `A 21 a m` (C m c 21) against Hall symbol `p -2 2a` (P m c 21), no listed operations; symbol used |
| 1010022 | rhombohedral, reverse setting `R -3 m HR`, stored as R -3 m in hexagonal axes |
| 1010918 | bassanite: `P 32 2 1 S`, operations with the origin shifted by c/3 |
| 2014900 | listed operations in a shifted-origin setting that gemmi does not match to the symbol |
| 1524392 | partially occupied Ca and Tl sites with images 0.2-0.5 A apart (split-atom decision) |
| 2101363 | hydrogen site whose symmetry images lie 0.09 A apart (near-special-position merge) |
| 2002730 | sites on special positions that the CIF states as general: images counted separately (merge_note) |
| 5910160 | Wyckoff compilation (SiP2O7, P a -3): every atom of the cell listed besides the space group, no Z; rows that are symmetry images of another row counted once (merged_images), the formula's proportions decide |
| 5910217 | Wyckoff compilation (Mo2BC, C m c m): as 5910160, but counted separately the content is still proportional to the formula; the impossible density (35 g/cm3) decides |
| 5910088 | Wyckoff compilation (scheelite CaWO4): every atom of the cell listed, coordinates in origin choice 1; correction `spacegroup` I 41/a:1, then the duplicate rows count once |
| 5910132 | Wyckoff compilation (KReO4): as 5910088, but the doubled content stays proportional to the formula; correction omits the duplicate rows |
| 5910161 | Wyckoff compilation (spinel NiMn2O4): no setting gives one cell (a stray row); whole structure replaced by a corrected one of the same publication (not checked against the article) |
| 9000012 | Sb2O3 (valentinite): a and b swapped against the coordinates (12 impossible contacts); whole structure replaced by a corrected one of the same publication (not checked against the article) |
| 1000408 | TlAlF4: 12 impossible contacts; whole structure replaced by a corrected one of the same publication (not checked against the article) |
| 1010542 | baryte: dummy atom (`dum`) that does not scatter; mineral name |
| 1101168 | physically invalid ADPs replaced (invalid-adp.md) |
| 1544611 | melanophlogite: correction with `sites` / `stand_in` species for the cage atoms; atom types that are not framework elements (atom-types.md) |
| 5910139 | correction `types`: element Cb (columbium) read as niobium (atom-types.md) |
| 7201135 | correction `types`: virtual atoms of a Rietveld refinement mapped to their atoms (atom-types.md) |
| 5910129 | cyanide group given as one atom CN (cyanide.md), correction |
| 6000055 | correction with a `cell` and `add` (the CIF gives no atoms) |
| 6000615 | correction with `add` and `species` (the CIF gives no atoms) |
| 6000092 | correction with a `superspace` key: modulated structure built from the article |
| 1549598 | mullite: incommensurate (3+1)-dimensional modulation, Fourier waves |
| 4302713 | incommensurate (3+1)-dimensional modulation with a crenel occupancy; metallic subfile |
| 4001801 | incommensurate (3+2)-dimensional modulation (slowest file, about 7 s) |
| 2310761 | commensurate modulation, `_jana_cell_commen_t_section_1` 0.125; crenel occupancy |
| 2311738 | commensurate modulation, t section 0.05 |
| 4331459 | commensurate modulation, t section 0.167 |
| 2310038 | commensurate modulation, t section -0.25; crenel occupancy; sites stated general on special positions |
| 4002590 | commensurate modulation, t section 0; sawtooth function (Na and Mn on one site) |
| 2101649 | skipped: atom type `M1` is not an element |
| 3000357 | skipped: no atom coordinates |
