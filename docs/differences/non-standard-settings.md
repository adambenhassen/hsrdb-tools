# Non-standard space-group settings

**COD24:** 29 of these 81 CIFs are included; 5 of them with a wrong space group.
**hsrdb-tools:** all 81 included, in the standard setting.

## Method

gemmi knows the space-group settings tabulated in International Tables. These CIFs use others: centred cells such as `B 1 21/d 1` or `F -4 d 2`, supercells written as `P 1 21/c 1 (2*c,2*a+c,b)`, origins shifted by a third, the reverse rhombohedral setting, or symbols with trailing text such as `P 32 2 1 S`.

hsrdb-tools expands the atoms with the CIF's own symmetry operations (or its Hall symbol, which can carry a change of basis) and identifies the space group of the result with spglib. The structure is accepted only if the operations form a group, preserve the cell's metric, and spglib finds exactly the listed operations; it is then transformed to the reference setting. `dev/recover_check.py` confirms for every CIF below that the standardised pattern equals the pattern of the CIF expanded in its own cell. The tests cover one example of each notation.

Rejected: [5000046](https://www.crystallography.net/cod/5000046.cif) (`C 1 1 21/d`), whose operations do not fit its cell.

## Operations that disagree with the symbol

A CIF can state its symmetry three ways: listed operations, a Hall symbol and a Hermann–Mauguin symbol. gemmi takes the first it can match to a tabulated setting, so when the operations are in a setting it does not know (a shifted origin, other axes), it silently falls back to the symbol. Earlier builds then expanded the atoms with the symbol's operations. `dev/ops_fallback_scan.py` finds 1,833 COD CIFs whose sources disagree or that gemmi cannot match, and `dev/sym_conflict_check.py` builds each of them and compares the unit-cell content with the CIF's own formula × Z (hydrogen left out, as it is often not located), or without a formula its density (`_exptl_crystal_density_diffrn`):

- **Operations with a shifted origin or other axes (102 CIFs, e.g. [2014900](https://www.crystallography.net/cod/2014900.cif), [9011493](https://www.crystallography.net/cod/9011493.cif), [1545190](https://www.crystallography.net/cod/1545190.cif)).** Expanded with the symbol's operations, atoms on special positions of the listed setting fall on general positions: the published build had up to twice the atoms (density 2.00 × the CIF's). Now the listed operations (else the Hall symbol) define the group, and all 102 agree with the CIF's formula or density. The same holds for a Hall symbol gemmi does not know next to a readable symbol ([2001219](https://www.crystallography.net/cod/2001219.cif), which is also [corrected from its article](../corrections.md)).
- **Operations and symbol naming different groups.** The operations still decide, unless they contradict the CIF's density or formula and the symbol does not: then the symbol is used and the build report says so (9 CIFs, e.g. [2003052](https://www.crystallography.net/cod/2003052.cif) lists only the two P2₁ operations of P2₁/c). When neither agrees, the CIF is skipped (8 CIFs; origin choice 1 does not fit either). Their articles resolve 7: the symmetry they state plus occupancies the CIF leaves out (disordered or split sites, site-multiplicity-scaled occupancies, a mixed O/F site) are [corrections](../corrections.md) (2003764, 2004899, 2005933, 2100940, 2100941, 2101194, 2101195); 2000974 is built in P2₁2₁2₁, which its article confirms, without its unpublished ethanol of crystallization, as in COD24 ([corrections.md](../corrections.md)).

No CIF that agreed with its density or formula before disagrees now.

## Where COD24 differs

- `R -3 m HR` and `R -3 HR` (1010022, 1010098) describe rhombohedral structures in the reverse setting. COD24 stores them as primitive `P -3 m 1` and `P -3`, which drops the rhombohedral centring and adds reflections the crystal does not have.
- `P 32 2 1 S` (1010918, 1011159, 1011176): the listed operations are those of P 3₂21 with the origin shifted by c/3. COD24 stores P 3₁21, the mirror-image space group. The powder pattern is the same, but the space group is not.

## All 81 CIFs

| COD | CIF symbol | hsrdb-tools | COD24 | Same as ours |
|---|---|---|---|---|
| [1010022](https://www.crystallography.net/cod/1010022.cif) | `R -3 m HR` | R -3 m (166) | P -3 m 1 | **different** |
| [1010025](https://www.crystallography.net/cod/1010025.cif) | `B 1 21/d 1` | P 1 21/c 1 (14) | P 1 21/c 1 | same |
| [1010098](https://www.crystallography.net/cod/1010098.cif) | `R -3 HR` | R -3 (148) | P -3 | **different** |
| [1010430](https://www.crystallography.net/cod/1010430.cif) | `P 31 2 1 S` | P 31 2 1 (152) | P 31 2 1 | same |
| [1010918](https://www.crystallography.net/cod/1010918.cif) | `P 32 2 1 S` | P 32 2 1 (154) | P 31 2 1 | **different** |
| [1011012](https://www.crystallography.net/cod/1011012.cif) | `B 1 21/d 1` | P 1 21/c 1 (14) | P 1 21/c 1 | same |
| [1011025](https://www.crystallography.net/cod/1011025.cif) | `P 63 m c S` | P 63 m c (186) | P 63 m c | same |
| [1011159](https://www.crystallography.net/cod/1011159.cif) | `P 32 2 1 S` | P 32 2 1 (154) | P 31 2 1 | **different** |
| [1011176](https://www.crystallography.net/cod/1011176.cif) | `P 32 2 1 S` | P 32 2 1 (154) | P 31 2 1 | **different** |
| [1100106](https://www.crystallography.net/cod/1100106.cif) | `C2:b1` | C 1 2 1 (5) | C 1 2 1 | same |
| [1528279](https://www.crystallography.net/cod/1528279.cif) | `(no symbol, operations only)` | R -3 m (166) | not included |  |
| [1537595](https://www.crystallography.net/cod/1537595.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [1542047](https://www.crystallography.net/cod/1542047.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [1551314](https://www.crystallography.net/cod/1551314.cif) | `P 1 21/m 1 (c,2*a+c,b)` | P 1 21/m 1 (11) | not included |  |
| [1573719](https://www.crystallography.net/cod/1573719.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [2003832](https://www.crystallography.net/cod/2003832.cif) | `F 41/a d c` | I 41/a c d (142) | not included |  |
| [2022089](https://www.crystallography.net/cod/2022089.cif) | `P 1 21/c 1 (c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2100940](https://www.crystallography.net/cod/2100940.cif) | `F d d d {origin @ -1 @ d d d}` | F d d d (70) | F d d d | same |
| [2100941](https://www.crystallography.net/cod/2100941.cif) | `F d d d {origin @ -1 @ d d d}` | F d d d (70) | F d d d | same |
| [2103867](https://www.crystallography.net/cod/2103867.cif) | `R 1 2/c 1` | C 1 2/c 1 (15) | not included |  |
| [2105506](https://www.crystallography.net/cod/2105506.cif) | `X c` | C 1 c 1 (9) | not included |  |
| [2105936](https://www.crystallography.net/cod/2105936.cif) | `Xmc21` | P m c 21 (26) | not included |  |
| [2107405](https://www.crystallography.net/cod/2107405.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [2108010](https://www.crystallography.net/cod/2108010.cif) | `(no symbol, operations only)` | P -1 (2) | not included |  |
| [2200614](https://www.crystallography.net/cod/2200614.cif) | `P b -3` | P a -3 (205) | not included |  |
| [2300039](https://www.crystallography.net/cod/2300039.cif) | `B 21/c` | P 1 21/c 1 (14) | not included |  |
| [2300547](https://www.crystallography.net/cod/2300547.cif) | `P 1 2/c 1 (a,2*b,c)` | P 1 2/c 1 (13) | not included |  |
| [2310735](https://www.crystallography.net/cod/2310735.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310736](https://www.crystallography.net/cod/2310736.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310737](https://www.crystallography.net/cod/2310737.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310738](https://www.crystallography.net/cod/2310738.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310739](https://www.crystallography.net/cod/2310739.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310740](https://www.crystallography.net/cod/2310740.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310741](https://www.crystallography.net/cod/2310741.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310742](https://www.crystallography.net/cod/2310742.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310743](https://www.crystallography.net/cod/2310743.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310744](https://www.crystallography.net/cod/2310744.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310745](https://www.crystallography.net/cod/2310745.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310746](https://www.crystallography.net/cod/2310746.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310747](https://www.crystallography.net/cod/2310747.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310748](https://www.crystallography.net/cod/2310748.cif) | `P 1 21/c 1 (2*c,2*a+c,b)` | P 1 21/c 1 (14) | not included |  |
| [2310762](https://www.crystallography.net/cod/2310762.cif) | `C m m 2 (2*c,a,b)` | C m m 2 (35) | not included |  |
| [2310763](https://www.crystallography.net/cod/2310763.cif) | `P 4 b m (a,b,2*c)` | P 4 b m (100) | not included |  |
| [2310764](https://www.crystallography.net/cod/2310764.cif) | `P 4 b m (a,b,2*c)` | P 4 b m (100) | not included |  |
| [2310775](https://www.crystallography.net/cod/2310775.cif) | `P m m a (2*b+1/4,c,a-1/3)` | P m m a (51) | not included |  |
| [2310776](https://www.crystallography.net/cod/2310776.cif) | `P m m a (2*b+1/4,c,a-1/3)` | P m m a (51) | not included |  |
| [2310777](https://www.crystallography.net/cod/2310777.cif) | `P m m a (2*b+1/4,c,a-1/3)` | P m m a (51) | not included |  |
| [2310778](https://www.crystallography.net/cod/2310778.cif) | `P m m a (2*b+1/4,c,a-1/3)` | P m m a (51) | not included |  |
| [2310779](https://www.crystallography.net/cod/2310779.cif) | `P m m a (2*b,c,a)` | P m m a (51) | not included |  |
| [2311739](https://www.crystallography.net/cod/2311739.cif) | `P 4 b m (a,b,2*c)` | P 4 b m (100) | not included |  |
| [2311740](https://www.crystallography.net/cod/2311740.cif) | `P 4 b m (a,b,2*c)` | P 4 b m (100) | not included |  |
| [3000225](https://www.crystallography.net/cod/3000225.cif) | `B 1 21/m 1` | P 1 21/m 1 (11) | P 1 21/m 1 | same |
| [4002628](https://www.crystallography.net/cod/4002628.cif) | `P m m m (2*a,2*b,c)` | P m m m (47) | not included |  |
| [7059153](https://www.crystallography.net/cod/7059153.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [7059154](https://www.crystallography.net/cod/7059154.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [7059155](https://www.crystallography.net/cod/7059155.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [7059156](https://www.crystallography.net/cod/7059156.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [7222424](https://www.crystallography.net/cod/7222424.cif) | `(no symbol, operations only)` | P 1 21/m 1 (11) | not included |  |
| [7222425](https://www.crystallography.net/cod/7222425.cif) | `(no symbol, operations only)` | P 1 21/m 1 (11) | not included |  |
| [7222426](https://www.crystallography.net/cod/7222426.cif) | `(no symbol, operations only)` | P 1 21/m 1 (11) | not included |  |
| [7222427](https://www.crystallography.net/cod/7222427.cif) | `(no symbol, operations only)` | P 1 21/m 1 (11) | not included |  |
| [7222428](https://www.crystallography.net/cod/7222428.cif) | `(no symbol, operations only)` | P 1 21/m 1 (11) | not included |  |
| [7701112](https://www.crystallography.net/cod/7701112.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [7703174](https://www.crystallography.net/cod/7703174.cif) | `P 1 21/c 1 (2*a+c,b,c)` | P 1 21/c 1 (14) | not included |  |
| [9000018](https://www.crystallography.net/cod/9000018.cif) | `B 1 21/d 1` | P 1 21/c 1 (14) | P 1 21/c 1 | same |
| [9000043](https://www.crystallography.net/cod/9000043.cif) | `P 1 2/c 1 (c,2*a+c,b)` | P 1 2/c 1 (13) | not included |  |
| [9000069](https://www.crystallography.net/cod/9000069.cif) | `B 1 21/m 1` | P 1 21/m 1 (11) | P 1 21/m 1 | same |
| [9000279](https://www.crystallography.net/cod/9000279.cif) | `F 1 2/d 1` | C 1 2/c 1 (15) | C 1 2/c 1 | same |
| [9004170](https://www.crystallography.net/cod/9004170.cif) | `F 1 2/d 1` | C 1 2/c 1 (15) | C 1 2/c 1 | same |
| [9007431](https://www.crystallography.net/cod/9007431.cif) | `F -4 d 2` | I -4 2 d (122) | I -4 2 d | same |
| [9009387](https://www.crystallography.net/cod/9009387.cif) | `F 1 1 2` | C 1 2 1 (5) | C 1 2 1 | same |
| [9009394](https://www.crystallography.net/cod/9009394.cif) | `F 1 1 2` | C 1 2 1 (5) | C 1 2 1 | same |
| [9009395](https://www.crystallography.net/cod/9009395.cif) | `F 1 1 2` | C 1 2 1 (5) | C 1 2 1 | same |
| [9009396](https://www.crystallography.net/cod/9009396.cif) | `F 1 1 2` | C 1 2 1 (5) | C 1 2 1 | same |
| [9009397](https://www.crystallography.net/cod/9009397.cif) | `F 1 1 2` | C 1 2 1 (5) | C 1 2 1 | same |
| [9009398](https://www.crystallography.net/cod/9009398.cif) | `F 1 1 2` | C 1 2 1 (5) | C 1 2 1 | same |
| [9010978](https://www.crystallography.net/cod/9010978.cif) | `F -4 d 2` | I -4 2 d (122) | I -4 2 d | same |
| [9011409](https://www.crystallography.net/cod/9011409.cif) | `B 1 21/d 1` | P 1 21/c 1 (14) | P 1 21/c 1 | same |
| [9011487](https://www.crystallography.net/cod/9011487.cif) | `C 1 1 21/d` | P 1 21/c 1 (14) | P 1 21/c 1 | same |
| [9012076](https://www.crystallography.net/cod/9012076.cif) | `B 1 21/d 1` | P 1 21/c 1 (14) | P 1 21/c 1 | same |
| [9013163](https://www.crystallography.net/cod/9013163.cif) | `F 1 1 2` | C 1 2 1 (5) | C 1 2 1 | same |
