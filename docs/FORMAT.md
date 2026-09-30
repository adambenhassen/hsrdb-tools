# The HighScore `.hsrdb` reference database format

Decoded from the official COD databases for HighScore (COD Oct 2014 for 3.x, COD24 for 4.x/5.x). Everything is little-endian. HighScore labels these files "HighScore(Plus) V3.X database" in both cases; the differences are listed per version.

## Container

A SQLite 3 database, page size 1024. `hsrdb_tools/schema.sql` (4.x) and `hsrdb_tools/schema_hs3.sql` (3.x) are the exact schemas. Rows are inserted through the `general` and `literature` views, whose `INSTEAD OF` triggers fill the underlying tables (including the FTS3 text indexes) and increment `Properties.NoOfPatterns`.

`Properties`:

| Key | Value |
|---|---|
| `NoOfPatterns` | number of patterns |
| `ReferenceCodePrefixes` | every `96-xxx` prefix present, one per line. HighScore rejects two active databases that share a prefix |
| `ID` | `CIFDB` for databases converted from CIF files |

`Subfiles`: 0 User Inorganic, 1 User Organic, 2 User Mineral, 3 User Metallic (then Amorphous, Ceramic, … unused by COD). `SubFileRef(pid, SubFileID)` assigns them; an entry can be in several.

## Pattern rows

| Field | Content |
|---|---|
| `ProductID` | reference code as text: `96` + (COD id + 1), e.g. COD 1000050 → `961000051`, shown as 96-100-0051 |
| `XTSLSYS` | A triclinic, M monoclinic, O orthorhombic, T tetragonal, H trigonal/hexagonal, C cubic |
| `A`…`GAMMA` | cell, float32 precision |
| `XTLASPECTN` | space group number |
| `SPGR` | Hermann-Mauguin symbol of the reference setting, e.g. `P 1 21/c 1`, `F d -3 m` (origin choice 2), `R -3 m` (hexagonal axes) |
| `DX` | calculated density; `DM` 0 (not measured); `D` = `DM` if > 0 else `DX`, set by the insert trigger |
| `IIC` | I/Ic (reference intensity ratio against corundum, 50:50 by weight) |
| `QUALFINAL` | `=` (calculated) |
| `Comment` | `Crystal color: …`, `Crystal description: …`, `Publication title: …`, `COD database code: N`, separated by CR |
| `chemicalformula` | unit-cell contents, `El1n.nn El2n.nn …` in order of first appearance |
| `e1`, `e2` | element bitmasks: bit Z in `e1` for Z < 64, bit Z−64 in `e2`; signed 64-bit |

Peak list, strongest line first, up to 203 lines:

| Blob | Per line |
|---|---|
| `Lines` | float32 d (Å), uint16 intensity (0–1000) |
| `HKL` | int16 h, k, l |
| `LinesI` (4.x only) | float32 intensity |

`pattern_strongestlines.Lines` is the first 60 bytes (10 lines) of `Lines`.

Lines are calculated for Cu Kα1 (1.540598 Å) up to 90° 2θ (140° if that gives fewer than 10 lines). A reflection within 0.005° 2θ of the previous one joins its line, so a merged line can span more than 0.005°. Intensity is multiplicity × |F|² × powder Lorentz-polarisation (1 + cos²2θ)/(sin²θ cos θ), with anomalous scattering (|F|² averaged over Friedel pairs) and the published displacement parameters (anisotropic tensors rotated with each symmetry image; B = 0.5 Å² where none is given, as in the official databases). The README lists the full calculation model.

## Structure records

`Phases(Name, Phase)` with `Name` = `PDF:96-xxx-xxxx`, linked from `PatternPhase(ReferenceCode, PhaseName)`.

Blob: uint16 `8`, uint32 CRC-32 of the payload, raw-deflated payload.

Payload, with "S" a string (int32 character count, then UTF-16LE in version 6 or Windows-1252 in version 3):

```
int32  record version            6 (HighScore 4.x+) or 3 (HighScore 3.x)
TRV    Scale Factor
int32  3 (v6) / 2 (v3)
TRV×6  Cell a, b, c [Å], alpha, beta, gamma [°]
u8     flag (1)
int32  crystal system: 1 triclinic, 3 monoclinic, 5 orthorhombic, 6 tetragonal, 8 hexagonal, 9 cubic, 10 trigonal
25     zero bytes
f64    -1.0                      (v6 only)
int32  2
int32  number of atoms
per atom:
  int32  4 (v6) / 2 (v3)
  TRV×3  X, Y, Z
  int32  0
  TRV    Biso, sof
  S      label, element
  int32  charge, site multiplicity
  S      Wyckoff symbol, e.g. "4e"
  TRV×6  B11, B22, B33, B12, B13, B23
  26     bytes; byte 20 = 1 for calculated (riding) atoms, 2 for dummy atoms
19 (v6) / 18 (v3) constant bytes
TRV×19 profile parameters (Caglioti U/V/W, peak shape, asymmetry, TOF …; v3 names differ, e.g. "Profile U Right")
24     zero bytes
TRV×7  preferred orientation ×2, overall B, extinction, absorption, porosity, roughness
       v6: spherical-harmonics block, fixed per Laue class, ending with int32 n, then n TRV coefficients
       v3: 25 fixed bytes, no coefficients
16     zero bytes
f64    R factor
S      phase name
int32  3 (v6) / 2 (v3)
int32  COD id
int32  0
S×14   audit date, audit method, -, -, systematic name, mineral name, common name, formula sum,
       structural formula, title, authors, -, journal, CODEN
int32×4 volume, year, first page, last page
S      -
S      comment (crystal colour / description, CR LF separated)
       106 bytes: journal issue (int32) at offset 24, display colour at offset 61, otherwise constant
       (about 3% of official records carry further optional fields here; hsrdb-tools writes none)
```

TRV (refinable value):

```
int32  2 (v6) / 1 (v3)
u8     refine flag
f64×7  value, esd, last shift, minimum, maximum, maximum − minimum, 0
7      flag bytes
S      "<phase name> <parameter>"
```

Limits, flags and the default profile values are the same in every official record of a version; `hsrdb_tools/templates.json` (v6) and `templates_hs3.json` (v3) hold them, extracted with `dev/extract_templates.py`. B values are 8π²U. Values are kept within the version's limits, as in the official records (e.g. occupancy at most 1).

The space group is not stored in the structure record; HighScore takes it from the pattern's `SPGR`.
