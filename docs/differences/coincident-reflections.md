# Coincident reflections

**COD24:** where two reflections that are not related by symmetry fall at exactly the same angle, the line intensity appears to include only one of them.
**hsrdb-tools:** adds all reflections that fall on a line.

## Why

A powder line is the sum of every reflection diffracting at that angle. In Laue classes −3, 6/m, 4/m and m−3, pairs such as (h k l) and (k h l) are not symmetry-equivalent (their intensities differ) but always have the same d-spacing. Pseudo-symmetric cells produce the same situation by accident.

## Evidence: COD 9009599 (P 6₃/m, Laue class 6/m)

`python dev/coincident_proof.py 9009599` (run in the repository root, with `COD24_HS4x.hsrdb` and the COD mirror present) lists every reflection behind the first lines:

| d (Å) | Reflections | Intensity of each (ours, absolute) | COD24 line | Ours | First reflection alone, on COD24's scale |
|---|---|---|---|---|---|
| 3.2809 | (0 1 2) | 1.55·10⁷ | 246.5 | 224.7 | 246.5 (reference line) |
| 3.3787 | (1 2 0) + (2 1 0) | 1.348·10⁷ + 1.086·10⁷ | 215.7 | 352.9 | 214.4 |
| 3.0472 | (1 2 1) + (2 1 1) | 4.49·10⁷ + 2.408·10⁷ | 711.1 | 1000.0 | 714 |

On lines with a single reflection, COD24 and hsrdb-tools agree after one common scale factor (COD24 is 1.097 times ours, because its strongest line is different). On the two lines with two reflections, COD24's value equals the first reflection alone: the second one is missing. Such a line is then up to 40% too weak, and because the strongest line of the pattern changes, all intensities are rescaled by about 9%.

## Extent

Rate of entries with a strong line more than 150/1000 away from COD24 (`dev/outlier_by_laue.py`): Laue class −3 38%, 6/m 35%, 4/m 33%, m−3 30%, −3m 22%, 2/m with β close to 90° 19%, other classes 3–8%. The classes where non-equivalent reflections coincide are exactly the ones with high rates.
