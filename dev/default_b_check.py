"""Which displacement parameter did the official COD24 patterns assume for atoms without one?

python dev/default_b_check.py LIST [N]
For up to N CIFs from LIST whose atoms carry no displacement parameters at all, computes the pattern with U = 0
and with B = 0.5 A^2 and reports which one matches the COD24 intensities (strong lines, matched by d) better.
The two differ most at high angle, where the Debye-Waller factor exp(-B s^2/4) matters.
"""

import sqlite3
import struct
import sys
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry, pattern, sites  # noqa: E402

db = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
limit = int(sys.argv[2]) if len(sys.argv) > 2 else 200


def official(cod_id):
    row = db.execute("select s.Lines, s.LinesI from general_stored s join general_indexed g on g.id = s.pid "
                     "where g.ProductID = ?", (f"96{cod_id + 1}",)).fetchone()
    if row is None:
        return None
    n = len(row[1]) // 4
    return np.array([struct.unpack_from("<f", row[0], 6 * k)[0] for k in range(n)]), np.array(struct.unpack(f"<{n}f", row[1]))


def error(lines, ref):
    d, i = np.array([x.d for x in lines]), np.array([x.intensity for x in lines])
    rd, ri = ref
    errs = []
    for dk, ik in zip(rd, ri):
        j = np.argmin(np.abs(d - dk))
        if abs(d[j] - dk) < 1e-3 * dk and ik >= 50:
            errs.append(abs(i[j] - ik))
    return np.median(errs) if errs else None


done = votes0 = votes05 = 0
for path in Path(sys.argv[1]).read_text().split():
    block = gemmi.cif.read(path).sole_block()
    if any(block.find_values(t) for t in ("_atom_site_U_iso_or_equiv", "_atom_site_B_iso_or_equiv",
                                          "_atom_site_aniso_U_11", "_atom_site_aniso_B_11", "_atom_site_aniso_beta_11")):
        continue
    try:
        e = entry.from_cif(path)
    except Exception:
        continue
    ref = official(e.cod_id)
    if ref is None:
        continue
    st = e.structure
    zero = error(e.lines, ref)
    for s in st.sites:
        s.u_iso = 0.5 / (8 * np.pi**2)
    b05 = error(pattern.stick_pattern(st, images=sites.unit_cell_images(st, site_orders=e.site_orders))[0], ref)
    if zero is None or b05 is None:
        continue
    done += 1
    votes0 += zero < b05
    votes05 += b05 < zero
    print(f"{e.cod_id}\tU=0 {zero:.1f}\tB=0.5 {b05:.1f}", flush=True)
    if done >= limit:
        break
print(f"CIFs {done}: U=0 closer in {votes0}, B=0.5 closer in {votes05}")
