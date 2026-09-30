"""Does COD24 keep physically invalid anisotropic tensors? Compares its lines with ours (tensor replaced by the
isotropic value) and with cctbx's reading of the CIF (tensor as given).

Runs in the cctbx environment, in the repository root with COD24_HS4x.hsrdb and the COD mirror:
  /opt/cctbx/bin/python dev/invalid_adp_check.py COD_ID...
"""

import sqlite3
import struct
import sys
from pathlib import Path

import numpy as np
from iotbx import cif as icif

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_cctbx as v  # noqa: E402
from hsrdb_tools import entry  # noqa: E402

db = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)


def official(cod):
    row = db.execute("select s.Lines, s.LinesI from general_stored s join general_indexed g on g.id = s.pid "
                     "where g.ProductID = ?", (f"96{cod + 1}",)).fetchone()
    if row is None:
        return None
    n = len(row[1]) // 4
    return (np.array([struct.unpack_from("<f", row[0], 6 * k)[0] for k in range(n)]),
            np.array(struct.unpack(f"<{n}f", row[1])))


def median_error(d, inten, ref):
    """Median intensity difference of COD24 lines of at least 20/1000 matched by d."""
    out = []
    for dk, ik in zip(*ref):
        j = np.argmin(abs(d - dk))
        if abs(d[j] - dk) < 1e-3 * dk and ik >= 20:
            out.append(abs(inten[j] - ik))
    return round(float(np.median(out)), 1) if out else None


v.init(None)
for cod in map(int, sys.argv[1:]):
    path = f"cod-mirror/cif/{str(cod)[0]}/{str(cod)[1:3]}/{str(cod)[3:5]}/{cod}.cif"
    e, ref = entry.from_cif(path), official(cod)
    if ref is None or not e.invalid_adps:
        print(cod, "not in COD24 or no invalid tensor")
        continue
    xs = next(iter(icif.reader(file_path=path).build_crystal_structures().values()))
    v.set_scattering(xs)
    v.default_b(xs, path)
    cd, ci, _ = v.cctbx_lines(xs, min(x.d for x in e.lines))
    ours = median_error(np.array([x.d for x in e.lines]), np.array([x.intensity for x in e.lines]), ref)
    print(f"{cod}\tinvalid tensors {e.invalid_adps}\treplaced {ours}\tas given {median_error(cd, ci, ref)}")
