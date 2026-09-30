"""Show the individual reflections behind one line where we and COD24 disagree (COD 9009599, P 63/m)."""
import sqlite3, struct, sys
from pathlib import Path
import gemmi, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry, pattern, sites

cod = int(sys.argv[1]) if len(sys.argv) > 1 else 9009599
f = Path(f"cod-mirror/cif/{str(cod)[0]}/{str(cod)[1:3]}/{str(cod)[3:5]}/{cod}.cif")
e = entry.from_cif(f)
st, im = e.structure, e.images
hkl = np.asarray(gemmi.make_miller_array(st.cell, st.spacegroup, pattern.d_for_two_theta(60)), dtype=np.int64)
d = 1 / np.linalg.norm(hkl @ np.array(st.cell.frac.mat.tolist()), axis=1)
inten = pattern._intensities(st, hkl, d, pattern.CU_KA1, im)
db = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
L, H, I = db.execute("select s.Lines, s.HKL, s.LinesI from general_indexed g join general_stored s on s.pid=g.id "
                     "where g.ProductID=?", (str(int(f"96{cod + 1}")),)).fetchone()
n = len(I) // 4
ref = [(struct.unpack_from("<f", L, 6 * i)[0], struct.unpack(f"<{n}f", I)[i], struct.unpack_from("<hhh", H, 6 * i))
       for i in range(n)]
ours_lines = {round(x.d, 4): x for x in e.lines}
# scale our absolute intensities to COD24's on a line where both have a single reflection
order = np.argsort(-d)
print(f"COD {cod}, {st.spacegroup.hm}, Laue {st.spacegroup.laue_str()}")
print(f"{'d':>8} {'hkl':>12} {'our |F|^2 LP m (single refl)':>30}")
groups = {}
for i in order:
    groups.setdefault(round(d[i], 4), []).append(i)
for dk, idx in list(groups.items())[:14]:
    refl = ", ".join(f"{tuple(int(v) for v in hkl[i])}: {inten[i]:.4g}" for i in idx)
    cod24 = [r for r in ref if abs(r[0] - dk) < 2e-4]
    ours = ours_lines.get(dk)
    print(f"{dk:8.4f}  COD24 {', '.join(f'{r[2]} I={r[1]:.1f}' for r in cod24) or '-':30}  ours I={ours.intensity if ours else 0:7.1f}   reflections: {refl}")
