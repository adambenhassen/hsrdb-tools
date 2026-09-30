"""The reference-setting transform must leave the powder pattern and Ueq unchanged."""
import sys, collections
from pathlib import Path
import gemmi, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import pattern, setting, sites
res = collections.Counter(); shown = 0
for p in sorted(Path("scratch/cif_cache").glob("*.cif")):
    st = gemmi.read_small_structure(str(p))
    if st.spacegroup is None or st.spacegroup.basisop == gemmi.Op("x,y,z"): continue
    if st.cell.volume > 6000: continue
    new, changed = setting.to_reference(st)
    a, _ = pattern.calculate(st, 60); b, _ = pattern.calculate(new, 60)
    ka = sorted((x.d, x.intensity) for x in a if x.d > 1.6); kb = sorted((x.d, x.intensity) for x in b if x.d > 1.6)
    ok = len(ka) == len(kb) and all(abs(u[0] - v[0]) < 1e-4 and abs(u[1] - v[1]) < 0.5 for u, v in zip(ka, kb))
    ueq = lambda s, c: np.trace(sites.ucif_to_cart(s.aniso, c)) / 3
    ok_u = all(abs(ueq(s, st.cell) - ueq(t, new.cell)) < 1e-6 for s, t in zip(st.sites, new.sites) if s.aniso.nonzero())
    key = (st.spacegroup.xhm(), "->", new.spacegroup.xhm(), "pattern " + ("same" if ok else "DIFF"), "Ueq " + ("same" if ok_u else "DIFF"))
    res[key] += 1
    if not ok and shown < 2:
        shown += 1; print(p.stem, ka[:5], kb[:5])
    if sum(res.values()) >= 150: break
for k, v in res.most_common(): print(v, k)
