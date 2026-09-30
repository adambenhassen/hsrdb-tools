"""Site-by-site view of one CIF for dev/verify_cctbx.py discrepancies (cctbx environment).

python dev/verify_detail.py CIF
"""
import sys
from pathlib import Path

import numpy as np
from iotbx import cif as icif

import verify_cctbx  # noqa: F401  (applies the site tolerance to the cctbx CIF builder)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from hsrdb_tools import entry  # noqa: E402
import verify_cctbx as v  # noqa: E402

path = sys.argv[1]
e = entry.from_cif(path)
st = e.structure
mult = e.images.multiplicities(len(st.sites))
print("ours:", st.spacegroup.xhm(), st.cell.parameters, "sites", len(st.sites), "images", len(e.images.site))
for i, s in enumerate(st.sites):
    a = s.aniso
    print(f"  {s.label:8} {s.element.name:3} occ={s.occ:.4f} uiso={s.u_iso:.4f} aniso={a.nonzero()} "
          f"U=({a.u11:.4f},{a.u22:.4f},{a.u33:.4f}) mult={mult[i]} xyz=({s.fract.x:.4f},{s.fract.y:.4f},{s.fract.z:.4f})")
xs = next(iter(icif.reader(file_path=path).build_crystal_structures().values()))
v.set_scattering(xs)
print("cctbx:", xs.space_group_info().symbol_and_number(), xs.unit_cell().parameters(), "sites", xs.scatterers().size())
for sc in xs.scatterers():
    u = sc.u_iso if not sc.flags.use_u_aniso() else sc.u_iso_or_equiv(xs.unit_cell())
    print(f"  {sc.label:8} {sc.scattering_type:3} occ={sc.occupancy:.4f} u={u:.4f} aniso={sc.flags.use_u_aniso()} "
          f"mult={sc.multiplicity()} xyz=({sc.site[0]:.4f},{sc.site[1]:.4f},{sc.site[2]:.4f})")
d_last = min(x.d for x in e.lines)
cd, ci, _ = v.cctbx_lines(xs, d_last)
od = np.array([x.d for x in e.lines])
oi = np.array([x.intensity for x in e.lines])
print("strongest lines ours:", [(round(d, 4), round(i, 1)) for d, i in sorted(zip(od, oi), key=lambda t: -t[1])[:6]])
print("strongest lines cctbx:", [(round(d, 4), round(i, 1)) for d, i in sorted(zip(cd, ci), key=lambda t: -t[1])[:6]])
