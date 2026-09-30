"""Count CIFs with physically invalid anisotropic displacement tensors (a negative eigenvalue).

python dev/scan_adp.py CIF_DIR [workers]
"""
import multiprocessing
import sys
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import sites  # noqa: E402


def check(path):
    try:
        st = gemmi.make_small_structure_from_block(gemmi.cif.read(str(path)).sole_block())
    except Exception:
        return None
    worst = 0.0
    for s in st.sites:
        u = s.aniso
        vals = (u.u11, u.u22, u.u33, u.u12, u.u13, u.u23)
        if not s.aniso.nonzero() or not all(np.isfinite(vals)):
            continue
        worst = min(worst, float(np.linalg.eigvalsh(sites.ucif_to_cart(u, st.cell)).min()))
    return (str(path), worst) if worst < 0 else None


if __name__ == "__main__":
    paths = sorted(Path(sys.argv[1]).rglob("*.cif"))
    bad = []
    with multiprocessing.Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as pool:
        for r in pool.imap_unordered(check, paths, chunksize=64):
            if r:
                bad.append(r)
    bad.sort(key=lambda t: t[1])
    for lim in (-1e-3, -1e-2, -1e-1, -1.0):
        print(f"min eigenvalue < {lim:g} A^2: {sum(1 for _, w in bad if w < lim)} CIFs")
    print(f"any negative: {len(bad)} of {len(paths)}")
    for p, w in bad[:10]:
        print(f"  {w:10.4f}  {p}")
    Path("adp_scan.tsv").write_text("".join(f"{p}\t{w}\n" for p, w in bad))
