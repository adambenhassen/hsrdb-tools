"""Check modulated structures (hsrdb_tools.superspace) against the structure factors their refinement program
lists in the CIF itself (JANA writes F_calc or F^2_calc for main and satellite reflections).

For each CIF: the entry as built (or why it is skipped), and, where the CIF lists calculated structure factors
with satellite indices, superspace.agreement: R over all, main and satellite reflections. A calculation that
reproduces the refinement gives R of a few 0.1% (rounding of the listed values); the build skips a CIF above
superspace.AGREEMENT_LIMIT.

    python3 dev/verify_modulated.py LIST [-j N] > verify/modulated.tsv
LIST: CIF paths, one per line.
"""

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry, superspace  # noqa: E402


def one(path):
    cod = Path(path).stem
    try:
        e = entry.from_cif(path)
    except entry.EntryError as exc:
        return f"{cod}\tskipped\t{exc}"
    except Exception as exc:  # report and continue
        return f"{cod}\tERROR\t{type(exc).__name__}: {exc}"
    # the unit-cell content against the CIF's formula x Z and its density (None where the CIF gives neither)
    have = {el: float(t[len(el):]) for el, t in zip(e.elements, e.formula.split())}
    measured = entry._number(e.block, "_exptl_crystal_density_diffrn")
    composition = (f"formula_matches={entry._formula_matches(e.block, have)}\tdensity_ratio="
                   f"{round(e.density / measured, 3) if measured else None}")
    if e.modulated is None:  # basic structure: not modulated, only the average given, or modulation not included
        status = ("basic" if "modulation not included" in e.modulation else
                  "average only" if e.modulation else "not modulated")
        return f"{cod}\t{status}\t{e.modulation}\t{composition}"
    sat = sum(1 for x in e.lines if x.hkl == (0, 0, 0))
    out = f"{cod}\tmodulated\t{e.modulation}\tlines={len(e.lines)}\tsatellite_lines={sat}\t{composition}"
    try:
        c = superspace.agreement(e.modulated, e.block)
    except Exception as exc:
        c = {"check_error": f"{type(exc).__name__}: {exc}"}
    if c:
        out += "\t" + "\t".join(f"{k}={v}" for k, v in c.items())
    return out


def main():
    paths = [line.strip() for line in open(sys.argv[1]) if line.strip()]
    jobs = int(sys.argv[sys.argv.index("-j") + 1]) if "-j" in sys.argv else 1
    with ProcessPoolExecutor(jobs) as pool:  # rows as they finish
        for row in (f.result() for f in as_completed([pool.submit(one, p) for p in paths])):
            print(row, flush=True)


if __name__ == "__main__":
    main()
