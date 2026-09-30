"""CIFs whose symmetry sources disagree, or that gemmi cannot match to a table setting.

python dev/ops_fallback_scan.py LIST [-j N] > out.tsv    (LIST: one CIF path per line)

A CIF can state its symmetry three ways: listed operations (S), a Hall symbol (H) and a Hermann-Mauguin symbol (M).
gemmi takes the first it can match to a table setting and ignores the others, e.g. a Hall symbol with a change of
basis that no table setting has, next to a standard H-M symbol. Per CIF where the sources give different operation
sets (compared exactly, translations modulo 1), or gemmi finds no group: the three sources, the number of
operations of each, whether the listed operations form a group, gemmi's choice, and whether the listed operations
are a known group with a shifted origin (setting.identify_shifted).
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import setting  # noqa: E402


def _keys(ops):
    den = gemmi.Op.DEN
    return frozenset((tuple(map(tuple, o.rot)), tuple(t % den for t in o.tran)) for o in ops)


def _hall(block):
    for tag in ("_space_group_name_Hall", "_symmetry_space_group_name_Hall"):
        v = block.find_value(tag)
        if v is not None and gemmi.cif.as_string(v).strip() not in ("", "?", "."):
            return gemmi.cif.as_string(v).strip()
    return ""


def check(path):
    try:
        block = gemmi.cif.read(path).sole_block()
        st = gemmi.make_small_structure_from_block(block)
    except Exception:  # unreadable CIF: not what this scan is about
        return None
    hall, hm = _hall(block), st.spacegroup_hm.strip()
    sets = {}
    if st.symops:
        try:
            ops = [gemmi.Op(o) for o in st.symops]
            sets["S"] = _keys(ops) if setting._is_group(ops) else "not_a_group"
        except RuntimeError:
            sets["S"] = "unparsed"
    if hall:
        try:
            sets["H"] = _keys(gemmi.symops_from_hall(hall))
        except RuntimeError:
            sets["H"] = "unparsed"
    if hm and hm not in ("?", "."):
        probe = gemmi.SmallStructure()
        probe.spacegroup_hm = setting.normalise_symbol(hm)
        probe.determine_and_set_spacegroup("12")
        sets["M"] = _keys(probe.spacegroup.operations()) if probe.spacegroup else "unknown"
    groups = {k: v for k, v in sets.items() if not isinstance(v, str)}
    agree = len(set(groups.values())) <= 1
    if st.spacegroup is not None and agree:
        return None
    chosen = _keys(st.spacegroup.operations()) if st.spacegroup else None
    source = "".join(k for k, v in groups.items() if v == chosen) or "-"
    shifted = None
    if isinstance(sets.get("S"), frozenset) and st.spacegroup is None:
        shifted = setting.identify_shifted(list(st.symops))

    def desc(k):
        v = sets.get(k)
        return "" if v is None else v if isinstance(v, str) else str(len(v))

    return "\t".join(map(str, (
        path, repr(hm), repr(hall), len(st.symops), desc("S"), desc("H"), desc("M"), int(agree),
        st.spacegroup.xhm() if st.spacegroup else "", source,
        f"{shifted[0].xhm()} {tuple(float(x) for x in shifted[1])}" if shifted else "")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("list")
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()
    paths = [s for s in Path(args.list).read_text().split() if s]
    print("path\thm\thall\tn_symops\tS\tH\tM\tagree\tgemmi\tgemmi_source\tshifted", flush=True)
    with ProcessPoolExecutor(args.j) as ex:
        for row in ex.map(check, paths, chunksize=256):
            if row:
                print(row, flush=True)


if __name__ == "__main__":
    main()
