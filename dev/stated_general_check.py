"""Atoms exactly on a special position that the CIF states as general: which reading agrees with the CIF.

python dev/stated_general_check.py LIST > out.txt    (LIST: one CIF path per line)
Builds each CIF with such sites twice, images merged (one atom per position) and counted separately
(sites.merge_tolerances separate_stated), and counts which reading agrees with the CIF's density or formula x Z
(entry._composition_matches). CIFs where both readings give the same formula are left out.
"""

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import entry  # noqa: E402


def status(e):
    m = entry._composition_matches(e)
    return "match" if m else "unknown" if m is None else "mismatch"


def main():
    counts, examples = collections.Counter(), collections.defaultdict(list)
    for path in Path(sys.argv[1]).read_text().split():
        try:
            merged = entry._from_cif(path)
            if not merged.stated_general:
                continue
            separate = entry._from_cif(path, separate_stated=True)
        except entry.EntryError:
            continue
        if merged.formula == separate.formula:
            continue
        key = (status(merged), status(separate))
        counts[key] += 1
        examples[key].append(Path(path).stem)
    for (m, s), n in counts.most_common():
        print(f"merged {m:8} separate {s:8} {n:5}  e.g. {' '.join(examples[(m, s)][:8])}")


if __name__ == "__main__":
    main()
