"""Markdown table of skipped CIFs, from the skipped report of a build.

python dev/skipped_md.py OUT.skipped.tsv [--recheck] > docs/SKIPPED.md
With --recheck, every listed CIF is run through the current code first and only those still skipped are listed
(useful before a rebuild). Notes about CIFs that were built are ignored.
"""

import argparse
import collections
import re
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hsrdb_tools import build, entry  # noqa: E402

REASONS = [  # (substring of the error, reason shown), first match wins
    ("no atom coordinates", "No atom coordinates"),
    ("do not match the structure", "Symmetry operations do not fit the cell or structure"),
    ("unknown space group", "No usable symmetry"),
    ("no valid unit cell", "No valid unit cell"),
    ("no X-ray scattering factors for element Cn", "Cyanide given as one pseudo-atom CN"),
    ("atom type 'CN' is not an element", "Cyanide given as one pseudo-atom CN"),
    ("is not an element", "Atom type that is not an element"),
    ("has no element", "Atom type that is not an element"),
]
SYMBOL_TAGS = ("_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt", "_space_group_name_Hall",
               "_symmetry_space_group_name_Hall")


def reason_of(error):
    return next((shown for key, shown in REASONS if key in error), error)


def details(path, reason, error):
    """The CIF's symmetry for symmetry problems, the type symbol for atom types, else nothing."""
    if reason.startswith("Atom type"):
        m = re.search(r"'([^']*)'", error)
        if not m:
            return ""
        return f"`{m.group(1)}`" if "is not an element" in error else f"atom `{m.group(1)}` without type"
    if "symmetry" not in reason.lower():
        return ""
    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except Exception:
        return ""
    ops = block.find_values("_space_group_symop_operation_xyz") or block.find_values("_symmetry_equiv_pos_as_xyz")
    symbols = [gemmi.cif.as_string(v).strip() for v in map(block.find_value, SYMBOL_TAGS) if v is not None]
    symbols = [s for s in symbols if s not in ("", "?", ".")]
    if not symbols and not ops:
        return "no symbol, no operations"
    text = " / ".join(f"`{s}`" for s in dict.fromkeys(symbols)) or "no symbol"
    return f"{text}, {len(ops)} operations" if ops else f"{text}, no operations"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("report")
    ap.add_argument("--recheck", action="store_true")
    args = ap.parse_args()
    rows = []
    for line in Path(args.report).read_text().splitlines()[1:]:
        path, error = line.split("\t", 1)
        if error.startswith("note:"):
            continue
        if args.recheck:
            try:
                entry.from_cif(path)
                continue  # built now
            except build.DATA_ERRORS as exc:
                error = f"{type(exc).__name__}: {exc}"
        reason = reason_of(error)
        rows.append((reason, int(Path(path).stem), details(path, reason, error)))
    counts = collections.Counter(r for r, _, _ in rows)
    print("# Skipped CIFs\n")
    print(f"{len(rows):,} COD CIFs from which no correct pattern can be calculated. The README explains each reason.\n")
    print("| Reason | CIFs |\n|---|---|")
    for reason, n in counts.most_common():
        print(f"| {reason} | {n:,} |")
    print("\n| COD id | Reason | Details | CIF |\n|---|---|---|---|")
    for reason, cod, info in sorted(rows, key=lambda r: (-counts[r[0]], r[0], r[1])):
        print(f"| {cod} | {reason} | {info.replace('|', '\\|')} | [{cod}.cif](https://www.crystallography.net/cod/{cod}.cif) |")


if __name__ == "__main__":
    main()
