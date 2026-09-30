"""Space-group information of the CIFs skipped with 'unknown space group': symbols and listed operations.

python dev/unknown_sg.py SKIPPED.tsv > unknown_sg.tsv
"""
import sys
from pathlib import Path

import gemmi

TAGS = ["_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt", "_symmetry_space_group_name_Hall",
        "_space_group_name_Hall", "_symmetry_Int_Tables_number", "_space_group_IT_number"]
print("path\tn_ops\t" + "\t".join(TAGS))
for line in Path(sys.argv[1]).read_text().splitlines()[1:]:
    path, reason = line.split("\t", 1)
    if "unknown space group" not in reason:
        continue
    block = gemmi.cif.read(path).sole_block()
    ops = block.find_values("_space_group_symop_operation_xyz") or block.find_values("_symmetry_equiv_pos_as_xyz")
    vals = [gemmi.cif.as_string(block.find_value(t) or "") for t in TAGS]
    print(f"{path}\t{len(ops)}\t" + "\t".join(v.replace("\t", " ") for v in vals))
