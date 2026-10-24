"""Regenerate tests/data/cod/expected.json from the current code and report what changed.

python dev/update_expected.py

Every tests/data/cod/*.cif is summarised with tests/cod_summary.summarize (the function the regression test uses).
The file is rewritten, and each id whose summary differs from the stored one (beyond the test's tolerances) is
printed with the fields that moved.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "tests")]
from cod_summary import differences, summarize  # noqa: E402

DATA = ROOT / "tests" / "data" / "cod"
TARGET = DATA / "expected.json"


def main():
    old = json.loads(TARGET.read_text(encoding="utf-8")) if TARGET.exists() else {}
    new = {p.stem: summarize(p) for p in sorted(DATA.glob("*.cif"))}
    changed = 0
    for cod_id in sorted(set(old) | set(new)):
        if cod_id not in old:
            print(f"{cod_id}: new")
        elif cod_id not in new:
            print(f"{cod_id}: CIF removed")
        else:
            diff = differences(old[cod_id], new[cod_id])
            if diff:
                print(f"{cod_id}: changed")
                for field, text in diff.items():
                    print(f"    {field}: {text}")
            else:
                continue
        changed += 1
    TARGET.write_text(json.dumps(new, indent=1) + "\n", encoding="utf-8")
    print(f"{len(new)} CIFs, {changed} changed; wrote {TARGET.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
