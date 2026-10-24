"""Which of our skipped CIFs does the official database contain, and with what?"""
import collections, sqlite3, sys
from pathlib import Path

report = Path(sys.argv[1])
ref = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
by_reason = collections.defaultdict(list)
for line in report.read_text().splitlines()[1:]:
    path, reason = line.split("\t", 1)
    if reason.startswith("note:"):
        continue
    by_reason[reason].append(int(Path(path).stem))
for reason, ids in by_reason.items():
    present = [i for i in ids if ref.execute("select 1 from general_indexed where ProductID=?", (str(int(f"96{i + 1}")),)).fetchone()]
    print(f"{reason}: {len(ids)} skipped, {len(present)} in COD24, e.g. {present[:5]}")
