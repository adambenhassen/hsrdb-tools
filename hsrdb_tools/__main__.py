"""hsrdb-tools: build HighScore reference databases (.hsrdb) from COD CIF files."""

import argparse
import collections
import os
import shutil
import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from . import build, entry, hsrdb

def _open(path, target, base, append):
    """Connection to the database to fill, and the file it is built in until the build succeeds."""
    if append:
        if not path.exists():
            sys.exit(f"{path} does not exist")
        db, work = sqlite3.connect(path), path
    else:
        if path.exists():
            sys.exit(f"{path} already exists")
        work = path.with_name(path.name + ".partial")
        work.unlink(missing_ok=True)
        if base:
            shutil.copyfile(base, work)
            db = sqlite3.connect(work)
        else:
            db = hsrdb.create(work, target)
    found = hsrdb.target_of(db)
    if found != target:
        db.close()
        if work != path:
            work.unlink()
        sys.exit(f"{base or path}: schema is for {found}, not {target}")
    return db, work


def cmd_add(args):
    outputs = {t: Path(p) for t, p in (("hs4", args.hs4), ("hs3", args.hs3)) if p}
    if not outputs:
        sys.exit("give --hs4 and/or --hs3 output paths")
    if args.base and args.append:
        sys.exit("use either --base or --append")
    if args.base and len(outputs) > 1:
        sys.exit("--base extends one database; give a single output with it")
    opened = {t: _open(p, t, args.base, args.append) for t, p in outputs.items()}
    dbs = {t: db for t, (db, _) in opened.items()}
    existing = {t: {int(c) for (c,) in db.execute("select ProductID from general_indexed")} for t, db in dbs.items()}
    report = next(iter(outputs.values())).with_suffix(".skipped.tsv")
    cifs = sorted(Path(p) for src in args.cifs
                  for p in (Path(src).rglob("*.cif") if Path(src).is_dir() else [Path(src)]))
    added = collections.Counter()
    present = internal = 0
    built = {}  # reference code -> CIF added in this run
    problems = collections.Counter()
    ok = False
    try:
        with ProcessPoolExecutor(args.jobs) as pool, open(report, "w") as rep:
            rep.write("file\tproblem\n")
            jobs = [(p, tuple(outputs)) for p in cifs]
            for i, (path, e, lines, blobs, err) in enumerate(pool.map(build.prepare, jobs, chunksize=16), 1):
                if err:
                    kind, detail = err
                    reason = detail if kind == "skipped" else "INTERNAL " + detail.strip().splitlines()[-1]
                    problems[reason] += 1
                    rep.write(f"{path}\t{reason}\n")
                    if kind == "internal error":
                        internal += 1
                        print(f"[{i}/{len(cifs)}] {path.name}: INTERNAL ERROR\n{detail}", flush=True)
                    continue
                if e.dropped_sites:
                    rep.write(f"{path}\tnote: ignored {e.dropped_sites} atom_site rows without element or coordinates\n")
                if e.invalid_adps:
                    rep.write(f"{path}\tnote: replaced {e.invalid_adps} physically invalid displacement parameters\n")
                if e.correction:
                    rep.write(f"{path}\tnote: corrected from the publication: {e.correction}\n")
                if e.symmetry_note:
                    rep.write(f"{path}\tnote: {e.symmetry_note}\n")
                if e.merge_note:
                    rep.write(f"{path}\tnote: {e.merge_note}\n")
                if e.modulation:
                    rep.write(f"{path}\tnote: {e.modulation}\n")
                if e.dummy_sites:
                    rep.write(f"{path}\tnote: {e.dummy_sites} dummy atoms do not scatter\n")
                masks = entry.element_masks(e.elements)
                new = [t for t in dbs if e.reference_code not in existing[t]]
                present += not new
                if e.reference_code in built:
                    rep.write(f"{path}\tnote: not added: same COD id as {built[e.reference_code]}\n")
                elif new:
                    built[e.reference_code] = path
                for t in new:
                    hsrdb.add(dbs[t], e, masks, lines, blobs[t], t)
                    existing[t].add(e.reference_code)
                    added[t] += 1
                if i % 1000 == 0:
                    print(f"[{i}/{len(cifs)}] {dict(added)} added, {sum(problems.values())} skipped", flush=True)
        for db in dbs.values():
            hsrdb.update_prefixes(db)
            db.commit()  # each database is filled in one transaction
        ok = True
    finally:
        for db, work in opened.values():
            if not ok:
                db.rollback()
            db.close()
            if not ok and not args.append:
                work.unlink(missing_ok=True)
    for t, (_, work) in opened.items():
        db = sqlite3.connect(work)
        db.execute("vacuum")
        db.close()
        if work != outputs[t]:
            os.replace(work, outputs[t])
    print(f"added {dict(added)}, already present {present}, skipped {sum(problems.values())} (listed in {report})")
    for reason, n in problems.most_common(15):
        print(f"  {n:7d}  {reason}")
    if internal:
        sys.exit(f"{internal} internal errors: bugs in hsrdb-tools, tracebacks above")


def main():
    ap = argparse.ArgumentParser(prog="hsrdb-tools", description=__doc__)
    sub = ap.add_subparsers(required=True)
    a = sub.add_parser("add", help="add CIFs to .hsrdb databases (new, a copy of --base, or --append in place)")
    a.add_argument("cifs", nargs="+", help="CIF files or directories (searched recursively)")
    a.add_argument("--hs4", metavar="PATH", help="database for HighScore 4.x and later")
    a.add_argument("--hs3", metavar="PATH", help="database for HighScore 3.x")
    a.add_argument("--base", help="existing .hsrdb to copy and extend (single output only)")
    a.add_argument("--append", action="store_true", help="add to the given databases in place")
    a.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="worker processes")
    a.set_defaults(func=cmd_add)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
