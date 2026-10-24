"""hsrdb-tools: build HighScore reference databases (.hsrdb) from COD CIF files."""

import argparse
import collections
import datetime
import hashlib
import json
import os
import platform
import shutil
import sqlite3
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from importlib import metadata
from pathlib import Path

from . import __version__, build, entry, hsrdb

CHECKPOINT = 2000  # CIFs between commits of a new build: the points an interrupted build resumes from
# private tables of an unfinished build, dropped before it is renamed: its run state, and the reference code and CIF
# of each entry it added (for the note on duplicate COD ids; in the first output only)
STATE, BUILT = "hsrdb_tools_resume", "hsrdb_tools_built"


def _state(db):
    """Run state saved by the last checkpoint, or None."""
    if not db.execute("select 1 from sqlite_master where type = 'table' and name = ?", (STATE,)).fetchone():
        return None
    row = db.execute(f"select state from {STATE}").fetchone()
    return json.loads(row[0]) if row else None


def _save_state(db, state):
    db.execute(f"create table if not exists {STATE}(state text not null)")
    db.execute(f"delete from {STATE}")
    db.execute(f"insert into {STATE}(state) values (?)", (json.dumps(state),))


def _open(path, target, base, append, resume):
    """Connection to the database to fill, and the file it is built in until the build succeeds."""
    if append:
        if not path.exists():
            sys.exit(f"{path} does not exist")
        db, work = sqlite3.connect(path), path
    else:
        if path.exists():
            sys.exit(f"{path} already exists")
        work = path.with_name(path.name + ".partial")
        if resume:
            if not work.exists():
                sys.exit(f"{work} does not exist: nothing to resume")
            db = sqlite3.connect(work)
            if _state(db) is None:
                db.close()
                sys.exit(f"{work} holds no interrupted build to resume")
        else:
            if work.exists():
                with open(work, "rb") as f:
                    sqlite = f.read(16) == b"SQLite format 3\0"
                if sqlite:  # anything else is no build of ours
                    old = sqlite3.connect(work)
                    try:
                        interrupted = _state(old) is not None
                    except sqlite3.DatabaseError as e:  # locked or damaged: it may hold a checkpoint
                        sys.exit(f"{work}: cannot be read ({e}); delete it to start a new build")
                    finally:
                        old.close()
                    if interrupted:
                        sys.exit(f"{work} holds an interrupted build: continue it with --resume, or delete it")
                work.unlink()
            if base:
                shutil.copyfile(base, work)
                db = sqlite3.connect(work)
            else:
                db = hsrdb.create(work, target)
    found = hsrdb.target_of(db)
    if found != target:
        db.close()
        if work != path and not resume:
            work.unlink()
        sys.exit(f"{base or path}: schema is for {found}, not {target}")
    return db, work


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def _source_sha256():
    """Hash of the package's code and data: identifies the code that built a database, with or without git."""
    digest = hashlib.sha256()
    for p in sorted(Path(__file__).parent.iterdir()):
        if p.suffix in (".py", ".json", ".sql"):
            digest.update(p.name.encode() + b"\0" + p.read_bytes())
    return digest.hexdigest()


def _git(package, *args):
    try:
        r = subprocess.run(["git", "-C", str(package), *args], capture_output=True, text=True, encoding="utf-8",
                           timeout=10)
    except (OSError, subprocess.SubprocessError):  # no git
        return None
    return r.stdout.strip() if r.returncode == 0 else None


def _provenance():
    """What built a database: tool version, the hash of its code and data, the git commit if the package runs from
    its own checkout, the platform and the versions of Python and of the libraries the patterns depend on."""
    package = Path(__file__).resolve().parent
    top = _git(package, "rev-parse", "--show-toplevel")
    # not the commit of another repository that happens to contain an installed copy (a venv inside a project)
    own = top is not None and (Path(top) / "hsrdb_tools").resolve() == package
    commit = _git(package, "describe", "--always", "--dirty", "--abbrev=12") if own else None
    versions = {}
    for name in ("gemmi", "numpy", "spglib"):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return {"version": __version__, "source_sha256": _source_sha256(), "git": commit,
            "python": platform.python_version(), "packages": versions, "platform": platform.platform()}


def cmd_add(args):
    outputs = {t: Path(p) for t, p in (("hs4", args.hs4), ("hs3", args.hs3)) if p}
    if not outputs:
        sys.exit("give --hs4 and/or --hs3 output paths")
    if args.base and args.append:
        sys.exit("use either --base or --append")
    if args.resume and args.append:
        sys.exit("--append builds in place and has nothing to resume")
    if args.base and len(outputs) > 1:
        sys.exit("--base extends one database; give a single output with it")
    missing = [s for s in args.cifs if not Path(s).exists()]
    if missing:
        sys.exit(f"not found: {', '.join(missing)}")
    opened = {}
    try:
        for t, p in outputs.items():
            opened[t] = _open(p, t, args.base, args.append, args.resume)
    except BaseException:  # a later output cannot be opened: leave nothing behind for the earlier ones
        for t, (db, work) in opened.items():
            db.close()
            if work != outputs[t] and not args.resume:
                work.unlink(missing_ok=True)
        raise
    dbs = {t: db for t, (db, _) in opened.items()}
    first = next(iter(outputs))
    report = outputs[first].with_suffix(".skipped.tsv")
    cifs = sorted(Path(p) for src in args.cifs
                  for p in (Path(src).rglob("*.cif") if Path(src).is_dir() else [Path(src)]))
    run = {"inputs": list(map(str, args.cifs)), "outputs": {t: str(p) for t, p in outputs.items()},
           "base": args.base, "cifs": len(cifs),
           "cifs_sha256": hashlib.sha256(b"\n".join(map(os.fsencode, cifs))).hexdigest(),
           "source_sha256": _source_sha256()}  # one build, one version of the code
    if args.resume:
        state = {t: _state(db) for t, db in dbs.items()}
        offset = state[first].get("report_offset")
        problem = ("--resume needs the same CIFs, outputs, --base and hsrdb-tools code as the interrupted build"
                   if state[first].get("run") != run else
                   f"{report}: the report of the interrupted build is missing or shorter than at its checkpoint"
                   if not (report.exists() and report.stat().st_size >= offset) else None)
        if problem:
            for db in dbs.values():
                db.close()
            sys.exit(problem)
        done, started, resumes = state[first]["done"], state[first]["started"], state[first]["resumes"] + 1
        added = collections.Counter({t: s["added"] for t, s in state.items()})
        present, internal = state[first]["present"], state[first]["internal"]
        problems = collections.Counter(state[first]["problems"])
        built = {code: path for code, path in dbs[first].execute(f"select code, path from {BUILT}")}
        rep = open(report, "r+", encoding="utf-8", errors="backslashreplace")  # file names that are not UTF-8
        rep.truncate(offset)  # lines written after the checkpoint come again
        rep.seek(0, os.SEEK_END)
    else:
        done, started, resumes = 0, _now(), 0
        added, problems = collections.Counter(), collections.Counter()
        present = internal = 0
        built = {}  # reference code -> first CIF added with it in this build
        rep = open(report, "w", encoding="utf-8", errors="backslashreplace")
        rep.write("file\tproblem\n")
    existing = {t: {int(c) for (c,) in db.execute("select ProductID from general_indexed")} for t, db in dbs.items()}
    new_built = []  # (reference code, CIF) added since the last checkpoint
    checkpoints = not args.append  # --append writes in place: nothing is committed until the end
    saved = args.resume  # whether the .partial files hold a checkpoint

    def record(i, path, e, lines, blobs, err):
        nonlocal present, internal
        if err:
            kind, detail = err
            reason = detail if kind == "skipped" else "INTERNAL " + detail.strip().splitlines()[-1]
            problems[reason] += 1
            rep.write(f"{path}\t{reason}\n")
            if kind == "internal error":
                internal += 1
                print(f"[{i}/{len(cifs)}] {path.name}: INTERNAL ERROR\n{detail}", flush=True)
            return
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
            built[e.reference_code] = str(path)
            new_built.append((e.reference_code, str(path)))
        for t in new:
            hsrdb.add(dbs[t], e, masks, lines, blobs[t], t)
            existing[t].add(e.reference_code)
            added[t] += 1

    def checkpoint(i):
        """Commit the first i CIFs; the first output last, so that its state is never ahead of the others'."""
        rep.flush()
        os.fsync(rep.fileno())
        shared = {"run": run, "started": started, "resumes": resumes, "present": present, "internal": internal,
                  "problems": dict(problems), "report_offset": os.fstat(rep.fileno()).st_size}
        dbs[first].execute(f"create table if not exists {BUILT}(code integer primary key, path text not null)")
        dbs[first].executemany(f"insert into {BUILT}(code, path) values (?, ?)", new_built)
        new_built.clear()
        for t in reversed(list(dbs)):
            _save_state(dbs[t], {"done": i, "added": added[t]} | (shared if t == first else {}))
            dbs[t].commit()

    ok = False
    finished = []  # outputs whose private tables are dropped
    try:
        with ProcessPoolExecutor(args.jobs) as pool, rep:
            jobs = [(p, tuple(outputs)) for p in cifs[done:]]
            for i, result in enumerate(pool.map(build.prepare, jobs, chunksize=16), done + 1):
                record(i, *result)
                if checkpoints and i % CHECKPOINT == 0:
                    checkpoint(i)
                    saved = True
                if i % 1000 == 0:
                    print(f"[{i}/{len(cifs)}] {dict(added)} added, {sum(problems.values())} skipped", flush=True)
            if checkpoints:  # the last CIFs too: a failure while committing them leaves a build that resumes
                checkpoint(len(cifs))
                saved = True
        for t in reversed(list(dbs)):  # each database is filled in one transaction, or since its last checkpoint
            db = dbs[t]
            if not db.in_transaction:
                db.execute("begin")
            db.execute(f"drop table if exists {STATE}")
            db.execute(f"drop table if exists {BUILT}")
            hsrdb.update_prefixes(db)
            db.commit()
            finished.append(t)
        patterns = {t: db.execute("select count(*) from general_indexed").fetchone()[0] for t, db in dbs.items()}
        ok = True
    finally:
        for t, (db, work) in opened.items():
            if not ok:
                db.rollback()
            db.close()
            if not ok and work != outputs[t] and not saved:
                work.unlink(missing_ok=True)
        if not ok and saved:
            kept = ", ".join(f"{w} ({w.stat().st_size / 1e6:,.0f} MB)" for _, w in opened.values() if w.exists())
            print(f"interrupted while finishing: kept {kept}, which cannot be resumed; delete them and build again"
                  if finished else
                  f"interrupted: kept {kept}; continue from the last checkpoint with the same command and --resume, "
                  "or delete them", file=sys.stderr, flush=True)
    not_compacted = {}
    for t, (_, work) in opened.items():
        if work != outputs[t]:  # complete: under its name before compacting, which an interruption must not undo
            os.replace(work, outputs[t])
        db = sqlite3.connect(outputs[t])
        try:
            db.execute("vacuum")  # atomic: interrupted or failed, the database stays as it was
        except sqlite3.OperationalError as e:  # e.g. disk full: VACUUM needs up to twice the database's size
            not_compacted[t] = str(e)
        finally:
            db.close()
    manifest = {**_provenance(), "command": args.argv, "started": started, "finished": _now(), "resumes": resumes,
                "inputs": run["inputs"], "cifs": len(cifs), "base": args.base, "append": args.append,
                "outputs": {t: {"path": str(p), "patterns": patterns[t], "added": added[t]} for t, p in outputs.items()},
                "already_present": present, "skipped": sum(problems.values()), "internal_errors": internal,
                "report": str(report), "not_compacted": not_compacted}
    outputs[first].with_suffix(".build.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    print(f"added {dict(added)}, already present {present}, skipped {sum(problems.values())} (listed in {report})")
    for reason, n in problems.most_common(15):
        print(f"  {n:7d}  {reason}")
    for t, error in not_compacted.items():
        print(f"{outputs[t]}: complete but not compacted ({error}); with disk space for up to twice its size, run: "
              f"sqlite3 {outputs[t]} VACUUM", file=sys.stderr)
    if internal:
        sys.exit(f"{internal} internal errors: bugs in hsrdb-tools, tracebacks above")
    if not_compacted:
        sys.exit(1)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="hsrdb-tools", description=__doc__)
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(required=True)
    a = sub.add_parser("add", help="add CIFs to .hsrdb databases (new, a copy of --base, or --append in place)")
    a.add_argument("cifs", nargs="+", help="CIF files or directories (searched recursively)")
    a.add_argument("--hs4", metavar="PATH", help="database for HighScore 4.x and later")
    a.add_argument("--hs3", metavar="PATH", help="database for HighScore 3.x")
    a.add_argument("--base", help="existing .hsrdb to copy and extend (single output only)")
    a.add_argument("--append", action="store_true", help="add to the given databases in place")
    a.add_argument("--resume", action="store_true",
                   help="continue an interrupted build (same arguments) from its last checkpoint")
    a.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="worker processes")
    a.set_defaults(func=cmd_add)
    args = ap.parse_args(argv)
    args.argv = sys.argv[1:] if argv is None else list(argv)
    args.func(args)


if __name__ == "__main__":
    main()
