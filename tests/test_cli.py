"""The command line: outputs, report, guards, internal errors, interrupted and resumed builds, build manifest."""

import functools
import json
import os
import sqlite3
import subprocess
import sys

import pytest

from hsrdb_tools import __main__ as cli
from hsrdb_tools import __version__, hsrdb
from test_edge_cases import NACL_ATOMS, cif


def run_cli(*args):
    # as in CI: file I/O that depends on the platform's default encoding (cp1252 on Windows) is an error
    env = os.environ | {"PYTHONWARNDEFAULTENCODING": "1", "PYTHONWARNINGS": "error::EncodingWarning"}
    return subprocess.run([sys.executable, "-m", "hsrdb_tools", *map(str, args)],
                          capture_output=True, text=True, encoding="utf-8", env=env)


def write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def private_tables(db):
    """The build's own tables (checkpoint state), which a finished database must not hold."""
    return db.execute("select name from sqlite_master where name like 'hsrdb_tools%'").fetchall()


def test_cli_skips_bad_cifs_and_reports(tmp_path):
    src = tmp_path / "cifs"
    src.mkdir()
    write(src / "1234567.cif", cif(NACL_ATOMS))
    write(src / "1234568.cif", cif(NACL_ATOMS, cod=1234568).replace("_cell_length_a 5.6402\n", ""))
    out = tmp_path / "db.hsrdb"
    r = run_cli("add", src, "--hs4", out, "--hs3", tmp_path / "db3.hsrdb", "-j", "1")
    assert r.returncode == 0, r.stderr
    assert "skipped 1" in r.stdout
    report = (tmp_path / "db.skipped.tsv").read_text(encoding="utf-8")
    assert "1234568.cif" in report and "no valid unit cell" in report
    for path in (out, tmp_path / "db3.hsrdb"):
        db = sqlite3.connect(path)
        assert db.execute("select count(*) from general_indexed").fetchone()[0] == 1
        assert db.execute("select Value from Properties where Key='NoOfPatterns'").fetchone()[0] == "1"
        assert not private_tables(db)
        db.close()
    assert not list(tmp_path.glob("*.partial"))
    manifest = json.loads((tmp_path / "db.build.json").read_text(encoding="utf-8"))
    assert manifest["version"] == __version__ and len(manifest["source_sha256"]) == 64
    assert set(manifest["packages"]) == {"gemmi", "numpy", "spglib"} and all(manifest["packages"].values())
    assert manifest["cifs"] == 2 and manifest["skipped"] == 1 and manifest["internal_errors"] == 0
    assert manifest["outputs"]["hs4"]["patterns"] == manifest["outputs"]["hs3"]["patterns"] == 1


def test_cli_append_duplicates_and_guards(tmp_path):
    src = tmp_path / "cifs"
    src.mkdir()
    write(src / "a.cif", cif(NACL_ATOMS))
    write(src / "b.cif", cif(NACL_ATOMS))  # same COD id
    out = tmp_path / "db.hsrdb"
    r = run_cli("add", src, "--hs4", out, "-j", "1")
    assert r.returncode == 0 and "already present 1" in r.stdout
    assert "b.cif\tnote: not added: same COD id as" in (tmp_path / "db.skipped.tsv").read_text(encoding="utf-8")
    r = run_cli("add", src, "--hs4", out, "--append", "-j", "1")
    assert r.returncode == 0 and "added {}" in r.stdout
    assert run_cli("add", src, "--hs4", out, "-j", "1").returncode != 0  # exists
    r = run_cli("add", src, "--hs3", out, "--append", "-j", "1")
    assert r.returncode == 1 and "schema is for hs4, not hs3" in r.stderr  # wrong schema, clean exit
    for args, message in (
            (("--hs4", tmp_path / "c.hsrdb", "--hs3", tmp_path / "d.hsrdb", "--base", out), "--base extends one"),
            (("--hs4", out, "--base", out, "--append"), "either --base or --append"),
            (("--hs4", out, "--append", "--resume"), "--append builds in place and has nothing to resume"),
            ((), "give --hs4 and/or --hs3")):
        r = run_cli("add", src, *args)
        assert r.returncode == 1 and message in r.stderr, (args, r.stderr)
    r = run_cli("add", tmp_path / "missing.cif", "--hs4", out, "--append")
    assert r.returncode == 1 and "not found: " in r.stderr
    r = run_cli("add", src, tmp_path / "missing", "--hs4", tmp_path / "h.hsrdb")
    assert r.returncode != 0 and "not found:" in r.stderr and not list(tmp_path.glob("h.*"))
    (tmp_path / "empty").mkdir()  # a wrong directory: an error, not an empty database
    r = run_cli("add", src, tmp_path / "empty", "--hs4", tmp_path / "i.hsrdb")
    assert r.returncode == 1 and "empty: no .cif files" in r.stderr and not list(tmp_path.glob("i.*"))
    assert not list(tmp_path.glob("*.partial"))


def test_cli_version():
    r = run_cli("--version")
    assert r.returncode == 0 and r.stdout.strip() == f"hsrdb-tools {__version__}"


def test_cli_leaves_no_partial_file_when_a_later_output_exists(tmp_path):
    src = write(tmp_path / "1234567.cif", cif(NACL_ATOMS))
    write(tmp_path / "db3.hsrdb", "")
    r = run_cli("add", src, "--hs4", tmp_path / "db4.hsrdb", "--hs3", tmp_path / "db3.hsrdb", "-j", "1")
    assert r.returncode != 0 and "already exists" in r.stderr
    assert not list(tmp_path.glob("*.partial")) and not (tmp_path / "db4.hsrdb").exists()


def test_cli_base_copies_and_extends(tmp_path):
    src = tmp_path / "cifs"
    src.mkdir()
    write(src / "1234567.cif", cif(NACL_ATOMS))
    base = tmp_path / "base.hsrdb"
    assert run_cli("add", src, "--hs4", base, "-j", "1").returncode == 0
    write(src / "1234568.cif", cif(NACL_ATOMS, cod=1234568))
    out = tmp_path / "out.hsrdb"
    r = run_cli("add", src, "--hs4", out, "--base", base, "-j", "1")
    assert r.returncode == 0 and "added {'hs4': 1}, already present 1" in r.stdout, r.stdout + r.stderr
    for path, n in ((base, 1), (out, 2)):  # the base is not changed
        db = sqlite3.connect(path)
        assert db.execute("select count(*) from general_indexed").fetchone()[0] == n
        db.close()


def test_cli_internal_error_is_reported_and_exits_1(tmp_path):
    """A correction that names an atom the CIF lacks (COD 1544611 corrects sites M12 and M14) is a bug in the
    corrections, not in the CIF: traceback, INTERNAL in the report, exit 1, and the other CIFs still built."""
    src = tmp_path / "cifs"
    src.mkdir()
    write(src / "1544611.cif", cif(NACL_ATOMS, cod=1544611))
    write(src / "1234567.cif", cif(NACL_ATOMS))
    out = tmp_path / "db.hsrdb"
    r = run_cli("add", src, "--hs4", out, "-j", "1")
    assert r.returncode == 1 and "1 internal errors" in r.stderr
    assert "1544611.cif: INTERNAL ERROR" in r.stdout and "Traceback" in r.stdout
    assert "1544611.cif\tINTERNAL KeyError" in (tmp_path / "db.skipped.tsv").read_text(encoding="utf-8")
    db = sqlite3.connect(out)
    assert db.execute("select ProductID from general_indexed").fetchall() == [("961234568",)]
    db.close()


def test_cli_notes_in_the_report(tmp_path):
    rows = "Fe1 Fe 0.113 0.271 0.349 1 0.011 d\nS1 S 0.362 0.083 0.187 1 0.013 d\nCg1 C 0.25 0.25 0.25 1 0 dum\n"
    text = cif(rows, sg="P 1 21/c 1", a=5.7, b=5.6, c=5.7, be=112).replace(
        "_atom_site_U_iso_or_equiv\n", "_atom_site_U_iso_or_equiv\n_atom_site_calc_flag\n")
    src = write(tmp_path / "1234567.cif", text)
    assert run_cli("add", src, "--hs4", tmp_path / "db.hsrdb", "-j", "1").returncode == 0
    assert f"{src}\tnote: 1 dummy atoms do not scatter" in (tmp_path / "db.skipped.tsv").read_text(encoding="utf-8")


def _sources(tmp_path, n=7):
    """n NaCl CIFs in this order: c0, a CIF with c0's COD id (already present), c1 without a cell (skipped), c2 ...,
    c4 with a row without element (a note), and last another CIF with c0's id: report lines and counts before and
    after the checkpoints."""
    src = tmp_path / "cifs"
    src.mkdir()
    for k in range(n):
        text = cif(NACL_ATOMS + ("Q1 ? 0.1 0.1 0.1 0 0.01\n" if k == 4 else ""), cod=1234560 + k)
        write(src / f"{1234560 + k}.cif", text.replace("_cell_length_a 5.6402\n", "") if k == 1 else text)
    write(src / "1234560a.cif", cif(NACL_ATOMS, cod=1234560))
    write(src / "1234569dup.cif", cif(NACL_ATOMS, cod=1234560))
    return src


def _content(path):
    """Everything in a database except timestamps, and its schema."""
    db = sqlite3.connect(path)
    rows = {t: db.execute(f"select * from {t} order by 1").fetchall()
            for t in ("Phases", "PatternPhase", "general_stored", "literature_map", "SubFileRef", "Properties")}
    rows["general"] = db.execute("select id, ProductID, A, IIC, DX from general_indexed order by id").fetchall()
    schema = db.execute("select type, name, sql from sqlite_master order by name").fetchall()
    db.close()
    return rows, schema


@pytest.mark.parametrize("crash", ["add", "checkpoint"])
def test_cli_interrupted_build_resumes_to_the_same_databases(tmp_path, monkeypatch, capsys, crash):
    """Interrupted while adding c4 (its hs4 entry written, not committed), or in the second checkpoint
    after hs3 committed and before hs4 did (the outputs one checkpoint apart): resumed, the databases and the report
    equal those of an uninterrupted build."""
    src = _sources(tmp_path)
    ref = tmp_path / "ref"
    ref.mkdir()
    cli.main(["add", str(src), "--hs4", str(ref / "db4.hsrdb"), "--hs3", str(ref / "db3.hsrdb"), "-j", "1"])

    monkeypatch.setattr(cli, "CHECKPOINT", 2)
    real, calls = (hsrdb.add if crash == "add" else cli._save_state), []

    def crashing(*a, **k):  # the hs3 entry of c4; the hs4 state of the second checkpoint
        calls.append(1)
        if len(calls) == (8 if crash == "add" else 4):
            raise KeyboardInterrupt
        return real(*a, **k)

    target = (hsrdb, "add") if crash == "add" else (cli, "_save_state")
    monkeypatch.setattr(*target, crashing)
    out = tmp_path / "out"
    out.mkdir()
    args = ["add", str(src), "--hs4", str(out / "db4.hsrdb"), "--hs3", str(out / "db3.hsrdb"), "-j", "1"]
    with pytest.raises(KeyboardInterrupt):
        cli.main(args)
    assert "continue from the last checkpoint with the same command and --resume" in capsys.readouterr().err
    assert sorted(p.name for p in out.iterdir()) == ["db3.hsrdb.partial", "db4.hsrdb.partial", "db4.skipped.tsv"]
    for name, done, n in (("db4.hsrdb.partial", 4, 2), ("db3.hsrdb.partial", 4, 2)) if crash == "add" else (
            ("db4.hsrdb.partial", 2, 1), ("db3.hsrdb.partial", 4, 2)):
        db = sqlite3.connect(out / name)
        assert cli._state(db)["done"] == done and db.execute("select count(*) from general_indexed").fetchone()[0] == n
        db.close()

    monkeypatch.setattr(*target, real)  # CHECKPOINT stays 2: the resumed build checkpoints too
    with pytest.raises(SystemExit, match="--resume"):  # a new build would discard the checkpoint
        cli.main(args)
    with pytest.raises(SystemExit, match="same CIFs"):
        cli.main(["add", str(src / "1234560.cif"), *args[2:], "--resume"])
    with pytest.raises(SystemExit, match="nothing to resume"):  # a typing error in the second output
        cli.main([*args[:5], str(out / "typo.hsrdb"), *args[6:], "--resume"])
    real_packages = cli._packages
    monkeypatch.setattr(cli, "_packages", lambda: real_packages() | {"numpy": "0.0"})
    with pytest.raises(SystemExit, match="numpy and spglib versions"):  # other numerics, another pattern
        cli.main(args + ["--resume"])
    monkeypatch.setattr(cli, "_packages", real_packages)
    changed = src / "1234566.cif"
    before, stat = changed.read_bytes(), changed.stat()
    changed.write_bytes(before + b"# revised\n")  # a mirror updated in between
    with pytest.raises(SystemExit, match="same CIFs, unchanged"):
        cli.main(args + ["--resume"])
    changed.write_bytes(before)
    os.utime(changed, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    real_add = hsrdb.add

    def crash_now(*a, **k):
        raise KeyboardInterrupt

    monkeypatch.setattr(hsrdb, "add", crash_now)
    with pytest.raises(KeyboardInterrupt):  # a resume interrupted before its own first checkpoint
        cli.main(args + ["--resume"])
    monkeypatch.setattr(hsrdb, "add", real_add)
    assert (out / "db4.hsrdb.partial").exists() and (out / "db3.hsrdb.partial").exists()  # kept by all of these
    cli.main(args[:-1] + ["2", "--resume"])  # -j may differ
    assert sorted(p.name for p in out.iterdir()) == ["db3.hsrdb", "db4.build.json", "db4.hsrdb", "db4.skipped.tsv"]
    for name in ("db4.hsrdb", "db3.hsrdb"):
        assert _content(out / name) == _content(ref / name)
        db = sqlite3.connect(out / name)
        assert not private_tables(db)
        db.close()
    assert (out / "db4.skipped.tsv").read_text(encoding="utf-8") == (ref / "db4.skipped.tsv").read_text(encoding="utf-8")
    manifest = json.loads((out / "db4.build.json").read_text(encoding="utf-8"))
    assert manifest["resumes"] == 1 and manifest["already_present"] == 2 and manifest["skipped"] == 1
    for t in ("hs4", "hs3"):
        assert manifest["outputs"][t] == {"path": str(out / f"db{t[2]}.hsrdb"), "patterns": 6, "added": 6}


def test_cli_interruption_before_the_first_checkpoint_leaves_nothing(tmp_path, monkeypatch):
    src = _sources(tmp_path, 3)

    def crash(*a, **k):
        raise RuntimeError("simulated")

    monkeypatch.setattr(hsrdb, "add", crash)
    with pytest.raises(RuntimeError, match="simulated"):
        cli.main(["add", str(src), "--hs4", str(tmp_path / "db.hsrdb"), "-j", "1"])
    assert not list(tmp_path.glob("db.hsrdb*"))


def test_cli_keeps_a_complete_database_that_cannot_be_compacted(tmp_path, monkeypatch, capsys):
    """VACUUM needs up to twice the database's size free: if it fails (disk full), the complete database is still
    renamed into place, the failure is reported and recorded, and the exit status is 1."""
    real = sqlite3.connect

    class Connection(sqlite3.Connection):
        def execute(self, sql, *a):
            if sql.strip().lower() == "vacuum":
                raise sqlite3.OperationalError("database or disk is full")
            return super().execute(sql, *a)

    monkeypatch.setattr(cli.sqlite3, "connect", functools.partial(real, factory=Connection))
    src = write(tmp_path / "1234567.cif", cif(NACL_ATOMS))
    out = tmp_path / "db.hsrdb"
    with pytest.raises(SystemExit) as exit_:
        cli.main(["add", str(src), "--hs4", str(out), "-j", "1"])
    assert exit_.value.code == 1 and "complete but not compacted (database or disk is full)" in capsys.readouterr().err
    db = real(out)
    assert db.execute("select count(*) from general_indexed").fetchone()[0] == 1
    db.close()
    manifest = json.loads((tmp_path / "db.build.json").read_text(encoding="utf-8"))
    assert manifest["not_compacted"] == {"hs4": "database or disk is full"}
    assert manifest["command"] == ["add", str(src), "--hs4", str(out), "-j", "1"]


def test_cli_resume_refuses_what_it_cannot_continue(tmp_path, monkeypatch):
    src = _sources(tmp_path, 3)
    out = tmp_path / "db.hsrdb"
    args = ["add", str(src), "--hs4", str(out), "-j", "1"]
    with pytest.raises(SystemExit, match="nothing to resume"):
        cli.main(args + ["--resume"])
    hsrdb.create(out.with_name("db.hsrdb.partial")).close()  # a database without checkpoint
    with pytest.raises(SystemExit, match="holds no interrupted build"):
        cli.main(args + ["--resume"])
    monkeypatch.setattr(cli, "CHECKPOINT", 1)
    real = hsrdb.add

    def crash(*a, **k):
        raise KeyboardInterrupt

    monkeypatch.setattr(hsrdb, "add", lambda *a, **k: (monkeypatch.setattr(hsrdb, "add", crash), real(*a, **k))[1])
    with pytest.raises(KeyboardInterrupt):  # the first CIF is added and checkpointed, the second crashes
        cli.main(args)
    monkeypatch.setattr(hsrdb, "add", real)
    report = tmp_path / "db.skipped.tsv"
    report.write_text("file\tproblem\n", encoding="utf-8")  # shorter than at the checkpoint
    with pytest.raises(SystemExit, match="report of the interrupted build is missing or shorter"):
        cli.main(args + ["--resume"])
    report.unlink()
    with pytest.raises(SystemExit, match="report of the interrupted build is missing or shorter"):
        cli.main(args + ["--resume"])
    assert out.with_name("db.hsrdb.partial").exists()  # kept through all the refusals


@pytest.mark.parametrize("content", [b"not a database", b""])
def test_cli_replaces_a_partial_file_that_is_no_interrupted_build(tmp_path, content):
    """A .partial without checkpoint (an empty database, or another file) is no build to resume: replaced."""
    src = write(tmp_path / "1234567.cif", cif(NACL_ATOMS))
    out = tmp_path / "db.hsrdb"
    out.with_name("db.hsrdb.partial").write_bytes(content)
    cli.main(["add", str(src), "--hs4", str(out), "-j", "1"])
    db = sqlite3.connect(out)
    assert db.execute("select count(*) from general_indexed").fetchone()[0] == 1
    db.close()


def test_cli_keeps_a_partial_database_it_cannot_read(tmp_path):
    """A damaged .partial may hold a checkpoint: a new build stops and leaves it, rather than delete it."""
    src = write(tmp_path / "1234567.cif", cif(NACL_ATOMS))
    out = tmp_path / "db.hsrdb"
    partial = out.with_name("db.hsrdb.partial")
    partial.write_bytes(b"SQLite format 3\0" + bytes(84) + b"damaged" * 100)
    with pytest.raises(SystemExit, match="cannot be read"):
        cli.main(["add", str(src), "--hs4", str(out), "-j", "1"])
    assert partial.exists() and not out.exists()


def test_cli_failure_while_finishing_is_not_offered_as_resumable(tmp_path, monkeypatch, capsys):
    """Interrupted after one output dropped its checkpoint and before the other did: the message says the kept
    files cannot be resumed."""
    src = _sources(tmp_path, 3)
    monkeypatch.setattr(cli, "CHECKPOINT", 1)
    real, calls = hsrdb.update_prefixes, []

    def crash_on_second(db):
        calls.append(1)
        if len(calls) == 2:
            raise KeyboardInterrupt
        return real(db)

    monkeypatch.setattr(hsrdb, "update_prefixes", crash_on_second)
    with pytest.raises(KeyboardInterrupt):
        cli.main(["add", str(src), "--hs4", str(tmp_path / "db4.hsrdb"), "--hs3", str(tmp_path / "db3.hsrdb"),
                  "-j", "1"])
    assert "interrupted while finishing" in capsys.readouterr().err


def test_cli_resumed_build_keeps_the_internal_errors_of_its_first_part(tmp_path, monkeypatch):
    """An internal error before the interruption still fails the resumed build: exit 1, counted in the manifest."""
    src = _sources(tmp_path, 3)
    write(src / "0.cif", cif(NACL_ATOMS, cod=1544611))  # first; the correction for COD 1544611 names absent sites
    monkeypatch.setattr(cli, "CHECKPOINT", 1)
    real, calls = hsrdb.add, []

    def crash_on_second(*a, **k):
        calls.append(1)
        if len(calls) == 2:
            raise KeyboardInterrupt
        return real(*a, **k)

    monkeypatch.setattr(hsrdb, "add", crash_on_second)
    args = ["add", str(src), "--hs4", str(tmp_path / "db.hsrdb"), "-j", "1"]
    with pytest.raises(KeyboardInterrupt):
        cli.main(args)
    monkeypatch.setattr(hsrdb, "add", real)
    with pytest.raises(SystemExit, match="1 internal errors"):
        cli.main(args + ["--resume"])
    manifest = json.loads((tmp_path / "db.build.json").read_text(encoding="utf-8"))
    assert manifest["internal_errors"] == 1 and manifest["resumes"] == 1


def test_cli_report_that_cannot_be_opened_leaves_no_partial_file(tmp_path):
    src = write(tmp_path / "1234567.cif", cif(NACL_ATOMS))
    (tmp_path / "db.skipped.tsv").mkdir()
    with pytest.raises(OSError):
        cli.main(["add", str(src), "--hs4", str(tmp_path / "db.hsrdb"), "-j", "1"])
    assert not list(tmp_path.glob("db.hsrdb*"))


def test_cli_manifest_that_cannot_be_written(tmp_path, capsys):
    """The databases are complete: kept, the manifest's content printed instead, exit status 1."""
    src = write(tmp_path / "1234567.cif", cif(NACL_ATOMS))
    (tmp_path / "db.build.json").mkdir()
    with pytest.raises(SystemExit) as exit_:
        cli.main(["add", str(src), "--hs4", str(tmp_path / "db.hsrdb"), "-j", "1"])
    err = capsys.readouterr().err
    assert exit_.value.code == 1 and "db.build.json: cannot be written" in err and '"source_sha256"' in err
    db = sqlite3.connect(tmp_path / "db.hsrdb")
    assert db.execute("select count(*) from general_indexed").fetchone()[0] == 1
    db.close()


@pytest.mark.skipif(sys.platform == "win32", reason="a shell script stands in for git")
def test_provenance_with_git_output_that_is_not_utf8(tmp_path, monkeypatch):
    """A checkout path that is not UTF-8 (Latin-1 on Linux) must not fail a finished build."""
    git = tmp_path / "git"
    git.write_bytes(b"#!/bin/sh\nprintf '/home/caf\\351/repo\\n'\n")
    git.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    assert cli._provenance()["git"] is None  # not the package's own checkout
