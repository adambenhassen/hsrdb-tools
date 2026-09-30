"""Create and fill HighScore reference databases (.hsrdb, SQLite) for 3.x and 4.x+; HighScore calls both "V3.X"."""

import math
import sqlite3
import struct
import time
from pathlib import Path

import numpy as np

from . import phase, structure

# target -> (phase record format, schema); HighScore 3.x databases have no LinesI column
TARGETS = {
    "hs4": (phase.HS4, Path(__file__).with_name("schema.sql")),
    "hs3": (phase.HS3, Path(__file__).with_name("schema_hs3.sql")),
}
SUBFILES = ["User Inorganic", "User Organic", "User Mineral", "User Metallic", "User Amorphous", "User Ceramic",
            "User Phases Mixture", "User First Derivative", "User Second Derivative", "User Spectrum",
            "User Reflectance", "User Log(1/R)"]


def create(path, target="hs4", db_id="CIFDB"):
    """New empty database with the same layout as the official COD .hsrdb files for the target."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    db = sqlite3.connect(path)
    db.execute("pragma page_size = 1024")
    db.executescript(TARGETS[target][1].read_text())
    db.executemany("insert into Properties(Key, Value) values (?, ?)",
                   [("NoOfPatterns", "0"), ("ReferenceCodePrefixes", ""), ("ID", db_id)])
    db.executemany("insert into Subfiles(ID, SubFile) values (?, ?)", list(enumerate(SUBFILES)))
    db.commit()
    return db


def target_of(db):
    """'hs4' or 'hs3', from the schema of an existing database."""
    return "hs4" if "LinesI" in [r[1] for r in db.execute("pragma table_info(general_stored)")] else "hs3"


def encode_lines(lines):
    """Lines, HKL and LinesI blobs: float32 d + uint16 intensity, int16 hkl, float32 intensity; strongest first."""
    for x in lines:
        if not (math.isfinite(x.d) and x.d > 0 and math.isfinite(x.intensity) and 0 <= x.intensity <= 1000.5):
            raise ValueError(f"invalid line d={x.d} I={x.intensity}")
        if max(map(abs, x.hkl)) > 32767:
            raise ValueError(f"hkl {x.hkl} outside int16")
    lines = sorted(lines, key=lambda x: -x.intensity)
    L = b"".join(struct.pack("<fH", x.d, int(round(x.intensity))) for x in lines)
    H = b"".join(struct.pack("<hhh", *x.hkl) for x in lines)
    I = struct.pack(f"<{len(lines)}f", *(x.intensity for x in lines))
    return L, H, I


def _f32(v):
    return float(np.float32(v))


def add(db, entry, masks, line_blobs, phase_blob, target="hs4"):
    """Insert one entry with its encode_lines() blobs and phase record (from phase_blob()). Returns its row id."""
    fmt = TARGETS[target][0]
    pid = (db.execute("select coalesce(max(id), 0) from general_indexed").fetchone()[0]) + 1
    doc = (db.execute("select coalesce(max(doc), 0) from literature_map").fetchone()[0]) + 1
    now = int(time.time())
    a, b, c, al, be, ga = (_f32(v) for v in entry.cell)
    L, H, I = line_blobs
    e1, e2 = masks
    lines_i = ", LinesI" if fmt is phase.HS4 else ""
    db.execute(
        f"""insert into general (id, ProductID, XTSLSYS, STATUS, A, B, C, ALPHA, BETA, GAMMA, Z, DM, DX, QUALFINAL,
               IIC, NumberOfElements, XTLASPECTN, XTLVOL, CTIME, MTIME, CELLED, SPGR, SPGRED, ANX, XTLSG, XTLSGED,
               Comment, CAS, Lines, HKL{lines_i}, chemicalformula, compoundname, mineralname, commonname,
               empiricalformula, color, e1, e2)
           values (?, ?, ?, 0, ?, ?, ?, ?, ?, ?, 0, 0.0, ?, '=', ?, ?, ?, ?, ?, ?, ' ', ?, NULL, '', NULL, NULL, ?,
                   '', ?, ?{", ?" if lines_i else ""}, ?, ?, ?, ?, '', '', ?, ?)""",
        (pid, str(entry.reference_code), entry.crystal_system, a, b, c, al, be, ga, entry.density,
         _f32(entry.iic), len(entry.elements), entry.sg_number, entry.volume, now, now, entry.sg_hm,
         entry.comment, L, H, *((I,) if lines_i else ()), entry.formula, entry.compound_name, entry.mineral_name,
         entry.compound_name, e1, e2))
    lit = entry.literature
    db.execute("insert into literature(doc, pid, volume, pages, year, authors, journal, referencetype) "
               "values (?, ?, ?, ?, ?, ?, ?, 2)",
               (doc, pid, lit["volume"], lit["pages"], lit["year"], lit["authors"], lit["journal"]))
    db.executemany("insert into SubFileRef(pid, SubFileID) values (?, ?)", [(pid, s) for s in entry.subfiles])
    code = str(entry.reference_code)
    name = f"PDF:{code[:2]}-{code[2:5]}-{code[5:]}"
    db.execute("insert into Phases(Name, Phase) values (?, ?)", (name, phase_blob))
    db.execute("insert into PatternPhase(ReferenceCode, PhaseName) values (?, ?)", (code, name))
    return pid


def phase_blob(entry, target="hs4", sites=None):
    """Compressed phase record of an entry for the target HighScore version.

    sites is structure.site_info(entry.structure, entry.images), to share between several targets.
    """
    return phase.encode_blob(phase.encode_payload(structure.build(entry, TARGETS[target][0], sites)))


def update_prefixes(db):
    """Add every '96-xxx' prefix present to ReferenceCodePrefixes (newline-terminated), keeping existing ones."""
    prefixes = {f"{code[:2]}-{code[2:5]}" for (code,) in db.execute("select ProductID from general_indexed")}
    old = db.execute("select Value from Properties where Key='ReferenceCodePrefixes'").fetchone()[0] or ""
    known = [p for p in old.split("\n") if p]
    merged = known + sorted(prefixes - set(known))
    db.execute("update Properties set Value=? where Key='ReferenceCodePrefixes'", ("\n".join(merged) + "\n",))
