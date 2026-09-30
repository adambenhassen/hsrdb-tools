"""HighScore phase (crystal structure) records, as stored in the Phases table of a .hsrdb.

Blob layout: uint16 8, uint32 CRC-32 of the payload, raw-deflated payload.
The payload is a fixed sequence of refinable values (TRV: value, esd, limits and a description), atom
blocks and a bibliographic trailer; every field is little-endian and strings are int32-length prefixed.

Two record versions exist: 3 (HighScore 3.x, single-byte strings) and 6 (HighScore 4.x and later, UTF-16).
"""

import struct
import unicodedata
import zlib
from dataclasses import dataclass, field

BLOB_VERSION = 8
CRYSTAL_SYSTEM_CODES = {"triclinic": 1, "monoclinic": 3, "orthorhombic": 5, "tetragonal": 6,
                        "hexagonal": 8, "cubic": 9, "trigonal": 10}
CELL_NAMES = ["Cell a [Å]", "Cell b [Å]", "Cell c [Å]", "Cell alpha [°]", "Cell beta [°]", "Cell gamma [°]"]
ATOM_POS = ["X", "Y", "Z"]
ATOM_ANISO = ["B11", "B22", "B33", "B12", "B13", "B23"]


@dataclass(frozen=True)
class Format:
    version: int
    encoding: str
    char_size: int
    trv_tag: int
    after_scale: int
    header_has_minus_one: bool
    atom_marker: int
    source: int


HS4 = Format(version=6, encoding="utf-16le", char_size=2, trv_tag=2, after_scale=3, header_has_minus_one=True,
             atom_marker=4, source=3)
HS3 = Format(version=3, encoding="cp1252", char_size=1, trv_tag=1, after_scale=2, header_has_minus_one=False,
             atom_marker=2, source=2)
FORMATS = {f.version: f for f in (HS3, HS4)}


@dataclass
class TRV:
    """One refinable value. values = (value, esd, last shift, min, max, max - min, 0)."""
    name: str
    values: tuple
    refine: int = 1
    flags: bytes = b"\x01\x00\x00\x00\x00\x00\x00"

    @property
    def value(self):
        return self.values[0]


@dataclass
class Atom:
    label: str
    element: str
    pos: list  # 3 TRV
    biso: TRV
    sof: TRV
    aniso: list  # 6 TRV
    multiplicity: int
    wyckoff: str
    charge: int = 0
    tail: bytes = bytes(26)  # after the aniso values; byte 20 is the calc flag


@dataclass
class Phase:
    prefix: str  # descriptions are "<prefix><parameter>"; prefix is the phase name plus a space
    scale: TRV
    cell: list  # 6 TRV
    crystal_system: int
    atoms: list
    after_atoms: bytes  # constant block between the last atom and the profile parameters
    profile: list  # TRV up to 'TOF Beta2'
    mid: bytes
    corrections: list  # TRV from preferred orientation to 'Roughness'
    harmonics_head: bytes  # preferred-orientation setup incl. the coefficient index list
    harmonics: list  # TRV 'Coefficient ...'
    trailer: dict = field(default_factory=dict)
    fmt: Format = HS4
    cell_flag: int = 1  # byte after the cell; 1 in almost all official records


class _R:
    def __init__(self, data, fmt):
        self.d = data
        self.i = 0
        self.fmt = fmt

    def i32(self):
        v = struct.unpack_from("<i", self.d, self.i)[0]
        self.i += 4
        return v

    def u8(self):
        v = self.d[self.i]
        self.i += 1
        return v

    def f64(self):
        v = struct.unpack_from("<d", self.d, self.i)[0]
        self.i += 8
        return v

    def raw(self, n):
        v = self.d[self.i:self.i + n]
        self.i += n
        return v

    def s(self):
        n = self.i32()
        return self.raw(self.fmt.char_size * n).decode(self.fmt.encoding)

    def trv(self, prefix):
        if self.i32() != self.fmt.trv_tag:
            raise ValueError(f"not a refinable value at {self.i - 4}")
        refine = self.u8()
        values = struct.unpack_from("<7d", self.d, self.i)
        self.i += 56
        flags = self.raw(7)
        name = self.s()
        if not name.startswith(prefix):
            raise ValueError(f"unexpected TRV name {name!r} at {self.i}")
        return TRV(name[len(prefix):], values, refine, flags)

    def peek_trv_name(self, prefix):
        save = self.i
        try:
            return self.trv(prefix).name
        except (ValueError, UnicodeDecodeError, struct.error):
            return None
        finally:
            self.i = save


class _W:
    def __init__(self, fmt):
        self.b = bytearray()
        self.fmt = fmt

    def i32(self, v):
        self.b += struct.pack("<i", v)

    def u8(self, v):
        self.b.append(v)

    def f64(self, v):
        self.b += struct.pack("<d", v)

    def raw(self, v):
        self.b += v

    def s(self, v):
        if self.fmt.char_size == 1:
            try:
                data = v.encode(self.fmt.encoding)
            except UnicodeEncodeError:
                # single-byte (HighScore 3.x) records: drop accents first, then mark what is left with '?'
                v = "".join(c for c in unicodedata.normalize("NFKD", v) if not unicodedata.combining(c))
                data = v.encode(self.fmt.encoding, errors="replace")
        else:
            data = v.encode(self.fmt.encoding)
        self.i32(len(data) // self.fmt.char_size)
        self.b += data

    def trv(self, t, prefix):
        if len(t.flags) != 7:
            raise ValueError(f"TRV {t.name!r} needs 7 flag bytes")
        self.i32(self.fmt.trv_tag)
        self.u8(t.refine)
        self.b += struct.pack("<7d", *t.values)
        self.b += t.flags
        self.s(prefix + t.name)


TRAILER_STRINGS_1 = ["date", "method", "s1", "s2", "systematic", "mineral", "common", "formula_sum",
                     "formula_struct", "title", "authors", "s3", "journal", "coden"]
TRAILER_INTS = ["volume", "year", "page_first", "page_last"]


def decode_payload(d):
    fmt = FORMATS.get(struct.unpack_from("<i", d, 0)[0])
    if fmt is None:
        raise ValueError("unknown record version")
    r = _R(d, fmt)
    r.i = 4
    # the phase name prefix comes from the scale factor's description
    scale_name = r.trv("").name
    prefix = scale_name[: -len("Scale Factor")]
    r.i = 4
    scale = r.trv(prefix)
    if r.i32() != fmt.after_scale:
        raise ValueError("unexpected value after scale")
    cell = [r.trv(prefix) for _ in range(6)]
    cell_flag = r.u8()
    cs = r.i32()
    if r.raw(25) != bytes(25) or (fmt.header_has_minus_one and r.f64() != -1.0) or r.i32() != 2:
        raise ValueError("unexpected cell trailer")
    natoms = r.i32()
    atoms = []
    for _ in range(natoms):
        if r.i32() != fmt.atom_marker:
            raise ValueError("expected atom marker")
        pos = [r.trv(prefix) for _ in range(3)]
        if r.i32() != 0:
            raise ValueError("expected 0 after position")
        biso, sof = r.trv(prefix), r.trv(prefix)
        label, element = r.s(), r.s()
        charge, mult = r.i32(), r.i32()
        wyckoff = r.s()
        aniso = [r.trv(prefix) for _ in range(6)]
        tail = r.raw(26)
        atoms.append(Atom(label=label, element=element, pos=pos, biso=biso, sof=sof, aniso=aniso,
                          multiplicity=mult, wyckoff=wyckoff, charge=charge, tail=tail))
    start = r.i
    while r.peek_trv_name(prefix) is None:
        r.i += 1
    after_atoms = d[start:r.i]
    profile = []
    while True:
        profile.append(r.trv(prefix))
        if profile[-1].name == "TOF Beta2":
            break
    start = r.i
    while r.peek_trv_name(prefix) is None:
        r.i += 1
    mid = d[start:r.i]
    corrections = []
    while True:
        corrections.append(r.trv(prefix))
        if corrections[-1].name == "Roughness":
            break
    start = r.i
    while r.peek_trv_name(prefix) is None and not _trailer_follows(d, r.i, prefix, fmt):
        r.i += 1
    harmonics_head = d[start:r.i]
    harmonics = []
    while r.peek_trv_name(prefix) is not None:
        harmonics.append(r.trv(prefix))
    t = {"pre": r.raw(16), "r_factor": r.f64(), "name": r.s(), "source": r.i32(), "cod_id": r.i32(), "i0": r.i32()}
    for k in TRAILER_STRINGS_1:
        t[k] = r.s()
    for k in TRAILER_INTS:
        t[k] = r.i32()
    t["s4"] = r.s()
    t["comment"] = r.s()
    t["rest"] = r.raw(len(d) - r.i)
    return Phase(prefix=prefix, scale=scale, cell=cell, crystal_system=cs, atoms=atoms, after_atoms=after_atoms,
                 profile=profile, mid=mid, corrections=corrections, harmonics_head=harmonics_head,
                 harmonics=harmonics, trailer=t, fmt=fmt, cell_flag=cell_flag)


def _trailer_follows(d, i, prefix, fmt):
    """The trailer starts with 16 zero bytes, a double and the phase name."""
    if d[i:i + 16] != bytes(16) or i + 28 > len(d):
        return False
    n = struct.unpack_from("<i", d, i + 24)[0]
    try:
        return d[i + 28:i + 28 + fmt.char_size * n].decode(fmt.encoding) + " " == prefix
    except UnicodeDecodeError:
        return False


def encode_payload(p):
    if len(p.cell) != 6 or any(len(a.pos) != 3 or len(a.aniso) != 6 or len(a.tail) != 26 for a in p.atoms):
        raise ValueError("phase record has wrongly sized cell or atom fields")
    fmt = p.fmt
    w = _W(fmt)
    w.i32(fmt.version)
    w.trv(p.scale, p.prefix)
    w.i32(fmt.after_scale)
    for t in p.cell:
        w.trv(t, p.prefix)
    w.u8(p.cell_flag)
    w.i32(p.crystal_system)
    w.raw(bytes(25))
    if fmt.header_has_minus_one:
        w.f64(-1.0)
    w.i32(2)
    w.i32(len(p.atoms))
    for a in p.atoms:
        w.i32(fmt.atom_marker)
        for t in a.pos:
            w.trv(t, p.prefix)
        w.i32(0)
        w.trv(a.biso, p.prefix)
        w.trv(a.sof, p.prefix)
        w.s(a.label)
        w.s(a.element)
        w.i32(a.charge)
        w.i32(a.multiplicity)
        w.s(a.wyckoff)
        for t in a.aniso:
            w.trv(t, p.prefix)
        w.raw(a.tail)
    w.raw(p.after_atoms)
    for t in p.profile:
        w.trv(t, p.prefix)
    w.raw(p.mid)
    for t in p.corrections:
        w.trv(t, p.prefix)
    w.raw(p.harmonics_head)
    for t in p.harmonics:
        w.trv(t, p.prefix)
    t = p.trailer
    w.raw(t["pre"])
    w.f64(t["r_factor"])
    w.s(t["name"])
    w.i32(t["source"])
    w.i32(t["cod_id"])
    w.i32(t["i0"])
    for k in TRAILER_STRINGS_1:
        w.s(t[k])
    for k in TRAILER_INTS:
        w.i32(t[k])
    w.s(t["s4"])
    w.s(t["comment"])
    w.raw(t["rest"])
    return bytes(w.b)


def decode_blob(blob):
    version, crc = struct.unpack_from("<HI", blob, 0)
    if version != BLOB_VERSION:
        raise ValueError(f"unknown phase blob version {version}")
    payload = zlib.decompress(blob[6:], -15)
    if zlib.crc32(payload) != crc:
        raise ValueError("phase blob CRC mismatch")
    return payload


def encode_blob(payload):
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    return struct.pack("<HI", BLOB_VERSION, zlib.crc32(payload)) + c.compress(payload) + c.flush()
