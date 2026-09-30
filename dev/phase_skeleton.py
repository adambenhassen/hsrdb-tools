"""Print the skeleton of a decompressed phase record: refinable values (TRV) by name, gaps as hex."""

import sqlite3
import struct
import sys
import zlib

TRV_FIXED = 4 + 1 + 7 * 8 + 7  # int32 tag, u8, 7 doubles, 7 bytes, then a length-prefixed UTF-16 name


def read_str(d, i):
    n = struct.unpack_from("<i", d, i)[0]
    return d[i + 4:i + 4 + 2 * n].decode("utf-16le"), i + 4 + 2 * n


def trv_at(d, i, prefix):
    if i + TRV_FIXED + 4 > len(d) or struct.unpack_from("<i", d, i)[0] != 2:
        return None
    j = i + TRV_FIXED
    n = struct.unpack_from("<i", d, j)[0]
    if not 0 < n < 300 or j + 4 + 2 * n > len(d):
        return None
    try:
        name = d[j + 4:j + 4 + 2 * n].decode("utf-16le")
    except UnicodeDecodeError:
        return None
    if not name.startswith(prefix):
        return None
    vals = struct.unpack_from("<7d", d, i + 5)
    flags = d[i + 4], d[i + 61:i + 68]
    return name[len(prefix):], vals, flags, j + 4 + 2 * n


def skeleton(d, prefix):
    out = []
    i = 0
    gap = bytearray()
    gap_start = 0
    while i < len(d):
        t = trv_at(d, i, prefix)
        if t:
            if gap:
                out.append(("gap", gap_start, bytes(gap)))
                gap = bytearray()
            name, vals, flags, i2 = t
            out.append(("trv", i, name, vals, flags))
            i = i2
            gap_start = i
        else:
            gap.append(d[i])
            i += 1
    if gap:
        out.append(("gap", gap_start, bytes(gap)))
    return out


def load(db, pid):
    name, blob = db.execute("select Name, Phase from Phases where ID=?", (pid,)).fetchone()
    return name, zlib.decompress(blob[6:], -15)


if __name__ == "__main__":
    db = sqlite3.connect("file:COD24_HS4x.hsrdb?mode=ro", uri=True)
    pid = int(sys.argv[1])
    name, d = load(db, pid)
    code = name.split("96-")[1].replace("-", "")
    prefix = str(int(code) - 1) + " "  # descriptions start with the phase name (the COD id for unnamed ones)
    first = trv_at(d, 4, "")
    if first:
        prefix = first[0].rsplit("Scale Factor", 1)[0]
    print(name, len(d), "prefix", repr(prefix))
    for item in skeleton(d, prefix):
        if item[0] == "trv":
            _, off, nm, vals, (b, fl) = item
            print(f"{off:6d} TRV {nm!r:40} b={b} v={vals} f={fl.hex()}")
        else:
            _, off, g = item
            print(f"{off:6d} GAP[{len(g)}] {g.hex()}")
