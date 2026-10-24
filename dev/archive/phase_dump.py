"""Heuristic token dump of a decompressed phase record: UTF-16 strings, doubles and raw bytes, with offsets."""

import struct
import sys


def is_str(d, i):
    if i + 4 > len(d):
        return 0
    n = struct.unpack_from("<i", d, i)[0]
    if not 0 < n < 400 or i + 4 + 2 * n > len(d):
        return 0
    s = d[i + 4:i + 4 + 2 * n]
    try:
        t = s.decode("utf-16le")
    except UnicodeDecodeError:
        return 0
    return n if all(c.isprintable() for c in t) else 0


def plausible_double(d, i):
    if i + 8 > len(d):
        return None
    v = struct.unpack_from("<d", d, i)[0]
    if v != v or v == 0:
        return None
    return v if 1e-6 < abs(v) < 1e7 else None


def dump(d, start=0, end=None):
    i = start
    end = end or len(d)
    raw = []

    def flush():
        if raw:
            hexs = ' '.join(f'{b:02x}' for _, b in raw)
            while ' 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00' in hexs:
                hexs = hexs.replace(' 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00', ' 0x16', 1)
            print(f"  raw {raw[0][0]:6d}: {hexs}")
            raw.clear()

    while i < end:
        n = is_str(d, i)
        if n:
            flush()
            print(f"{i:6d} str[{n}] {d[i + 4:i + 4 + 2 * n].decode('utf-16le')!r}")
            i += 4 + 2 * n
            continue
        v = plausible_double(d, i)
        if v is not None:
            flush()
            print(f"{i:6d} f64 {v!r}")
            i += 8
            continue
        raw.append((i, d[i]))
        if len(raw) == 64:
            flush()
        i += 1
    flush()


if __name__ == "__main__":
    data = open(sys.argv[1], "rb").read()
    a = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    b = int(sys.argv[3]) if len(sys.argv) > 3 else None
    dump(data, a, b)
