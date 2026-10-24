"""List ANSI (int32-length) strings in a version-3 phase record, with the bytes between them."""
import struct, sys, re
d = open(sys.argv[1], "rb").read()
i = 0; last = 0; out = []
while i < len(d) - 4:
    n = struct.unpack_from("<i", d, i)[0]
    if 2 <= n < 300 and i + 4 + n <= len(d):
        s = d[i + 4:i + 4 + n]
        if all(32 <= c < 127 for c in s) and re.search(rb"[A-Za-z]{2}", s):
            out.append((last, d[last:i].hex(), s.decode()))
            i += 4 + n; last = i; continue
    i += 1
out.append((last, d[last:].hex(), None))
lim = int(sys.argv[2]) if len(sys.argv) > 2 else 60
for off, gap, s in out[:lim]:
    print(f"{off:6d} gap[{len(gap)//2}] {gap[:260]}  ->  {s!r}")
print("...", len(out), "strings")
