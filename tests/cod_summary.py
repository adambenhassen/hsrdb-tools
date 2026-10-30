"""Summary of one real COD CIF, and its comparison with the stored one (tests/data/cod/expected.json).

Shared by tests/test_cod_regression.py and dev/update_expected.py.
"""

import dataclasses
import math

from hsrdb_tools import entry, hsrdb, phase, structure

REL = 1e-9  # floats: cell, d, density, I/Ic, volume, structure record values
PHASE_ABS = 1e-9  # structure record values near zero
LINE_ABS = 1e-6  # line intensities, on the 0-1000 scale
TARGETS = ("hs4", "hs3")
NOTES = ("correction", "symmetry_note", "merge_note", "modulation", "setting_changed", "dropped_sites",
         "invalid_adps", "dummy_sites", "symbol_conflict", "stated_general", "separate_stated", "merged_images")


def _plain(x):
    """JSON-ready copy: dataclasses as dicts, refinable values as their number lists, bytes as hex."""
    if isinstance(x, phase.TRV):
        return list(x.values)  # name and flags come from the templates, not from the entry
    if dataclasses.is_dataclass(x):
        return {f.name: _plain(getattr(x, f.name)) for f in dataclasses.fields(x) if f.name != "fmt"}
    if isinstance(x, dict):
        return {k: _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, bytes):
        return x.hex()
    return x


def _record(e, target, sites):
    """The decoded structure record for the target; the parts that are constant templates are left out."""
    p = phase.decode_payload(phase.decode_blob(hsrdb.phase_blob(e, target, sites)))
    out = _plain(p)
    for k in ("after_atoms", "profile", "mid", "corrections", "harmonics_head", "harmonics"):
        del out[k]
    return out


def summarize(path):
    """Expected-output record of one CIF: {"skip": reason} or the entry's fields, lines and structure records."""
    try:
        e = entry.from_cif(path)
    except entry.EntryError as exc:
        return {"skip": str(exc)}
    sites = structure.site_info(e.structure, e.images)
    out = {"sg_hm": e.sg_hm, "sg_number": e.sg_number, "crystal_system": e.crystal_system, "cell": list(e.cell),
           "volume": e.volume, "density": e.density, "iic": e.iic, "formula": e.formula,
           "elements": list(e.elements), "subfiles": list(e.subfiles), "comment": e.comment,
           "compound_name": e.compound_name, "mineral_name": e.mineral_name,
           "notes": {k: getattr(e, k) for k in NOTES},
           "lines": [[x.d, x.intensity, *x.hkl] for x in e.lines]}
    out["phase"] = {t: _record(e, t, sites) for t in TARGETS}
    return out


def _tolerance(path):
    if path[0] == "lines" and path[-1] == 1:
        return 0.0, LINE_ABS
    if path[0] == "phase":
        return REL, PHASE_ABS
    return REL, 0.0


def _walk(exp, act, path, found):
    if isinstance(exp, dict) and isinstance(act, dict):
        for k in sorted(set(exp) | set(act)):
            if k not in exp or k not in act:
                found.append((path + (k,), exp.get(k, "<missing>"), act.get(k, "<missing>")))
            else:
                _walk(exp[k], act[k], path + (k,), found)
    elif isinstance(exp, list) and isinstance(act, list):
        if len(exp) != len(act):
            found.append((path + ("<length>",), len(exp), len(act)))
        for i, (a, b) in enumerate(zip(exp, act)):
            _walk(a, b, path + (i,), found)
    elif isinstance(exp, float) and isinstance(act, float):
        rel, ab = _tolerance(path)
        if not math.isclose(exp, act, rel_tol=rel, abs_tol=ab):
            found.append((path, exp, act))
    elif type(exp) is not type(act) or exp != act:
        found.append((path, exp, act))


def differences(expected, actual):
    """{top-level field: description} of where actual differs from expected, empty if none (within tolerance)."""
    found = []
    _walk(expected, actual, (), found)
    by_field = {}
    for path, a, b in found:
        by_field.setdefault(path[0], []).append((path, a, b))
    out = {}
    for field, items in by_field.items():
        if field == "lines":
            out[field] = f"{len(items)} values differ, e.g. {_fmt(*items[0])}"
        else:
            out[field] = "; ".join(_fmt(*i) for i in items[:3]) + (f"; ... {len(items)} in all" if len(items) > 3 else "")
    return out


def _fmt(path, a, b):
    name = ".".join(str(p) for p in path)
    return f"{name}: {_clip(a)} -> {_clip(b)}"


def _clip(v, n=70):
    s = repr(v)
    return s if len(s) <= n else s[:n] + "..."
