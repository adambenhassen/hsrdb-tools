"""Per-CIF work for the parallel build, in its own module so worker processes can import it."""

import traceback

from . import entry, hsrdb, structure

# problems in the CIF itself: the entry is skipped and reported; anything else is a bug in this tool
DATA_ERRORS = (entry.EntryError,)


def prepare(job):
    """Worker: parse one CIF and encode everything that goes into the databases."""
    path, targets = job
    try:
        e = entry.from_cif(path)
        sites = structure.site_info(e.structure, e.images)
        blobs = {t: hsrdb.phase_blob(e, t, sites) for t in targets}
        lines = hsrdb.encode_lines(e.lines)
    except DATA_ERRORS as exc:
        return path, None, None, None, ("skipped", f"{type(exc).__name__}: {exc}")
    except Exception:  # a bug: report it with its traceback instead of hiding it as a bad CIF
        return path, None, None, None, ("internal error", traceback.format_exc())
    e.structure = e.block = e.modulated = e.images = e.scattering = e.site_orders = None  # gemmi objects do not cross process boundaries
    return path, e, lines, blobs, None
