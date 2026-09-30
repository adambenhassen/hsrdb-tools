"""Count CIF features that need special reading: anisotropic B_ij / beta_ij, dummy atoms, Cartesian coordinates.

python dev/cif_features.py LIST [-j N] > features.tsv   (one row per CIF that has any of them)
"""
import argparse
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

TAGS = {"aniso_B": re.compile(rb"^_atom_site_aniso_B_11\b", re.M),
        "aniso_beta": re.compile(rb"^_atom_site_aniso_beta_11\b", re.M),
        "calc_flag": re.compile(rb"^_atom_site_calc_flag\b", re.M),
        "cartn": re.compile(rb"^_atom_site_Cartn_x\b", re.M),
        "fract": re.compile(rb"^_atom_site_fract_x\b", re.M)}
DUM = re.compile(rb"\sdum\s", re.I)


def scan(path):
    text = Path(path).read_bytes()
    found = {k: bool(p.search(text)) for k, p in TAGS.items()}
    found["dum"] = found.pop("calc_flag") and bool(DUM.search(text))
    found["cartn_only"] = found.pop("cartn") and not found.pop("fract")
    found.pop("fract", None)
    return path, found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("list")
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()
    paths = Path(args.list).read_text().split()
    print("path\taniso_B\taniso_beta\tdum\tcartn_only")
    with ProcessPoolExecutor(args.j) as ex:
        for path, f in ex.map(scan, paths, chunksize=256):
            if any(f.values()):
                print(f"{path}\t{int(f['aniso_B'])}\t{int(f['aniso_beta'])}\t{int(f['dum'])}\t{int(f['cartn_only'])}")


if __name__ == "__main__":
    main()
