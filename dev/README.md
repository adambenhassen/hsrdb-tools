# dev scripts

Scripts used to decode the HighScore format, check built databases and investigate differences from COD24. They are not part of the installed package.

Run from the repository root with the package and its dependencies installed (`pip install -e .`) so that `python dev/NAME.py` finds `hsrdb_tools`. Usage is in each script's docstring or its argument handling. Scripts that read `scratch/cif_cache` or the COD mirror expect those folders in the working directory.

Inputs used below:

- **COD24**: the official `COD24_HS4x.hsrdb` in the working directory.
- **cif_cache**: `scratch/cif_cache/*.cif`, CIFs of entries in the official database.
- **mirror**: the COD CIF mirror in `cod-mirror/cif`.
- **LIST**: a text file with one CIF path per line.
- **cctbx**: the cctbx environment of the docker image `hsrdb-dev` (`docker build --target dev -t hsrdb-dev .`), run as `docker run --rm -v $PWD:/w hsrdb-dev /opt/cctbx/bin/python dev/NAME.py ...`.

## Release

`release.sh` runs the release pipeline (build, audit, compare, verify, pack; no publishing) from the repository root on a docker host. Settings are environment variables, listed at the top of the script. Steps and exit statuses go to `$OUT/post.log`; the first failure stops the run.

```
DATE=2026_10_15 VERIFY_LIST=out/verify.lst sh dev/release.sh
```

A failed build keeps its `.hsrdb.partial` files (up to a full database each, about 7 GB for all of COD) once it has passed its first checkpoint, and a new build refuses to start over them: rerun with `RESUME=1` to continue, or delete them.

## Scripts in dev/

| Script | What it does | Inputs | Release |
|---|---|---|---|
| `release.sh` | Release pipeline, see above | docker, `hsrdb-dev` image, mirror; COD24 and `VERIFY_LIST` optional | pipeline |
| `audit_database.py` | Checks every entry of a built `.hsrdb`: pattern blobs, metadata, structure records | the built database | step |
| `compare_databases.py` | Coverage and pattern/metadata agreement of a built database against the official one; writes outliers to `OURS.compare.tsv` | built database, COD24 | step |
| `verify_cctbx.py` | Compares calculated patterns with cctbx, one TSV row per CIF | LIST or CIFs, cctbx, optional `--db` | step |
| `skipped_md.py` | Markdown table of the CIFs a build skipped; `--recheck` re-runs them first | `*.skipped.tsv` from a build; `--recheck` needs the CIFs | step |
| `verify_explain.py` | Assigns a cause to every discrepancy in `verify_cctbx.py` output | `verify_cctbx.py` TSV, optional `metric_scan.py` output | |
| `verify_modulated.py` | Checks modulated structures against the structure factors listed in their CIFs | CIFs | |
| `metric_scan.py` | CIFs whose cell does not fit the space group | LIST | |
| `cif_features.py` | CIFs with B_ij/beta_ij tensors, dummy atoms or Cartesian coordinates | LIST | |
| `element_formula_scan.py` | CIFs with atoms read as elements their formula does not contain | LIST | |
| `label_element_scan.py` | Atoms whose element is read from an unusual label (no type symbol column) | LIST | |
| `near_special.py` | Sites whose symmetry images are close but do not coincide | LIST | |
| `near_special_summary.py` | Summarises what CIFs state about those sites | `near_special.py` output | |
| `ops_fallback_scan.py` | CIFs whose symmetry sources disagree or that gemmi cannot match to a table setting | LIST | |
| `sym_conflict_check.py` | For those CIFs, which symmetry reading agrees with the CIF's composition and density | `ops_fallback_scan.py` output | |
| `unknown_sg.py` | Symbols and operations of CIFs skipped as "unknown space group" | `*.skipped.tsv` | |
| `recover_check.py` | Checks CIFs recovered from their operations (`setting.from_operations`) | `unknown_sg.py` output | |
| `setting_check.py` | Reference-setting transform must leave pattern and Ueq unchanged (runs at import, no arguments) | cif_cache | |
| `stated_general_check.py` | Atoms on special positions that the CIF states as general: which reading fits | LIST | |
| `default_b_check.py` | Which displacement the official patterns assume for atoms without one (U = 0 or B = 0.5) | LIST, COD24 | |
| `invalid_adp_check.py` | Whether COD24 keeps physically invalid anisotropic tensors | COD ids, COD24, mirror, cctbx | |
| `coincident_proof.py` | Reflections behind one line where we and COD24 disagree (COD 9009599 by default) | COD24, mirror | |
| `outlier_by_laue.py` | Rate of strong-line outliers per Laue class | `compare_databases.py` outliers file | |
| `equivalence.py` | Official phase record vs ours: same structure, compared by powder pattern | COD24, cif_cache | |
| `extract_templates.py` | Derives the constant parts of phase records and writes `hsrdb_tools/templates.json` | COD24 | |
| `update_expected.py` | Regenerates `tests/data/cod/expected.json`, the stored output of the regression test on real COD CIFs, and prints which CIFs and fields changed | `tests/data/cod/*.cif` | |

## Scripts in dev/archive/

One-off investigations, kept for reference. Paths in them are relative to the repository root.

| Script | What it does | Inputs |
|---|---|---|
| `compare_fields.py` | Rebuilds official entries from their CIFs and diffs every metadata field | COD24, cif_cache |
| `compare_patterns.py` | Compares calculated patterns with stored ones for a sample of entries | COD24, cif_cache, network (fetches CIFs missing from the cache) |
| `compare_phase.py` | Builds phase records for official entries and diffs them field by field | COD24, cif_cache |
| `roundtrip_phase.py` | Decodes and re-encodes official phase records; reports mismatches | COD24 |
| `phase_dump.py` | Token dump of a decompressed phase record (strings, doubles, bytes) | a decompressed record file |
| `phase_skeleton.py` | Skeleton of a decompressed phase record: values by name, gaps as hex | COD24 |
| `v3_scan.py` | Lists ANSI strings in a version-3 phase record | a decompressed record file |
| `skipped_in_official.py` | Which skipped CIFs the official database contains | `*.skipped.tsv`, COD24 |
| `skipped_detail.py` | For those: what the official entry holds and what the CIF says now | `*.skipped.tsv`, COD24 |
| `outlier_causes.py` | Samples large intensity outliers and classifies their causes | built database, COD24, mirror |
| `outlier_detail.py` | Line-by-line pattern comparison for outlier entries with unchanged CIFs | built database, COD24, mirror |
| `scan_adp.py` | Counts CIFs with invalid anisotropic displacement tensors | a folder of CIFs |
| `near_special_formula.py` | Which reading of near-special sites matches the CIF formula | `near_special.py` output, cctbx |
| `verify_detail.py` | Site-by-site view of one CIF for a `verify_cctbx.py` discrepancy | a CIF, cctbx |
| `verify_summary.py` | Summary of `verify_cctbx.py` output | `verify_cctbx.py` TSV |
| `verify_batches.sh` | Runs `verify_cctbx.py` over a LIST in resumable batches | LIST, docker, `hsrdb-dev` image |
