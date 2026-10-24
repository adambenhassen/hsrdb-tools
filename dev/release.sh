#!/bin/sh
# Release pipeline: build, audit, compare, verify, pack. No publishing.
#
# Run from the repo root on a host with docker and the hsrdb-dev image (docker build --target dev -t hsrdb-dev .).
# The repo is mounted at /w in the container, so MIRROR, OUT, COD24 and VERIFY_LIST must be relative paths inside it.
#
# Settings (environment, defaults in brackets):
#   DATE         date stamp in file names [today, YYYY_MM_DD]
#   MIRROR       directory of COD CIFs [cod-mirror/cif]
#   OUT          output directory [out]
#   JOBS         worker processes [4]
#   IMAGE        docker image [hsrdb-dev]
#   COD24        official COD24 HS4 database; the comparison is skipped if the file is missing [COD24_HS4x.hsrdb]
#   VERIFY_LIST  file listing CIF paths to cross-check with cctbx (all of COD, or only new and changed CIFs);
#                the cctbx check is skipped if unset []
#   RESUME       non-empty: continue an interrupted build from its last checkpoint (hsrdb-tools add --resume) []
#
# A failed build keeps $OUT/<name>.hsrdb.partial (up to a full database each, about 7 GB for all of COD) when it got
# past its first checkpoint; a new build refuses to start over them. Rerun with RESUME=1, or delete them.
#
# Each step appends "name exit-status" to $OUT/post.log; the first failure appends FAILED and stops the run.
set -eu

DATE=${DATE:-$(date +%Y_%m_%d)}
MIRROR=${MIRROR:-cod-mirror/cif}
OUT=${OUT:-out}
JOBS=${JOBS:-4}
IMAGE=${IMAGE:-hsrdb-dev}
COD24=${COD24:-COD24_HS4x.hsrdb}
VERIFY_LIST=${VERIFY_LIST:-}
RESUME=${RESUME:-}

D=COD_$DATE
LOG=$OUT/post.log

mkdir -p "$OUT"
: > "$LOG"

dk() { docker run --rm -e OMP_NUM_THREADS=1 -v "$PWD:/w" "$IMAGE" "$@"; }

# step NAME COMMAND...: run, log the exit status, stop on failure. Redirections on the call apply to the command.
step() {
    name=$1
    shift
    rc=0
    "$@" || rc=$?
    echo "$name $rc" >> "$LOG"
    [ "$rc" = 0 ] || { echo FAILED >> "$LOG"; exit 1; }
}

[ -n "$RESUME" ] || : > "$OUT/build.log"  # a resumed build appends: the first part's tracebacks stay
# shellcheck disable=SC2086  # RESUME expands to --resume or to nothing
step build dk python -m hsrdb_tools add "$MIRROR" --hs4 "$OUT/${D}_HS4x.hsrdb" --hs3 "$OUT/${D}_HS3x.hsrdb" -j "$JOBS" \
    ${RESUME:+--resume} >> "$OUT/build.log" 2>&1

step audit_hs4 dk python dev/audit_database.py "$OUT/${D}_HS4x.hsrdb" "$JOBS" > "$OUT/audit_hs4.txt" 2>&1
step audit_hs3 dk python dev/audit_database.py "$OUT/${D}_HS3x.hsrdb" "$JOBS" > "$OUT/audit_hs3.txt" 2>&1

if [ -f "$COD24" ]; then
    step compare_cod24 dk python dev/compare_databases.py "$OUT/${D}_HS4x.hsrdb" "$COD24" all "$JOBS" \
        > "$OUT/compare_cod24.txt" 2>&1
else
    echo "compare_cod24 skipped ($COD24 not found)" | tee -a "$LOG"
fi

if [ -n "$VERIFY_LIST" ]; then
    step verify dk /opt/cctbx/bin/python dev/verify_cctbx.py "$VERIFY_LIST" -j "$JOBS" --db "$OUT/${D}_HS4x.hsrdb" \
        > "$OUT/verify.tsv" 2> "$OUT/verify.err"
else
    echo "verify skipped (VERIFY_LIST not set)" | tee -a "$LOG"
fi

step skipped_md dk python dev/skipped_md.py "$OUT/${D}_HS4x.skipped.tsv" > "$OUT/SKIPPED.md" 2> "$OUT/skipped_md.err"

# 7z is only in the image; it packs from the output directory.
pack() { docker run --rm -v "$PWD/$OUT:/o" -w /o "$IMAGE" 7z "$@"; }

for v in HS4x HS3x; do
    rm -f "$OUT/${D}_$v.7z"  # 7z a adds to an existing archive: a rerun must not keep old contents
    step "pack $v" pack a -mx=7 "${D}_$v.7z" "${D}_$v.hsrdb" > "$OUT/pack_$v.log" 2>&1
    step "test $v" pack t -slt "${D}_$v.7z" >> "$OUT/pack_$v.log" 2>&1
done

sums() (
    cd "$OUT"
    if command -v sha256sum > /dev/null 2>&1; then
        sha256sum "${D}_HS4x.7z" "${D}_HS3x.7z"
    else
        shasum -a 256 "${D}_HS4x.7z" "${D}_HS3x.7z"
    fi
)
step sha sums > "$OUT/SHA256SUMS_$DATE.txt"

echo DONE >> "$LOG"
