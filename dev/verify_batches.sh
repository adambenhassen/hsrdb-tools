#!/bin/sh
# Run dev/verify_cctbx.py over a CIF list in batches, resumably: verify/batches/NNN.tsv per finished batch (NNN = aa, ab, ...).
# usage (repo root, host with docker): sh dev/verify_batches.sh LIST BATCH_SIZE WORKERS [--db DB]
set -eu
list=$1 size=$2 jobs=$3
shift 3
mkdir -p verify/batches
split -l "$size" -a 2 "$list" verify/batches/list.
for part in verify/batches/list.*; do
  n=${part##*.}
  [ -f "verify/batches/$n.tsv" ] && continue
  docker run --rm -v "$PWD:/w" hsrdb-dev /opt/cctbx/bin/python dev/verify_cctbx.py "$part" -j "$jobs" "$@" \
    > "verify/batches/$n.partial" 2> "verify/batches/$n.err"
  mv "verify/batches/$n.partial" "verify/batches/$n.tsv"
  echo "batch $n done $(date -u +%H:%M)"
done
echo "all batches done"
