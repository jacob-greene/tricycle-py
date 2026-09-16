#!/usr/bin/env bash
# Runtime and peak memory for both implementations, measured the same way.
#
# Each measurement is a whole separate process wrapped in /usr/bin/time -v, so
# "peak memory" means peak resident set size of that process. A "load" run
# measures the same process doing everything except the function under test,
# which gives the baseline the other rows should be read against.
#
# Usage:
#   equivalence/benchmark.sh <input.rds> <input.h5ad> <outdir> [replicate_factor]
set -u -o pipefail

RDS="$1"; H5AD="$2"; OUT="$3"; REPS="${4:-1}"
HERE="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-python}"
RSCRIPT="${RSCRIPT:-Rscript}"

mkdir -p "$OUT"
RESULT="$OUT/benchmark.tsv"
printf 'side\ttarget\treplicates\tcells\tseconds\tpeak_rss_kb\n' > "$RESULT"

measure() {   # side command...
    local side="$1"; shift
    local target="$1"; shift
    local log="$OUT/${side}_${target}.log"
    /usr/bin/time -v -o "$log.time" "$@" > "$log" 2>&1
    local cells seconds rss
    cells=$(awk -F'\t' '$1=="cells"{print $2}' "$log")
    seconds=$(awk -F'\t' '$1=="seconds"{print $2}' "$log")
    rss=$(awk '/Maximum resident set size/{print $NF}' "$log.time")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$side" "$target" "$REPS" "${cells:-NA}" "${seconds:-NA}" "${rss:-NA}" \
        >> "$RESULT"
}

for target in load position schwabe loess_direct loess_interpolate; do
    measure R "$target" "$RSCRIPT" "$HERE/bench_r.R" "$RDS" "$target" "$REPS"
done
# PYTHONPATH is cleared for the Python runs. An environment-module R often
# exports a PYTHONPATH of its own, and it shadows the virtual environment's
# packages with older ones, which fails at import time.
for target in load position schwabe loess_direct; do
    measure Python "$target" env -u PYTHONPATH "$PYTHON" "$HERE/bench_py.py" \
        "$H5AD" "$target" "$REPS"
done

column -t -s $'\t' "$RESULT"
