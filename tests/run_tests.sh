#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
bin="$root/target/release/polish_gfa"
if [[ ! -x "$bin" ]]; then
    (cd "$root" && cargo build --release)
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cp "$here/data/ont_vntr_30x/reads.fq.gz" \
   "$here/data/ont_vntr_30x/base.bp.p_ctg.gfa" \
   "$here/data/ont_vntr_30x/degraded.p_ctg.gfa" "$tmp/"
cd "$tmp"

run() {
    if ! "$bin" "$@" > run.log 2>&1; then
        cat run.log
        exit 1
    fi
}

check() {
    local got=$1 exp=$2
    if ! cmp -s "$got" "$here/expected/$exp"; then
        echo "FAIL: $got differs from expected/$exp"
        exit 1
    fi
    echo "ok: $got"
}

# The TSV's first column is the GFA path as given, so compare the body only
# when the same assembly is supplied under a different filename.
check_body() {
    local got=$1 exp=$2
    if ! cmp -s <(cut -f2- "$got") <(cut -f2- "$here/expected/$exp"); then
        echo "FAIL: $got body differs from expected/$exp"
        exit 1
    fi
    echo "ok: $got (body)"
}

run base reads.fq.gz --gfa base.bp.p_ctg.gfa -t 8
if ! tail -n 1 run.log | grep -q "peak RSS"; then
    echo "FAIL: last log line is not the peak RSS report"
    exit 1
fi
echo "ok: peak RSS is the last log line"
check base.bp.p_ctg.poa.tsv base.poa.tsv
check base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.gfa
check base.bp.p_ctg.poa.fa base.bp.p_ctg.poa.fa

# gzipped GFA input
gzip -c base.bp.p_ctg.gfa > base.bp.p_ctg.gfa.gz
rm -f base.bp.p_ctg.poa.tsv base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.fa
run base reads.fq.gz --gfa base.bp.p_ctg.gfa.gz -t 8
check_body base.bp.p_ctg.poa.tsv base.poa.tsv
check base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.gfa
check base.bp.p_ctg.poa.fa base.bp.p_ctg.poa.fa

rm -f base.bp.p_ctg.poa.tsv base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.fa
run degraded reads.fq.gz --gfa degraded.p_ctg.gfa -t 8 --index-chunk 2
check degraded.p_ctg.poa.tsv degraded.poa.tsv
check degraded.p_ctg.poa.gfa degraded.p_ctg.poa.gfa
check degraded.p_ctg.poa.fa degraded.p_ctg.poa.fa

rm -f base.bp.p_ctg.poa.tsv base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.fa
run base reads.fq.gz --gfa base.bp.p_ctg.gfa -t 8 --decompressor zlib
check base.bp.p_ctg.poa.tsv base.poa.tsv
check base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.gfa

rm -f base.bp.p_ctg.poa.tsv base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.fa
run base reads.fq.gz --gfa base.bp.p_ctg.gfa -t 8 --decompressor rapidgzip --gzip-threads 4
check base.bp.p_ctg.poa.tsv base.poa.tsv
check base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.gfa

echo "all tests passed"
