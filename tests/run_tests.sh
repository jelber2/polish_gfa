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

run base reads.fq.gz -t 8
check base.poa.tsv base.poa.tsv
check base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.gfa
check base.bp.p_ctg.poa.fa base.bp.p_ctg.poa.fa

rm -f base.poa.tsv base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.fa
run degraded reads.fq.gz -t 8 --index-chunk 2
check degraded.poa.tsv degraded.poa.tsv
check degraded.p_ctg.poa.gfa degraded.p_ctg.poa.gfa
check degraded.p_ctg.poa.fa degraded.p_ctg.poa.fa

if command -v rapidgzip >/dev/null 2>&1; then
    rm -f base.poa.tsv base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.fa
    run base reads.fq.gz -t 8 --decompressor rapidgzip --gzip-threads 4
    check base.poa.tsv base.poa.tsv
    check base.bp.p_ctg.poa.gfa base.bp.p_ctg.poa.gfa
fi

echo "all tests passed"
