# Tests

`run_tests.sh` runs the release binary on the included 30× simulated ONT
dataset and compares the `.poa.tsv`, `.poa.gfa` and `.poa.fa` outputs against
`expected/`. The expected files were produced by the original C++ driver
(`poa-test/polish_gfa` in hifiasm-fast) before it was ported to Rust; the Rust
implementation reproduces them byte-for-byte. The script also runs the damaged
(`degraded.p_ctg.gfa`) assembly with `--index-chunk 2`, and exercises both the
built-in rapidgzip-core backend and the zlib backend.

```bash
./run_tests.sh
```

## Regenerating the data

The dataset was generated with:

```bash
python3 sim_ont_vntr.py \
    --outdir data/ont_vntr_30x \
    --coverage 30 --err-rate 0.045 --hp-del 3 --hp-ins 2 \
    --n50 12000 --seed 5 --panel xl
```

`degrade_contig.py` produced `degraded.p_ctg.gfa` by deleting 8% of every
array. `eval.py` computes the per-locus copy number and identity used in the
main README.

The base GFA (`base.bp.p_ctg.gfa`) comes from running hifiasm on
`reads.fq.gz`:

```bash
hifiasm --ont -t 8 -o base reads.fq.gz
```
