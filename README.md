# polish_gfa

Standalone Rust implementation of the poa-consensus VNTR/tandem-repeat
polishing pass. It takes an existing hifiasm contig GFA and the reads that
produced it, rebuilds tandem arrays from the spanning reads with
[poa-consensus](https://github.com/Psy-Fer/poa-consensus) single-allele flanked
mode, and writes a polished GFA/FASTA plus a per-locus TSV report. The input
GFA is never modified.

This is a Rust-only port of the `poa-test/polish_gfa` tool that was developed
inside hifiasm-fast. It calls the `poa-consensus` crate directly (no C++ and no
C ABI shim) and is independent of that tree: hifiasm-fast keeps its in-tree
integration and `poa-test/polish_gfa`, while this repository provides the
standalone Rust implementation with byte-identical outputs.

## Build

Requires Rust ≥ 1.88 (both `poa-consensus` and `rapidgzip-core` are edition
2024). The first build fetches the pinned crates from GitHub.

```bash
cargo build --release
# -> target/release/polish_gfa
```

Gzip input is decoded by the built-in
[rapidgzip-core](https://github.com/COMBINE-lab/rapidgzip-rust) crate (parallel
DEFLATE decoding of a single gzip file; no external `rapidgzip` binary is
needed). `--decompressor zlib` selects the slower single-threaded flate2 path.

## Usage

```bash
polish_gfa [options] --gfa <assembly.gfa> <reads.fq[.gz]> [reads2.fq.gz ...]

# multiple GFAs (repeat --gfa, or use a comma-separated list)
polish_gfa --gfa asm.bp.p_ctg.gfa --gfa asm.bp.hap1.p_ctg.gfa reads.fq.gz
polish_gfa --gfa asm.bp.p_ctg.gfa,asm.bp.hap1.p_ctg.gfa reads.fq.gz
```

`--gfa` takes the GFA file directly, so any hifiasm contig graph can be
polished without knowing the assembly prefix. Outputs per GFA:

```
<gfa-without-.gfa>.poa.gfa   patched contig graph, A-line offsets shifted
<gfa-without-.gfa>.poa.fa    patched contigs as FASTA
<first-gfa-without-.gfa>.poa.tsv
                             one row per locus: period, copies, lengths,
                             reads, k-mer support before/after, status
```

Use `--tsv <FILE>` to choose the report path when polishing several GFAs.

### Options

| Option | Default | Meaning |
|---|---|---|
| `--gfa GFA` | required | hifiasm contig GFA to polish; repeatable or comma-separated |
| `--tsv FILE` | `<first-gfa>.poa.tsv` | TSV report path |
| `-t, --threads INT` | all CPUs | stream/POA worker threads |
| `-f, --flank INT` | 120 | flank length used to anchor a repeat (≥40) |
| `--min-reads INT` | 6 | minimum spanning reads to rewrite a locus (effective min 3) |
| `--max-period INT` | 200 | longest repeat unit considered |
| `--seed-step INT` | 2 | check every Nth read k-mer (1–64) |
| `--edge-len INT` | 60 | flank edge length used for seeding (20–120) |
| `--max-candidate-bases INT` | 4,000,000,000 | cap on trimmed reads held in memory; `0` = unlimited |
| `--max-locus-len INT` | 100,000 | skip arrays longer than this; `0` = no limit |
| `--min-array INT` | 60 | minimum array length considered; use 150–200 for WGS |
| `--index-chunk INT` | 1,000,000 | build the seed index for this many loci at a time; `0` = all at once |
| `--decompressor STR` | auto | `auto`/`rapidgzip` use built-in rapidgzip-core; `zlib` uses flate2 |
| `--gzip-threads INT` | 0 (auto) | rapidgzip-core decoder threads per input file |

Examples:

```bash
# small test assembly
polish_gfa --gfa asm.bp.p_ctg.gfa reads.fq.gz

# WGS: fewer short-STR loci, 20 Gbp of candidate reads, parallel gzip
polish_gfa -t32 --min-array 150 --max-candidate-bases 20000000000 \
    --decompressor rapidgzip --gzip-threads 16 \
    --gfa asm.bp.p_ctg.gfa reads.fq.gz

# keep short STRs but bound index memory: 1M loci per chunk
polish_gfa -t32 --index-chunk 1000000 --max-candidate-bases 20000000000 \
    --gfa asm.bp.p_ctg.gfa reads.fq.gz
```

## How it works

1. **Find tandem arrays** in every contig (period 2–`--max-period`, length
   ≥`--min-array`, ≥2.5 copies), with tolerant phase drift and merging of
   overlapping calls.
2. **Index only the flanks.** All 17-mers of four 60 bp flank edges (both
   strands) are indexed per locus; low-complexity seeds are skipped. The index
   is built in chunks of `--index-chunk` loci so memory does not scale with the
   total number of arrays.
3. **Stream the reads once per chunk.** A producer/worker pipeline
   decompresses/parses the input, looks up every 2nd k-mer, and keeps only
   reads that reach the distal end of the left flank, span the repeat, and
   extend past the right flank. Cut-down reads are capped at 128 per locus and
   `--max-candidate-bases` globally.
4. **Call poa-consensus** `consensus_flanked()` on the trimmed reads in
   parallel across loci, with the crate's ONT defaults (homopolymer and tandem
   gap costs).
5. **Splice only when credible**: both flanks must be anchored in the full
   consensus, the length must be within 0.4–2.5×, k-mer support must not drop,
   and the contig length must be an outlier against the median spanning-read
   length and move closer to it. Loci whose length already agrees with the
   reads are skipped before the POA call.

The original GFA is untouched; `A`-line offsets are shifted by the splice
deltas so read-to-contig annotations stay consistent.

## Performance and memory

Read recruitment and POA are parallel across `-t` threads. On the included
30× benchmark (6 arrays, 3 GFAs):

* already-correct assembly: ~0.1 s stream + no POA per GFA;
* controlled 8%-deletion repair (6 loci): ~3 s total, dominated by the 8.1 kb
  array.

The seed index costs roughly 20–25 KB per detected array at a 60 bp edge.
`--min-array 60` can produce millions of short-STR loci in a human genome and
hundreds of GB of index; for WGS use `--min-array 150` (or 200) and/or chunk the
index. Each chunk costs one full read pass per GFA. For example, 8M loci at
1M per chunk is 8 passes, ~22 GB peak index instead of ~180 GB.

`--max-candidate-bases` is a safety limit on the trimmed reads held in memory
(~1.1 bytes/base). It is global across chunks; if reached, the remaining chunks
are skipped. For a 34× 3.2 Gbp run (~108.8 Gbp per GFA) use ~20 Gbp
(≈22 GB), and the default rapidgzip-core backend to avoid the single-threaded
gzip floor.

FASTA input does **not** reduce index memory; it only cuts decompression and
parsing time. `kseq`-style FASTA and FASTQ are both accepted.

## Tests

```bash
./tests/run_tests.sh
```

The included 30× simulated ONT dataset and golden outputs are described in
`tests/README.md`. The Rust outputs are byte-identical to the original C++
driver on both the correct and 8%-deleted assemblies.

## License

MIT OR Apache-2.0.
