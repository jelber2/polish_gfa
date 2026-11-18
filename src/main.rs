mod index;
mod io;
mod model;
mod polish;
mod recruit;
mod tr;
mod util;

use anyhow::{bail, Context, Result};
use clap::Parser;

use model::{Decompressor, Opts};

#[derive(Parser, Debug)]
#[command(
    name = "polish_gfa",
    version,
    about = "Standalone poa-consensus VNTR/tandem-repeat polishing for hifiasm contig GFAs",
    long_about = None
)]
struct Args {
    /// Read files (FASTQ/FASTA, optionally gzip-compressed).
    #[arg(value_name = "READS", num_args = 0..)]
    paths: Vec<String>,

    /// Hifiasm contig GFA to polish (optionally gzip/bgzip-compressed).
    /// Repeatable (or comma-separated).
    #[arg(long = "gfa", value_name = "GFA", required = true, action = clap::ArgAction::Append, value_delimiter = ',')]
    gfa: Vec<String>,

    /// TSV report path [default: <first-gfa-without-.gfa>.poa.tsv].
    #[arg(long = "tsv", value_name = "FILE")]
    tsv: Option<String>,

    /// Stream/POA worker threads [default: all available CPUs].
    #[arg(short = 't', long = "threads", value_name = "INT")]
    threads: Option<usize>,

    /// Flank length used to anchor a repeat.
    #[arg(long = "flank", short = 'f', default_value_t = 120, value_name = "INT")]
    flank: i32,

    /// Minimum spanning reads to rewrite a locus.
    #[arg(long = "min-reads", default_value_t = 6, value_name = "INT")]
    min_reads: i32,

    /// Longest repeat unit considered.
    #[arg(long = "max-period", default_value_t = 200, value_name = "INT")]
    max_period: i32,

    /// Check every Nth read k-mer when recruiting reads (1-64).
    #[arg(long = "seed-step", default_value_t = 2, value_name = "INT")]
    seed_step: i32,

    /// Flank edge length used for seeding (20-120).
    #[arg(long = "edge-len", default_value_t = 60, value_name = "INT")]
    edge_len: i32,

    /// Cap on recruited read bases held in memory; 0 = unlimited.
    #[arg(
        long = "max-candidate-bases",
        default_value_t = 4_000_000_000i64,
        value_name = "INT"
    )]
    max_candidate_bases: i64,

    /// Skip arrays longer than this many bp; 0 = no limit.
    #[arg(
        long = "max-locus-len",
        default_value_t = 100_000i64,
        value_name = "INT"
    )]
    max_locus_len: i64,

    /// Minimum array length considered; use 150-200 for WGS.
    #[arg(long = "min-array", default_value_t = 60, value_name = "INT")]
    min_array: i32,

    /// Build the seed index for this many loci at a time; 0 = all at once.
    #[arg(
        long = "index-chunk",
        default_value_t = 1_000_000usize,
        value_name = "INT"
    )]
    index_chunk: usize,

    /// Read stream backend: auto/rapidgzip use built-in rapidgzip-core.
    #[arg(
        long = "decompressor",
        default_value = "auto",
        value_name = "auto|zlib|rapidgzip"
    )]
    decompressor: String,

    /// rapidgzip-core decoder threads per input file; 0 = auto.
    #[arg(long = "gzip-threads", default_value_t = 0usize, value_name = "INT")]
    gzip_threads: usize,
}

fn main() {
    if let Err(e) = real_main() {
        eprintln!("polish_gfa: {e:#}");
        std::process::exit(1);
    }
}

fn real_main() -> Result<()> {
    let args = Args::parse();

    let reads = args.paths.clone();
    if reads.is_empty() {
        bail!("missing read file(s); usage: polish_gfa --gfa <GFA> <reads...>");
    }
    if args.gfa.is_empty() {
        bail!("missing --gfa <assembly.gfa>");
    }

    let threads = args
        .threads
        .unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        })
        .max(1);

    if args.flank < 40 {
        bail!("--flank must be >= 40 (got {})", args.flank);
    }
    if args.min_reads < 1 {
        bail!("--min-reads must be >= 1 (got {})", args.min_reads);
    }
    if args.max_period < 2 {
        bail!("--max-period must be >= 2 (got {})", args.max_period);
    }
    if !(1..=64).contains(&args.seed_step) {
        bail!(
            "--seed-step must be between 1 and 64 (got {})",
            args.seed_step
        );
    }
    if !(20..=120).contains(&args.edge_len) {
        bail!(
            "--edge-len must be between 20 and 120 (got {})",
            args.edge_len
        );
    }
    if args.max_candidate_bases < 0 {
        bail!(
            "--max-candidate-bases must be >= 0 (got {})",
            args.max_candidate_bases
        );
    }
    if args.max_locus_len < 0 {
        bail!("--max-locus-len must be >= 0 (got {})", args.max_locus_len);
    }
    if args.min_array < 40 {
        bail!("--min-array must be >= 40 (got {})", args.min_array);
    }

    let decompressor = match args.decompressor.as_str() {
        "auto" => Decompressor::Auto,
        "zlib" => Decompressor::Zlib,
        "rapidgzip" => Decompressor::Rapidgzip,
        other => bail!("--decompressor must be auto, zlib or rapidgzip (got '{other}')"),
    };

    let gzip_threads = if args.gzip_threads > 0 {
        args.gzip_threads
    } else {
        (threads / 2).max(1) / reads.len().clamp(1, 4)
    };

    let opts = Opts {
        gfas: args.gfa.clone(),
        tsv: args.tsv.clone(),
        reads,
        threads,
        flank_len: args.flank,
        min_reads: args.min_reads,
        max_period: args.max_period,
        seed_step: args.seed_step,
        edge_len: args.edge_len,
        max_cand_bases: args.max_candidate_bases,
        max_locus_len: args.max_locus_len,
        min_array: args.min_array,
        index_chunk: args.index_chunk,
        decompressor,
        gzip_threads,
        debug_dir: std::env::var("POA_POLISH_DEBUG_DIR").ok(),
    };

    let code =
        polish::run(&opts).with_context(|| format!("polishing {} GFA(s)", opts.gfas.len()))?;
    std::process::exit(code);
}
