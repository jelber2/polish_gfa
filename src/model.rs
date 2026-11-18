#[derive(Clone, Copy, Debug)]
pub struct TrRegion {
    pub start: i64,
    pub end: i64,
    pub period: i32,
}

#[derive(Clone, Debug)]
pub struct Locus {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub flank: i32,
    pub left: Vec<u8>,
    pub right: Vec<u8>,
    pub cut_reads: Vec<Vec<u8>>,
    pub spans: Vec<i64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Decompressor {
    Auto,
    Zlib,
    Rapidgzip,
}

#[derive(Clone, Debug)]
pub struct Opts {
    pub gfas: Vec<String>,
    pub tsv: Option<String>,
    pub reads: Vec<String>,
    pub threads: usize,
    pub flank_len: i32,
    pub min_reads: i32,
    pub max_period: i32,
    pub seed_step: i32,
    pub edge_len: i32,
    pub max_cand_bases: i64,
    pub max_locus_len: i64,
    pub min_array: i32,
    pub index_chunk: usize,
    pub decompressor: Decompressor,
    pub gzip_threads: usize,
    pub debug_dir: Option<String>,
}

#[derive(Clone, Copy, Default, Debug)]
pub struct RecruitStats {
    pub reads: u64,
    pub bases: u64,
    pub stored_reads: u64,
    pub stored_bases: u64,
    pub capped: bool,
}

#[derive(Clone, Debug)]
pub struct FoundHit {
    pub locus: u32,
    pub pos: i64,
    pub offset: i32,
    pub strand: u8,
    pub side: u8,
}

#[derive(Debug)]
pub struct CutRecord {
    pub key: u64,
    pub locus: usize,
    pub strand: u8,
    pub seq: Vec<u8>,
    pub span: i64,
}

#[derive(Clone, Debug, Default)]
pub struct LocusStat {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub period: i32,
    pub old_len: i64,
    pub new_len: i64,
    pub n_reads: i32,
    pub cont_old: f64,
    pub cont_new: f64,
    pub span_med: i64,
    pub span_min: i64,
    pub span_max: i64,
    pub status: String,
}

#[derive(Clone, Debug, Default)]
pub struct PatchResult {
    pub patch: bool,
    pub seconds: f64,
    pub repeat: Vec<u8>,
    pub stat: LocusStat,
}
