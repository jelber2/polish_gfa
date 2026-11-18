use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::thread;

use crossbeam_channel::{bounded, unbounded, Sender};

use crate::index::{SeedIndex, NONE};
use crate::io::{open_source, FastxReader};
use crate::model::{CutRecord, FoundHit, Locus, Opts, RecruitStats};
use crate::util::{
    base_code, find_subslice, median, revcomp, rfind_subslice, MAX_READS_PER_LOCUS, RECRUIT_K,
};

struct Batch {
    reads: Vec<Vec<u8>>,
    keys: Vec<u64>,
    bases: usize,
}

impl Batch {
    fn new() -> Self {
        Self {
            reads: Vec::new(),
            keys: Vec::new(),
            bases: 0,
        }
    }
    fn take(&mut self) -> Batch {
        std::mem::take(self)
    }
}

impl Default for Batch {
    fn default() -> Self {
        Self::new()
    }
}

struct RecruitCtx {
    found: Vec<FoundHit>,
    scratch: Vec<u8>,
}

fn send_batch(tx: &Sender<Batch>, batch: &mut Batch) {
    if batch.reads.is_empty() {
        return;
    }
    let b = batch.take();
    let _ = tx.send(b);
}

fn producer(
    opts: &Opts,
    next_file: &AtomicUsize,
    reads: &AtomicU64,
    bases: &AtomicU64,
    soft_capped: &AtomicBool,
    tx: &Sender<Batch>,
) {
    loop {
        if soft_capped.load(Ordering::Relaxed) {
            break;
        }
        let f = next_file.fetch_add(1, Ordering::Relaxed);
        if f >= opts.reads.len() {
            break;
        }
        let src = match open_source(&opts.reads[f], opts.decompressor, opts.gzip_threads) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("[polish_gfa] warning: cannot open {}: {e}", opts.reads[f]);
                continue;
            }
        };
        let mut fx = FastxReader::new(src);
        let mut batch = Batch::new();
        let mut ridx: u64 = 0;
        loop {
            if soft_capped.load(Ordering::Relaxed) {
                break;
            }
            let seq = match fx.next_seq() {
                Ok(Some(s)) => s,
                Ok(None) => break,
                Err(e) => {
                    eprintln!("[polish_gfa] warning: read error in {}: {e}", opts.reads[f]);
                    break;
                }
            };
            let len = seq.len() as u64;
            batch.keys.push(((f as u64) << 40) | ridx);
            ridx += 1;
            batch.bases += seq.len();
            batch.reads.push(seq);

            let nr = reads.fetch_add(1, Ordering::Relaxed) + 1;
            bases.fetch_add(len, Ordering::Relaxed);
            if nr % 1_000_000 == 0 {
                eprintln!(
                    "[polish_gfa] streamed {} reads ({:.2} Gbp)",
                    nr,
                    bases.load(Ordering::Relaxed) as f64 / 1e9
                );
            }
            if batch.reads.len() >= 256 || batch.bases >= 2 * 1024 * 1024 {
                send_batch(tx, &mut batch);
            }
        }
        send_batch(tx, &mut batch);
    }
}

#[allow(clippy::too_many_arguments)]
fn recruit_read(
    read: &[u8],
    key: u64,
    loci: &[Locus],
    idx: &SeedIndex,
    ctx: &mut RecruitCtx,
    out: &mut Vec<CutRecord>,
    edge_len: i32,
    seed_step: i32,
    cap_remaining: Option<u64>,
    soft_bases: &AtomicU64,
    soft_capped: &AtomicBool,
) {
    let n = read.len() as i64;
    ctx.found.clear();
    let mask = (1u64 << (2 * RECRUIT_K)) - 1;
    let mut code = 0u64;
    let mut run = 0usize;

    for p in 0..read.len() {
        let b = base_code(read[p]);
        if b < 0 {
            code = 0;
            run = 0;
            continue;
        }
        code = ((code << 2) | b as u64) & mask;
        run += 1;
        if run < RECRUIT_K {
            continue;
        }
        let start = (p + 1 - RECRUIT_K) as i64;
        if seed_step > 1 && start % seed_step as i64 != 0 {
            continue;
        }
        let Some(&h0) = idx.head.get(&code) else {
            continue;
        };
        let mut h = h0;
        let mut chain = 0;
        while h != NONE {
            chain += 1;
            if chain > 512 {
                break;
            }
            let q = idx.hits[h as usize].packed;
            ctx.found.push(FoundHit {
                locus: (q >> 15) as u32,
                strand: ((q >> 14) & 1) as u8,
                side: ((q >> 12) & 3) as u8,
                offset: (q & 0xFFF) as i32,
                pos: start,
            });
            if ctx.found.len() > 500_000 {
                ctx.found.clear();
                return;
            }
            h = idx.hits[h as usize].next;
        }
    }

    if ctx.found.is_empty() {
        return;
    }
    ctx.found
        .sort_by(|a, b| (a.locus, a.strand, a.pos).cmp(&(b.locus, b.strand, b.pos)));

    let mut have_rc = false;
    let mut i = 0usize;
    while i < ctx.found.len() {
        let mut j = i + 1;
        while j < ctx.found.len()
            && ctx.found[j].locus == ctx.found[i].locus
            && ctx.found[j].strand == ctx.found[i].strand
        {
            j += 1;
        }
        let strand = ctx.found[i].strand;
        let lid = ctx.found[i].locus as usize;
        let loc = &loci[lid];
        let region_len = loc.end - loc.start;
        let e = edge_len.min(loc.flank);

        let mut ls_v: Vec<i64> = Vec::new();
        let mut le_v: Vec<i64> = Vec::new();
        let mut rs_v: Vec<i64> = Vec::new();
        let mut re_v: Vec<i64> = Vec::new();
        let mut l_start_min_off = i32::MAX;
        let mut r_end_max_off: i32 = -1;
        let mut n_l_end = 0;
        let mut n_r_start = 0;

        for f in &ctx.found[i..j] {
            let mut pos = f.pos;
            let mut off = f.offset;
            if strand == 1 {
                pos = n - pos - RECRUIT_K as i64;
                off = e - off - RECRUIT_K as i32;
            }
            match f.side {
                0 => {
                    ls_v.push(pos - off as i64);
                    l_start_min_off = l_start_min_off.min(off);
                }
                1 => {
                    le_v.push(pos + (e - off) as i64);
                    n_l_end += 1;
                }
                2 => {
                    rs_v.push(pos - off as i64);
                    n_r_start += 1;
                }
                _ => {
                    re_v.push(pos + (e - off) as i64);
                    r_end_max_off = r_end_max_off.max(off);
                }
            }
        }
        i = j;

        if ls_v.is_empty() || re_v.is_empty() || le_v.is_empty() || rs_v.is_empty() {
            continue;
        }
        if n < (0.7 * region_len as f64) as i64 {
            continue;
        }
        if l_start_min_off > 20 || r_end_max_off < e - RECRUIT_K as i32 - 20 {
            continue;
        }
        if n_l_end < 1 || n_r_start < 1 {
            continue;
        }

        let ls = median(&mut ls_v);
        let le = median(&mut le_v);
        let rs = median(&mut rs_v);
        let re = median(&mut re_v);
        if !(ls < le && le <= rs && rs < re) {
            continue;
        }
        let mut span = rs - le;

        if strand == 1 && !have_rc {
            ctx.scratch = revcomp(read);
            have_rc = true;
        }
        let rd: &[u8] = if strand == 0 { read } else { &ctx.scratch };

        let mut lprobe: i64 = -1;
        let mut loff: i32 = -1;
        let mut rprobe: i64 = -1;
        let mut roff: i32 = -1;
        for off in [
            loc.flank - 16,
            loc.flank - 24,
            loc.flank - 32,
            loc.flank - 40,
        ] {
            if off < 0 || (off + 16) as usize > loc.left.len() {
                continue;
            }
            if let Some(p) = find_subslice(rd, &loc.left[off as usize..(off + 16) as usize]) {
                lprobe = p as i64;
                loff = off;
                break;
            }
        }
        for off in [0, 8, 16, 24] {
            if (off + 16) as usize > loc.right.len() {
                continue;
            }
            if let Some(p) = rfind_subslice(rd, &loc.right[off as usize..(off + 16) as usize]) {
                rprobe = p as i64;
                roff = off;
                break;
            }
        }

        let mut a = ls - 10;
        let mut b = re + 10;
        if lprobe >= 0 && rprobe > lprobe {
            let pspan = (rprobe - roff as i64) - (lprobe - loff as i64 + loc.flank as i64);
            if pspan >= (0.2 * region_len as f64) as i64 {
                span = pspan;
                a = lprobe - loff as i64 - 10;
                b = rprobe + loc.flank as i64 - roff as i64 + 10;
            }
        }
        if span > (5.0 * region_len as f64) as i64 + 2000 {
            continue;
        }
        if a < 0 {
            a = 0;
        }
        if b > n {
            b = n;
        }
        if b - a < span + 2 * loc.flank as i64 - 40 {
            continue;
        }
        if b - a < 100 {
            continue;
        }

        let len = (b - a) as u64;
        if let Some(cap) = cap_remaining {
            let cur = soft_bases.load(Ordering::Relaxed);
            if cur + len > cap {
                soft_capped.store(true, Ordering::Relaxed);
                return;
            }
            soft_bases.fetch_add(len, Ordering::Relaxed);
        }
        out.push(CutRecord {
            key,
            locus: lid,
            strand,
            seq: rd[a as usize..b as usize].to_vec(),
            span,
        });
    }
}

fn worker(
    opts: &Opts,
    loci: &[Locus],
    idx: &SeedIndex,
    cap_remaining: Option<u64>,
    soft_bases: &AtomicU64,
    soft_capped: &AtomicBool,
    rx: crossbeam_channel::Receiver<Batch>,
    tx: Sender<Vec<CutRecord>>,
) {
    let mut ctx = RecruitCtx {
        found: Vec::new(),
        scratch: Vec::new(),
    };
    while let Ok(batch) = rx.recv() {
        let mut local = Vec::new();
        for (read, &key) in batch.reads.iter().zip(batch.keys.iter()) {
            recruit_read(
                read,
                key,
                loci,
                idx,
                &mut ctx,
                &mut local,
                opts.edge_len,
                opts.seed_step,
                cap_remaining,
                soft_bases,
                soft_capped,
            );
        }
        if !local.is_empty() {
            let _ = tx.send(local);
        }
    }
}

fn merge_records(
    records: &mut Vec<CutRecord>,
    loci: &mut [Locus],
    cap_remaining: Option<u64>,
    stats: &mut RecruitStats,
) {
    records.sort_by(|a, b| (a.key, a.locus, a.strand).cmp(&(b.key, b.locus, b.strand)));
    let mut counts = vec![0u32; loci.len()];
    let mut total: u64 = 0;
    for r in records.drain(..) {
        if counts[r.locus] as usize >= MAX_READS_PER_LOCUS {
            continue;
        }
        if let Some(cap) = cap_remaining {
            if total + r.seq.len() as u64 > cap {
                stats.capped = true;
                break;
            }
        }
        total += r.seq.len() as u64;
        counts[r.locus] += 1;
        stats.stored_bases += r.seq.len() as u64;
        stats.stored_reads += 1;
        loci[r.locus].spans.push(r.span);
        loci[r.locus].cut_reads.push(r.seq);
    }
}

pub fn stream_recruit(
    opts: &Opts,
    loci: &mut Vec<Locus>,
    idx: &SeedIndex,
    cap_total: Option<u64>,
    already_stored: u64,
) -> RecruitStats {
    let mut stats = RecruitStats::default();
    let cap_remaining = match cap_total {
        Some(c) if already_stored >= c => {
            stats.capped = true;
            return stats;
        }
        Some(c) => Some(c - already_stored),
        None => None,
    };

    let nthreads = opts.threads.max(1);
    let nprod = opts.reads.len().clamp(1, 4);
    let soft_bases = AtomicU64::new(0);
    let soft_capped = AtomicBool::new(false);
    let reads = AtomicU64::new(0);
    let bases = AtomicU64::new(0);
    let next_file = AtomicUsize::new(0);

    let mut records = {
        let loci_ref: &[Locus] = loci;
        thread::scope(|scope| {
            let (batch_tx, batch_rx) = bounded::<Batch>((4 * nthreads).max(8));
            let (res_tx, res_rx) = unbounded::<Vec<CutRecord>>();

            let next_file_ref = &next_file;
            let reads_ref = &reads;
            let bases_ref = &bases;
            let capped_ref = &soft_capped;
            let mut prod_handles = Vec::new();
            for _ in 0..nprod {
                let tx = batch_tx.clone();
                prod_handles.push(scope.spawn(move || {
                    producer(opts, next_file_ref, reads_ref, bases_ref, capped_ref, &tx);
                }));
            }
            let soft_bases_ref = &soft_bases;
            let capped_ref2 = &soft_capped;
            let mut work_handles = Vec::new();
            for _ in 0..nthreads {
                let rx = batch_rx.clone();
                let tx = res_tx.clone();
                work_handles.push(scope.spawn(move || {
                    worker(
                        opts,
                        loci_ref,
                        idx,
                        cap_remaining,
                        soft_bases_ref,
                        capped_ref2,
                        rx,
                        tx,
                    );
                }));
            }
            drop(batch_tx);
            drop(batch_rx);
            drop(res_tx);
            for h in prod_handles {
                let _ = h.join();
            }
            for h in work_handles {
                let _ = h.join();
            }

            let mut all = Vec::new();
            for v in res_rx.iter() {
                all.extend(v);
            }
            all
        })
    };

    stats.reads = reads.load(Ordering::Relaxed);
    stats.bases = bases.load(Ordering::Relaxed);
    merge_records(&mut records, loci, cap_remaining, &mut stats);
    stats.capped |= soft_capped.load(Ordering::Relaxed);

    if stats.capped {
        let capv = cap_total.unwrap_or(0) as f64;
        if capv >= 1e9 {
            eprintln!(
                "[polish_gfa] candidate read cap reached ({:.2} Gbp); raise --max-candidate-bases to scan the remaining reads",
                capv / 1e9
            );
        } else {
            eprintln!(
                "[polish_gfa] candidate read cap reached ({:.2} Mbp); raise --max-candidate-bases to scan the remaining reads",
                capv / 1e6
            );
        }
    }
    stats
}
