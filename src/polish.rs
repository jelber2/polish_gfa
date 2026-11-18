use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use poa_consensus::{consensus, consensus_flanked, PoaConfig};
use rustc_hash::FxHashSet;

use crate::index::{build_index, SeedIndex};
use crate::model::{Locus, LocusStat, Opts, PatchResult, RecruitStats, TrRegion};
use crate::recruit::stream_recruit;
use crate::tr::TrDetector;
use crate::util::{base_code, find_subslice, strip_gfa_suffix, HIT_KMER, RECRUIT_K};

const POA_BUDGET: usize = 65536;

fn locus_flank(l: i64, tr: &TrRegion, opts: &Opts) -> i32 {
    let mut flank = if opts.flank_len > 0 {
        opts.flank_len
    } else {
        120
    };
    if flank < 40 {
        flank = 40;
    }
    flank.min((tr.start.min(l - tr.end)) as i32)
}

fn parse_s_line(line: &str) -> Option<(String, Vec<u8>, String)> {
    let mut it = line.splitn(4, '\t');
    let tag = it.next()?;
    if tag != "S" {
        return None;
    }
    let name = it.next()?.to_string();
    let seq = it.next()?.as_bytes().to_vec();
    let tags = it.next().unwrap_or("").to_string();
    Some((name, seq, tags))
}

fn collect_loci(gfa: &str, opts: &Opts, detector: &mut TrDetector) -> Result<Vec<Locus>> {
    let f = File::open(gfa).with_context(|| format!("open {gfa}"))?;
    let reader = BufReader::with_capacity(1 << 20, f);
    let mut loci = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let Some((name, seq, _tags)) = parse_s_line(&line) else {
            continue;
        };
        if seq.len() < opts.min_array.max(0) as usize {
            continue;
        }
        let trs = detector.detect(&seq, opts.max_period, opts.min_array);
        for tr in trs {
            if opts.max_locus_len > 0 && tr.end - tr.start > opts.max_locus_len {
                continue;
            }
            let flank = locus_flank(seq.len() as i64, &tr, opts);
            if flank < 40 {
                continue;
            }
            let f = flank as usize;
            loci.push(Locus {
                contig: name.clone(),
                start: tr.start,
                end: tr.end,
                flank,
                left: seq[(tr.start as usize - f)..(tr.start as usize)].to_vec(),
                right: seq[(tr.end as usize)..(tr.end as usize + f)].to_vec(),
                cut_reads: Vec::new(),
                spans: Vec::new(),
            });
        }
    }
    Ok(loci)
}

fn kmer_set(s: &[u8], k: usize) -> FxHashSet<u64> {
    let mut out = FxHashSet::default();
    if s.len() < k {
        return out;
    }
    let mask = (1u64 << (2 * k)) - 1;
    let mut code = 0u64;
    let mut run = 0usize;
    for &c in s {
        let b = base_code(c);
        if b < 0 {
            run = 0;
            code = 0;
            continue;
        }
        code = ((code << 2) | b as u64) & mask;
        run += 1;
        if run >= k {
            out.insert(code);
        }
    }
    out
}

fn containment(seq: &[u8], kmers: &FxHashSet<u64>, k: usize) -> f64 {
    if seq.len() < k {
        return -1.0;
    }
    let mine = kmer_set(seq, k);
    if mine.is_empty() {
        return -1.0;
    }
    let hit = mine.iter().filter(|k| kmers.contains(k)).count();
    hit as f64 / mine.len() as f64
}

fn compute_locus(seq: &[u8], tr: &TrRegion, loc: &Locus, opts: &Opts) -> PatchResult {
    let mut out = PatchResult::default();
    let stat: &mut LocusStat = &mut out.stat;
    stat.contig = loc.contig.clone();
    stat.start = tr.start;
    stat.end = tr.end;
    stat.period = tr.period;
    stat.old_len = tr.end - tr.start;
    stat.new_len = stat.old_len;
    stat.cont_old = -1.0;
    stat.cont_new = -1.0;
    stat.n_reads = loc.cut_reads.len() as i32;
    let flank = loc.flank;
    let region_len = tr.end - tr.start;

    if opts.max_locus_len > 0 && region_len > opts.max_locus_len {
        stat.status = "skip:array_too_long".into();
        return out;
    }

    stat.span_med = region_len;
    let min_reads = opts.min_reads.max(3);
    if loc.cut_reads.len() < min_reads as usize {
        stat.status = "skip:insufficient_reads".into();
        return out;
    }

    let mut span_med = region_len;
    let mut span_tol = 40i64;
    if !loc.spans.is_empty() {
        let mut sorted = loc.spans.clone();
        sorted.sort_unstable();
        span_med = sorted[sorted.len() / 2];
        stat.span_min = sorted[0];
        stat.span_max = *sorted.last().unwrap();
        let mut dev: Vec<i64> = sorted.iter().map(|&x| (x - span_med).abs()).collect();
        dev.sort_unstable();
        let mad = dev[dev.len() / 2];
        span_tol = 40i64.max((3 * mad).max((0.02 * span_med as f64) as i64));
    }
    stat.span_med = span_med;

    if (region_len - span_med).abs() <= span_tol {
        stat.status = "unchanged_read_length".into();
        return out;
    }

    let mut read_kmers: FxHashSet<u64> = FxHashSet::default();
    for r in &loc.cut_reads {
        let s = kmer_set(r, HIT_KMER);
        read_kmers.extend(s);
    }
    let old_repeat = &seq[tr.start as usize..tr.end as usize];
    stat.cont_old = containment(old_repeat, &read_kmers, HIT_KMER);

    let n_rec = loc.cut_reads.len();
    let total_bases: usize = loc.cut_reads.iter().map(|r| r.len()).sum();
    let mut n_use = n_rec;
    if total_bases > POA_BUDGET {
        n_use = ((n_rec as f64 * POA_BUDGET as f64) / total_bases as f64) as usize;
        n_use = n_use.max(min_reads as usize).min(n_rec);
    }
    let selected: Vec<&[u8]> = (0..n_use)
        .map(|j| {
            let mut idx = (j as f64 * n_rec as f64 / n_use as f64) as usize;
            if idx >= n_rec {
                idx = n_rec - 1;
            }
            loc.cut_reads[idx].as_slice()
        })
        .collect();

    if let Some(dir) = &opts.debug_dir {
        let fn_ = format!(
            "{dir}/locus_{}_{}_{}.reads.fa",
            stat.contig, stat.start, stat.end
        );
        if let Ok(mut f) = File::create(&fn_) {
            for (i, r) in loc.cut_reads.iter().enumerate() {
                let _ = writeln!(f, ">r{i}");
                let _ = f.write_all(r);
                let _ = f.write_all(b"\n");
            }
        }
    }

    let mut config = PoaConfig::default();
    config.min_reads = min_reads as usize;

    let repeat_res = consensus_flanked(&selected, &loc.left, &loc.right, &config);
    let new_repeat = match repeat_res {
        Ok(c) => c.sequence,
        Err(e) => {
            stat.status = format!("skip:{e}");
            return out;
        }
    };
    let whole = consensus(&selected, 0, &config).ok().map(|c| c.sequence);

    if let Some(dir) = &opts.debug_dir {
        let fn_ = format!(
            "{dir}/locus_{}_{}_{}.cons.fa",
            stat.contig, stat.start, stat.end
        );
        if let Ok(mut f) = File::create(&fn_) {
            let _ = writeln!(f, ">repeat_{}", new_repeat.len());
            let _ = f.write_all(&new_repeat);
            let _ = f.write_all(b"\n");
            if let Some(w) = &whole {
                let _ = writeln!(f, ">whole_{}", w.len());
                let _ = f.write_all(w);
                let _ = f.write_all(b"\n");
            }
        }
    }

    let mut anchored = false;
    if let Some(w) = &whole {
        if w.len() >= new_repeat.len() + (1.2 * flank as f64) as usize {
            for probe in [20usize, 16, 12] {
                if probe > flank as usize {
                    continue;
                }
                let pl = find_subslice(w, &loc.left[(flank as usize - probe)..flank as usize]);
                let pr = find_subslice(w, &loc.right[..probe]);
                if let (Some(pl), Some(pr)) = (pl, pr) {
                    if pr > pl {
                        anchored = true;
                        break;
                    }
                }
            }
        }
    }
    if !anchored {
        stat.status = "skip:flank_not_anchored".into();
        return out;
    }

    let ratio = if region_len > 0 {
        new_repeat.len() as f64 / region_len as f64
    } else {
        0.0
    };
    if new_repeat.is_empty() || !(0.4..=2.5).contains(&ratio) {
        stat.status = "skip:length_outlier".into();
        return out;
    }

    stat.cont_new = containment(&new_repeat, &read_kmers, HIT_KMER);
    if stat.cont_new >= 0.0 && stat.cont_old >= 0.0 && stat.cont_new + 0.05 < stat.cont_old {
        stat.status = "skip:lower_support".into();
        return out;
    }

    let delta = new_repeat.len() as i64 - region_len;
    if new_repeat == old_repeat || (delta.abs() <= 2 && stat.cont_new <= stat.cont_old + 0.005) {
        stat.status = "unchanged".into();
        return out;
    }

    let old_err = (region_len - span_med).abs();
    let new_err = (new_repeat.len() as i64 - span_med).abs();
    if old_err <= span_tol || new_err >= old_err {
        stat.status = "unchanged_read_length".into();
        return out;
    }

    out.patch = true;
    out.repeat = new_repeat;
    stat.old_len = region_len;
    stat.new_len = out.repeat.len() as i64;
    stat.status = "patched".into();
    out
}

struct Segment {
    name: String,
    seq: Vec<u8>,
    tags: String,
    a_lines: Vec<String>,
    trs: Vec<TrRegion>,
    ids: Vec<i32>,
}

#[derive(Default)]
struct Counters {
    n_loci: i64,
    n_patched: i64,
    bases_added: i64,
    n_seg: i64,
    n_seg_patched: i64,
    poa_secs: f64,
}

fn process_segment(
    opts: &Opts,
    gfa: &str,
    loci: &[Locus],
    seg: Segment,
    tsv: &mut dyn Write,
    out: &mut dyn Write,
    fa: &mut dyn Write,
    c: &mut Counters,
) -> Result<()> {
    let l = seg.trs.len();
    let mut results: Vec<Option<PatchResult>> = (0..l).map(|_| None).collect();
    let threads = opts.threads.max(1).min(l.max(1));

    if l > 0 {
        let seq = &seg.seq;
        let trs = &seg.trs;
        let ids = &seg.ids;
        if threads > 1 && l > 1 {
            let t0 = Instant::now();
            let locals = std::thread::scope(|scope| {
                let mut handles = Vec::new();
                for ti in 0..threads {
                    handles.push(scope.spawn(move || {
                        let mut local = Vec::new();
                        for k in (ti..l).step_by(threads) {
                            if ids[k] < 0 {
                                continue;
                            }
                            let q0 = Instant::now();
                            let mut pr = compute_locus(seq, &trs[k], &loci[ids[k] as usize], opts);
                            pr.seconds = q0.elapsed().as_secs_f64();
                            local.push((k, pr));
                        }
                        local
                    }));
                }
                let mut all = Vec::new();
                for h in handles {
                    if let Ok(v) = h.join() {
                        all.extend(v);
                    }
                }
                all
            });
            c.poa_secs += t0.elapsed().as_secs_f64();
            for (k, pr) in locals {
                results[k] = Some(pr);
            }
        } else {
            let t0 = Instant::now();
            for k in 0..l {
                if ids[k] < 0 {
                    continue;
                }
                let q0 = Instant::now();
                let mut pr = compute_locus(seq, &trs[k], &loci[ids[k] as usize], opts);
                pr.seconds = q0.elapsed().as_secs_f64();
                results[k] = Some(pr);
            }
            c.poa_secs += t0.elapsed().as_secs_f64();
        }
    }

    // Build the patched sequence once, left to right.
    let mut newseq = Vec::with_capacity(seg.seq.len());
    let mut shifts: Vec<(i64, i64)> = Vec::new();
    let mut prev = 0usize;
    for k in 0..l {
        if let Some(res) = &results[k] {
            if res.patch {
                newseq.extend_from_slice(&seg.seq[prev..seg.trs[k].start as usize]);
                newseq.extend_from_slice(&res.repeat);
                prev = seg.trs[k].end as usize;
                shifts.push((seg.trs[k].end, res.stat.new_len - res.stat.old_len));
            }
        }
    }
    newseq.extend_from_slice(&seg.seq[prev..]);
    let patched = !shifts.is_empty();

    // Report in descending order, matching the original tool.
    for k in (0..l).rev() {
        let Some(res) = &results[k] else { continue };
        let st = &res.stat;
        c.n_loci += 1;
        if res.patch {
            c.n_patched += 1;
            c.bases_added += st.new_len - st.old_len;
        }
        eprintln!(
            "[polish_gfa] {} {}-{} period={} copies={:.1} reads={} support={:.3}->{:.3} med={} [{}..{}] {} ({:.2}s)",
            st.contig,
            st.start,
            st.end,
            st.period,
            (seg.trs[k].end - seg.trs[k].start) as f64 / seg.trs[k].period as f64,
            st.n_reads,
            st.cont_old,
            st.cont_new,
            st.span_med,
            st.span_min,
            st.span_max,
            st.status,
            res.seconds
        );
        writeln!(
            tsv,
            "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}",
            gfa,
            st.contig,
            st.start,
            st.end,
            st.period,
            (seg.trs[k].end - seg.trs[k].start) as f64 / seg.trs[k].period as f64,
            st.old_len,
            st.new_len,
            st.n_reads,
            st.cont_old,
            st.cont_new,
            st.status
        )?;
    }
    if patched {
        c.n_seg_patched += 1;
    }
    c.n_seg += 1;

    // S line, updating LN:i: and preserving all other tags.
    write!(out, "S\t{}\t", seg.name)?;
    out.write_all(&newseq)?;
    let mut kept: Vec<&str> = Vec::new();
    let mut had_ln = false;
    if !seg.tags.is_empty() {
        for t in seg.tags.split('\t') {
            if t.is_empty() {
                continue;
            }
            if t.starts_with("LN:i:") {
                had_ln = true;
                continue;
            }
            kept.push(t);
        }
    }
    if !had_ln && seg.seq != b"*" {
        let ln = format!("LN:i:{}", newseq.len());
        write!(out, "\t{ln}")?;
    }
    for t in kept {
        write!(out, "\t{t}")?;
    }
    writeln!(out)?;

    // A lines with upstream offsets shifted by the splices.
    for a in &seg.a_lines {
        if shifts.is_empty() {
            writeln!(out, "{a}")?;
            continue;
        }
        let mut fields: Vec<String> = a.split('\t').map(|s| s.to_string()).collect();
        if fields.len() >= 3 {
            if let Ok(off) = fields[2].parse::<i64>() {
                let mut shift = 0i64;
                for &(end, delta) in &shifts {
                    if off >= end {
                        shift += delta;
                    }
                }
                fields[2] = (off + shift).to_string();
            }
        }
        writeln!(out, "{}", fields.join("\t"))?;
    }

    if seg.seq != b"*" {
        writeln!(fa, ">{}", seg.name)?;
        for chunk in newseq.chunks(60) {
            fa.write_all(chunk)?;
            fa.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn polish_one_gfa(gfa: &str, opts: &Opts, tsv: &mut dyn Write, c: &mut Counters) -> Result<()> {
    let t_start = Instant::now();
    let mut detector = TrDetector::new();
    let mut loci = collect_loci(gfa, opts, &mut detector)?;
    let t_detect = t_start.elapsed().as_secs_f64();

    let mut t_index_total = 0.0f64;
    let mut t_stream_total = 0.0f64;
    let mut nchunks_done = 0usize;
    let mut stored: u64 = 0;
    let mut last_bases: u64 = 0;

    if !loci.is_empty() {
        let seeds_per_locus = 8 * (opts.edge_len.max(20) as usize - RECRUIT_K + 1);
        let mut chunk = if opts.index_chunk > 0 {
            opts.index_chunk
        } else {
            loci.len()
        };
        if chunk == 0 {
            chunk = loci.len();
        }
        let nchunks = loci.len().div_ceil(chunk);
        let full_gb = loci.len() as f64 * seeds_per_locus as f64 * 72.0 / 1e9;
        let chunk_gb = chunk.min(loci.len()) as f64 * seeds_per_locus as f64 * 72.0 / 1e9;
        if loci.len() > 500_000 || full_gb > 8.0 {
            eprintln!(
                "[polish_gfa] warning: {} loci; full seed index ~{:.1} GB (edge {} bp); {} chunk(s) of <= {} loci -> ~{:.1} GB per chunk, {} read pass(es)",
                loci.len(),
                full_gb,
                opts.edge_len,
                nchunks,
                chunk.min(loci.len()),
                chunk_gb,
                nchunks
            );
            if nchunks == 1 {
                eprintln!(
                    "[polish_gfa] consider --min-array 150/200 or --index-chunk 250000 to bound memory"
                );
            }
        }

        let cap_total = if opts.max_cand_bases > 0 {
            Some(opts.max_cand_bases as u64)
        } else {
            None
        };

        for ci in 0..nchunks {
            let begin = ci * chunk;
            let end = (begin + chunk).min(loci.len());
            let mut idx = SeedIndex::new(0);
            let ti = Instant::now();
            build_index(&loci, &mut idx, begin, end, opts.edge_len.max(20) as usize);
            t_index_total += ti.elapsed().as_secs_f64();

            eprintln!(
                "[polish_gfa] {}: chunk {}/{}: loci {}-{} ({} loci, {} seed k-mers); streaming reads...",
                gfa,
                ci + 1,
                nchunks,
                begin,
                end,
                end - begin,
                idx.head.len()
            );
            let ts = Instant::now();
            let part: RecruitStats = stream_recruit(opts, &mut loci, &idx, cap_total, stored);
            let dt = ts.elapsed().as_secs_f64();
            t_stream_total += dt;
            stored += part.stored_bases;
            last_bases = part.bases;
            nchunks_done += 1;

            let kept: usize = (begin..end).map(|i| loci[i].cut_reads.len()).sum();
            eprintln!(
                "[polish_gfa] {}: chunk {}/{}: {} reads ({} bases), kept {} ({:.2} Mbp), {:.2}s ({:.1} Mbp/s)",
                gfa,
                ci + 1,
                nchunks,
                part.reads,
                part.bases,
                kept,
                part.stored_bases as f64 / 1e6,
                dt,
                if dt > 0.001 {
                    part.bases as f64 / 1e6 / dt
                } else {
                    0.0
                }
            );
            if part.capped {
                eprintln!(
                    "[polish_gfa] candidate cap reached; skipping remaining {} chunk(s)",
                    nchunks - ci - 1
                );
                break;
            }
        }

        let kept_all: usize = loci.iter().map(|l| l.cut_reads.len()).sum();
        eprintln!(
            "[polish_gfa] {}: {} loci, {} chunk(s), kept {} reads ({:.2} Mbp); detect {:.2}s, index {:.2}s, stream {:.2}s ({:.1} Mbp/s aggregate)",
            gfa,
            loci.len(),
            nchunks_done,
            kept_all,
            stored as f64 / 1e6,
            t_detect,
            t_index_total,
            t_stream_total,
            if t_stream_total > 0.001 && nchunks_done > 0 {
                last_bases as f64 * nchunks_done as f64 / 1e6 / t_stream_total
            } else {
                0.0
            }
        );
    }

    // Pass 2: patch and write.
    let seg_before = c.n_seg;
    let patched_before = c.n_seg_patched;
    let poa_before = c.poa_secs;
    let base = strip_gfa_suffix(gfa);
    let out_path = format!("{base}.poa.gfa");
    let fa_path = format!("{base}.poa.fa");
    let in_f = File::open(gfa).with_context(|| format!("open {gfa}"))?;
    let reader = BufReader::with_capacity(1 << 20, in_f);
    let mut out = BufWriter::new(File::create(&out_path)?);
    let mut fa = BufWriter::new(File::create(&fa_path)?);

    let mut next_locus = 0usize;
    let mut current: Option<Segment> = None;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("S\t") {
            if let Some(seg) = current.take() {
                process_segment(opts, gfa, &loci, seg, tsv, &mut out, &mut fa, c)?;
            }
            let Some((name, seq, tags)) = parse_s_line(&line) else {
                continue;
            };
            let mut trs = Vec::new();
            let mut ids = Vec::new();
            if seq != b"*" {
                trs = detector.detect(&seq, opts.max_period, opts.min_array);
                for tr in &trs {
                    if opts.max_locus_len > 0 && tr.end - tr.start > opts.max_locus_len {
                        eprintln!(
                            "[polish_gfa] {} {}-{} array is {} bp (> --max-locus-len {}); skipped",
                            name,
                            tr.start,
                            tr.end,
                            tr.end - tr.start,
                            opts.max_locus_len
                        );
                        ids.push(-1);
                        continue;
                    }
                    let flank = locus_flank(seq.len() as i64, tr, opts);
                    if flank < 40 {
                        ids.push(-1);
                        continue;
                    }
                    if next_locus >= loci.len() {
                        eprintln!(
                            "[polish_gfa] warning: locus mismatch in {gfa} at {name}:{}",
                            tr.start
                        );
                        ids.push(-1);
                        continue;
                    }
                    let loc = &loci[next_locus];
                    if loc.contig != name || loc.start != tr.start || loc.end != tr.end {
                        eprintln!(
                            "[polish_gfa] warning: locus mismatch in {gfa} at {name}:{}",
                            tr.start
                        );
                        ids.push(-1);
                        next_locus += 1;
                        continue;
                    }
                    ids.push(next_locus as i32);
                    next_locus += 1;
                }
            }
            current = Some(Segment {
                name,
                seq,
                tags,
                a_lines: Vec::new(),
                trs,
                ids,
            });
        } else if line.starts_with("A\t") {
            if let Some(seg) = current.as_mut() {
                seg.a_lines.push(line);
            }
        } else {
            if let Some(seg) = current.take() {
                process_segment(opts, gfa, &loci, seg, tsv, &mut out, &mut fa, c)?;
            }
            writeln!(out, "{line}")?;
        }
    }
    if let Some(seg) = current.take() {
        process_segment(opts, gfa, &loci, seg, tsv, &mut out, &mut fa, c)?;
    }

    if next_locus != loci.len() {
        eprintln!(
            "[polish_gfa] warning: {gfa}: used {} of {} loci",
            next_locus,
            loci.len()
        );
    }
    out.flush()?;
    fa.flush()?;

    let t_end = t_start.elapsed().as_secs_f64();
    eprintln!(
        "[polish_gfa] {}: {} contig(s), {} polished; timing: detect {:.2}s, index {:.2}s, stream {:.2}s, patch {:.2}s (poa {:.2}s), total {:.2}s",
        gfa,
        c.n_seg - seg_before,
        c.n_seg_patched - patched_before,
        t_detect,
        t_index_total,
        t_stream_total,
        t_end - t_detect - t_index_total - t_stream_total,
        c.poa_secs - poa_before,
        t_end
    );
    Ok(())
}

pub fn run(opts: &Opts) -> Result<i32> {
    let cap = if opts.max_cand_bases > 0 {
        format!("{:.2} Gbp", opts.max_cand_bases as f64 / 1e9)
    } else {
        "unlimited".to_string()
    };
    let backend = match opts.decompressor {
        crate::model::Decompressor::Auto => "rapidgzip-core (auto)",
        crate::model::Decompressor::Zlib => "zlib",
        crate::model::Decompressor::Rapidgzip => "rapidgzip-core",
    };
    eprintln!("[polish_gfa] read backend: {backend}; plain inputs bypass it");
    eprintln!(
        "[polish_gfa] standalone poa-consensus tandem-repeat polishing (v{}); streaming recruitment (flank {} bp, edge {} bp, seed step {}, candidate cap {}, min array {}, max locus {})",
        env!("CARGO_PKG_VERSION"),
        opts.flank_len,
        opts.edge_len,
        opts.seed_step,
        cap,
        opts.min_array,
        opts.max_locus_len
    );

    for g in &opts.gfas {
        if !Path::new(g).is_file() {
            bail!("GFA not found: {g}");
        }
    }
    let tsv_path = opts
        .tsv
        .clone()
        .unwrap_or_else(|| format!("{}.poa.tsv", strip_gfa_suffix(&opts.gfas[0])));
    let mut tsv =
        BufWriter::new(File::create(&tsv_path).with_context(|| format!("create {tsv_path}"))?);
    writeln!(
        tsv,
        "#gfa\tcontig\tstart\tend\tperiod\tcopies\told_len\tnew_len\treads\tsupport_old\tsupport_new\tstatus"
    )?;

    let mut c = Counters::default();
    for gfa in &opts.gfas {
        if let Err(e) = polish_one_gfa(gfa, opts, &mut tsv, &mut c) {
            eprintln!("[polish_gfa] {gfa}: {e:#}");
        }
    }
    tsv.flush()?;

    eprintln!(
        "[polish_gfa] done: {} locus/loci found, {} patched, {:+} bp; report: {}",
        c.n_loci, c.n_patched, c.bases_added, tsv_path
    );
    Ok(0)
}
