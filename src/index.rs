use rustc_hash::{FxHashMap, FxHashSet};

use crate::model::Locus;
use crate::util::{revcomp, RECRUIT_K};

pub const NONE: u32 = u32::MAX;

#[derive(Clone, Copy)]
pub struct SeedHit {
    pub packed: u64,
    pub next: u32,
}

pub struct SeedIndex {
    pub head: FxHashMap<u64, u32>,
    pub hits: Vec<SeedHit>,
}

impl SeedIndex {
    pub fn new(expected: usize) -> Self {
        Self {
            head: FxHashMap::with_capacity_and_hasher(expected, Default::default()),
            hits: Vec::with_capacity(expected),
        }
    }

    #[inline]
    pub fn add(&mut self, kmer: u64, locus: u64, strand: u8, side: u8, offset: i32) {
        let h = self.head.entry(kmer).or_insert(NONE);
        let idx = self.hits.len() as u32;
        self.hits.push(SeedHit {
            packed: (locus << 15)
                | ((strand as u64) << 14)
                | ((side as u64) << 12)
                | (offset as u64 & 0xFFF),
            next: *h,
        });
        *h = idx;
    }
}

pub fn flank_edge(loc: &Locus, side: u8, edge_len: usize) -> Vec<u8> {
    let f = loc.flank as usize;
    let e = edge_len.min(f);
    let src = match side {
        0 => &loc.left,
        1 => &loc.left,
        2 => &loc.right,
        _ => &loc.right,
    };
    match side {
        0 | 2 => src[..e].to_vec(),
        _ => src[f - e..f].to_vec(),
    }
}

pub fn low_complexity_kmer(code: u64, k: usize) -> bool {
    let mut prev: i32 = -1;
    let mut run = 1;
    let mut maxrun = 1;
    let mut seen = 0u32;
    let mut distinct = 0;
    for j in 0..k {
        let b = ((code >> (2 * (k - 1 - j))) & 3) as i32;
        if seen & (1u32 << b) == 0 {
            seen |= 1u32 << b;
            distinct += 1;
        }
        if b == prev {
            run += 1;
            maxrun = maxrun.max(run);
        } else {
            run = 1;
            prev = b;
        }
    }
    distinct <= 2 || maxrun >= 10
}

pub fn build_index(loci: &[Locus], idx: &mut SeedIndex, begin: usize, end: usize, edge_len: usize) {
    let n = end.saturating_sub(begin);
    let seeds_per_locus = 8 * (edge_len.saturating_sub(RECRUIT_K) + 1);
    let expected = n * seeds_per_locus + 16;
    *idx = SeedIndex::new(expected);

    let mut seen: FxHashSet<u64> = FxHashSet::default();
    for (li, loc) in loci[begin..end].iter().enumerate() {
        let locus = (begin + li) as u64;
        for side in 0..4u8 {
            let edge = flank_edge(loc, side, edge_len);
            for strand in 0..2u8 {
                let s = if strand == 0 {
                    edge.clone()
                } else {
                    revcomp(&edge)
                };
                if s.len() < RECRUIT_K {
                    continue;
                }
                seen.clear();
                let mut code = 0u64;
                let mut run = 0usize;
                for (p, &ch) in s.iter().enumerate() {
                    let b = crate::util::base_code(ch);
                    if b < 0 {
                        run = 0;
                        code = 0;
                        continue;
                    }
                    code = ((code << 2) | b as u64) & ((1u64 << (2 * RECRUIT_K)) - 1);
                    run += 1;
                    if run < RECRUIT_K {
                        continue;
                    }
                    if low_complexity_kmer(code, RECRUIT_K) {
                        continue;
                    }
                    if seen.insert(code) {
                        idx.add(code, locus, strand, side, (p + 1 - RECRUIT_K) as i32);
                    }
                }
            }
        }
    }
}

