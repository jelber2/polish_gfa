use crate::model::TrRegion;

pub const TR_KMER: usize = 11;
const MIN_COPIES: f64 = 2.5;

pub struct TrDetector {
    last: Vec<i32>,
}

impl TrDetector {
    pub fn new() -> Self {
        Self {
            last: vec![-1; 1usize << (2 * TR_KMER)],
        }
    }

    #[inline]
    fn code(c: u8) -> i32 {
        match c {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => -1,
        }
    }

    fn periodic_at(&self, s: &[u8], i: i64, p: i32, p_used: &mut i32) -> bool {
        let l = s.len() as i64;
        let mut best = 0.0f64;
        let mut best_p = p;
        for d in -3..=3 {
            let pp = p + d;
            if pp < 2 || i + pp as i64 >= l {
                continue;
            }
            let mut m = 0i32;
            let mut n = 0i32;
            for t in 0..12i64 {
                if i + t + pp as i64 >= l {
                    break;
                }
                m += (s[(i + t) as usize] == s[(i + t + pp as i64) as usize]) as i32;
                n += 1;
            }
            if n < 6 {
                continue;
            }
            let f = m as f64 / n as f64;
            if f > best {
                best = f;
                best_p = pp;
            }
        }
        *p_used = best_p;
        best >= 0.60
    }

    pub fn detect(&mut self, s: &[u8], max_period: i32, min_len: i32) -> Vec<TrRegion> {
        let mut regions = Vec::new();
        let l = s.len() as i64;
        if l < min_len as i64 {
            return regions;
        }

        self.last.fill(-1);
        let tbl = 1usize << (2 * TR_KMER);
        let mask = (tbl - 1) as u32;
        let mut code = 0u32;
        let mut run = 0usize;
        let mut next_allowed = 0i64;

        for i in 0..l {
            let bc = Self::code(s[i as usize]);
            if bc < 0 {
                run = 0;
                code = 0;
                continue;
            }
            code = ((code << 2) | bc as u32) & mask;
            run += 1;
            if run < TR_KMER {
                continue;
            }

            let prev = self.last[code as usize];
            self.last[code as usize] = i as i32;
            if prev < 0 || i < next_allowed {
                continue;
            }

            let p = i - prev as i64;
            if p < 2 || p > max_period as i64 {
                continue;
            }
            let p = p as i32;

            let w = std::cmp::min(2 * p + 50, 500);
            let mut m = 0i32;
            let mut n = 0i32;
            for t in 0..w as i64 {
                if i + t + p as i64 >= l {
                    break;
                }
                m += (s[(i + t) as usize] == s[(i + t + p as i64) as usize]) as i32;
                n += 1;
            }
            if n < std::cmp::min(2 * p, w - 5) || m * 100 < 72 * n {
                continue;
            }

            let mut a = i;
            let mut pa = p;
            while a > 0 {
                let mut np = pa;
                if !self.periodic_at(s, a - 1, pa, &mut np) {
                    break;
                }
                pa = np;
                a -= 1;
            }
            let mut b = i;
            let mut pb = p;
            while b + (pb as i64) < l {
                let mut np = pb;
                if !self.periodic_at(s, b, pb, &mut np) {
                    break;
                }
                pb = np;
                b += 1;
            }
            let e = std::cmp::min(l, b + pb as i64);

            if e - a >= min_len as i64 && (e - a) as f64 >= MIN_COPIES * p as f64 {
                regions.push(TrRegion {
                    start: a,
                    end: e,
                    period: p,
                });
                next_allowed = e;
            }
        }

        regions.sort_by(|x, y| {
            if x.start != y.start {
                x.start.cmp(&y.start)
            } else {
                (y.end - y.start).cmp(&(x.end - x.start))
            }
        });
        let mut merged: Vec<TrRegion> = Vec::new();
        for r in regions {
            if let Some(m) = merged.last_mut() {
                if r.start <= m.end {
                    let ml = m.end - m.start;
                    let rl = r.end - r.start;
                    if (r.period - m.period).abs() <= 2 {
                        m.start = m.start.min(r.start);
                        m.end = m.end.max(r.end);
                        if rl > ml {
                            m.period = r.period;
                        }
                    } else if rl > ml {
                        *m = r;
                    }
                } else if r.start - m.end <= 2 * std::cmp::max(r.period, m.period) as i64 + 30
                    && (r.period - m.period).abs() <= 2
                {
                    m.end = r.end;
                } else {
                    merged.push(r);
                }
            } else {
                merged.push(r);
            }
        }
        merged
    }
}
