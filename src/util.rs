pub const RECRUIT_K: usize = 17;
pub const HIT_KMER: usize = 15;
pub const MAX_READS_PER_LOCUS: usize = 128;

#[inline]
pub fn base_code(c: u8) -> i32 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => -1,
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

pub fn find_subslice(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    let first = needle[0];
    let mut start = 0usize;
    while start + needle.len() <= hay.len() {
        let p = memchr::memchr(first, &hay[start..])?;
        let i = start + p;
        if i + needle.len() > hay.len() {
            return None;
        }
        if hay[i..].starts_with(needle) {
            return Some(i);
        }
        start = i + 1;
    }
    None
}

pub fn rfind_subslice(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(hay.len());
    }
    let mut end = hay.len();
    let first = needle[0];
    while end >= needle.len() {
        let p = memchr::memrchr(first, &hay[..end])?;
        if p + needle.len() <= hay.len() && hay[p..].starts_with(needle) {
            return Some(p);
        }
        if p == 0 {
            return None;
        }
        end = p;
    }
    None
}

pub fn median(v: &mut [i64]) -> i64 {
    if v.is_empty() {
        return 0;
    }
    v.sort_unstable();
    v[v.len() / 2]
}

pub fn strip_gfa_suffix(fn_: &str) -> &str {
    for suf in [".gfa.gz", ".gfa.bgz", ".gfa", ".bgz", ".gz"] {
        if let Some(stem) = fn_.strip_suffix(suf) {
            return stem;
        }
    }
    fn_
}

/// Peak resident set size of this process, in bytes.
#[cfg(target_os = "linux")]
fn peak_rss_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            let kib: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kib * 1024);
        }
    }
    None
}

#[cfg(all(unix, not(target_os = "linux")))]
fn peak_rss_bytes() -> Option<u64> {
    unsafe {
        let mut ru: libc::rusage = std::mem::zeroed();
        if libc::getrusage(libc::RUSAGE_SELF, &mut ru) == 0 {
            let n = ru.ru_maxrss as u64;
            #[cfg(target_os = "macos")]
            {
                Some(n)
            }
            #[cfg(all(unix, not(target_os = "macos"), not(target_os = "linux")))]
            {
                Some(n * 1024)
            }
        } else {
            None
        }
    }
}

#[cfg(not(unix))]
fn peak_rss_bytes() -> Option<u64> {
    None
}

/// Print the process peak RSS. Intended as the last log line.
pub fn log_peak_rss() {
    if let Some(bytes) = peak_rss_bytes() {
        let gib = bytes as f64 / (1u64 << 30) as f64;
        let mib = bytes / (1u64 << 20);
        eprintln!("[polish_gfa] peak RSS: {gib:.2} GiB ({mib} MiB)");
    }
}
