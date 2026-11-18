use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use std::process::{Child, Command, Stdio};

use flate2::read::MultiGzDecoder;

use crate::model::Decompressor;

pub fn command_available(name: &str) -> bool {
    let Some(path) = env::var_os("PATH") else {
        return false;
    };
    env::split_paths(&path).any(|dir| {
        let p = dir.join(name);
        p.is_file() && is_executable(&p)
    })
}

#[cfg(unix)]
fn is_executable(p: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    p.metadata()
        .map(|m| m.permissions().mode() & 0o111 != 0)
        .unwrap_or(false)
}

#[cfg(not(unix))]
fn is_executable(p: &Path) -> bool {
    p.is_file()
}

pub fn path_is_gzip(path: &str) -> bool {
    if path.ends_with(".gz") || path.ends_with(".bgz") {
        return true;
    }
    let Ok(mut f) = File::open(path) else {
        return false;
    };
    let mut m = [0u8; 2];
    f.read_exact(&mut m).is_ok() && m == [0x1f, 0x8b]
}

pub struct ReadSource {
    reader: Option<Box<dyn BufRead + Send>>,
    child: Option<Child>,
}

impl ReadSource {
    fn read_until(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        self.reader.as_mut().unwrap().read_until(b'\n', buf)
    }
}

impl Drop for ReadSource {
    fn drop(&mut self) {
        // Drop the pipe/stdout before waiting, or a blocked writer deadlocks.
        self.reader.take();
        if let Some(mut child) = self.child.take() {
            match child.wait() {
                Ok(st) if !st.success() => {
                    eprintln!("[polish_gfa] warning: rapidgzip exited with {st}");
                }
                _ => {}
            }
        }
    }
}

pub fn open_source(
    path: &str,
    decompressor: Decompressor,
    gzip_threads: usize,
) -> io::Result<ReadSource> {
    if !path_is_gzip(path) {
        let f = File::open(path)?;
        return Ok(ReadSource {
            reader: Some(Box::new(BufReader::with_capacity(1 << 20, f))),
            child: None,
        });
    }

    let want_rapidgzip = decompressor == Decompressor::Rapidgzip
        || (decompressor == Decompressor::Auto && command_available("rapidgzip"));
    if want_rapidgzip {
        match Command::new("rapidgzip")
            .args(["-d", "-c", "-P", &gzip_threads.max(1).to_string(), path])
            .stdout(Stdio::piped())
            .spawn()
        {
            Ok(mut child) => {
                let stdout = child.stdout.take().expect("piped stdout");
                return Ok(ReadSource {
                    reader: Some(Box::new(BufReader::with_capacity(1 << 20, stdout))),
                    child: Some(child),
                });
            }
            Err(e) => {
                if decompressor == Decompressor::Rapidgzip {
                    eprintln!("[polish_gfa] warning: could not start rapidgzip ({e}); using zlib");
                }
            }
        }
    }

    let f = File::open(path)?;
    let dec = MultiGzDecoder::new(f);
    Ok(ReadSource {
        reader: Some(Box::new(BufReader::with_capacity(1 << 20, dec))),
        child: None,
    })
}

pub struct FastxReader {
    src: ReadSource,
    pending: Option<Vec<u8>>,
}

impl FastxReader {
    pub fn new(src: ReadSource) -> Self {
        Self { src, pending: None }
    }

    fn read_raw(&mut self, out: &mut Vec<u8>) -> io::Result<bool> {
        out.clear();
        let n = self.src.read_until(out)?;
        if n == 0 {
            return Ok(false);
        }
        while out.last().is_some_and(|&c| c == b'\n' || c == b'\r') {
            out.pop();
        }
        Ok(true)
    }

    /// Returns the next sequence, handling FASTA (wrapped) and FASTQ.
    pub fn next_seq(&mut self) -> io::Result<Option<Vec<u8>>> {
        loop {
            let header = match self.pending.take() {
                Some(h) => h,
                None => {
                    let mut line = Vec::new();
                    loop {
                        if !self.read_raw(&mut line)? {
                            return Ok(None);
                        }
                        if !line.is_empty() {
                            break;
                        }
                    }
                    line
                }
            };

            if header.first() == Some(&b'>') {
                let mut seq = Vec::new();
                let mut line = Vec::new();
                loop {
                    if !self.read_raw(&mut line)? {
                        break;
                    }
                    if line.first() == Some(&b'>') {
                        self.pending = Some(line);
                        break;
                    }
                    if line.first() == Some(&b';') {
                        continue;
                    }
                    seq.extend_from_slice(&line);
                }
                if !seq.is_empty() {
                    return Ok(Some(seq));
                }
            } else if header.first() == Some(&b'@') {
                let mut seq = Vec::new();
                let mut line = Vec::new();
                loop {
                    if !self.read_raw(&mut line)? {
                        break;
                    }
                    if line.first() == Some(&b'+') {
                        break;
                    }
                    seq.extend_from_slice(&line);
                }
                let mut qlen = 0usize;
                while qlen < seq.len() {
                    if !self.read_raw(&mut line)? {
                        break;
                    }
                    qlen += line.len();
                }
                if !seq.is_empty() {
                    return Ok(Some(seq));
                }
            }
            // Unknown header: skip to the next record.
        }
    }
}
