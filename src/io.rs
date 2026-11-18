use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

use flate2::read::MultiGzDecoder;
use rapidgzip_core::Decoder;

use crate::model::Decompressor;

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
    reader: Box<dyn BufRead + Send>,
}

impl ReadSource {
    fn read_until(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        self.reader.read_until(b'\n', buf)
    }
}

/// Open a read file. `rapidgzip` uses the built-in `rapidgzip-core` crate
/// (parallel DEFLATE decoding); `zlib` uses flate2. `auto` prefers
/// rapidgzip-core. Plain files bypass both.
pub fn open_source(
    path: &str,
    decompressor: Decompressor,
    gzip_threads: usize,
) -> io::Result<ReadSource> {
    if !path_is_gzip(path) {
        let f = File::open(path)?;
        return Ok(ReadSource {
            reader: Box::new(BufReader::with_capacity(1 << 20, f)),
        });
    }

    if decompressor != Decompressor::Zlib {
        let opened = match Decoder::builder()
            .decoder_threads(gzip_threads.max(1))
            .build()
        {
            Ok(decoder) => decoder.open(path).map_err(io::Error::other),
            Err(e) => Err(io::Error::other(e)),
        };
        match opened {
            Ok(reader) => {
                return Ok(ReadSource {
                    reader: Box::new(BufReader::with_capacity(1 << 20, reader)),
                });
            }
            Err(e) => {
                eprintln!(
                    "[polish_gfa] warning: rapidgzip-core failed for {path} ({e}); using zlib"
                );
            }
        }
    }

    let f = File::open(path)?;
    let dec = MultiGzDecoder::new(f);
    Ok(ReadSource {
        reader: Box::new(BufReader::with_capacity(1 << 20, dec)),
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
