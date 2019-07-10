use multiple_alignment_format::parser::next_maf_item;
use multiple_alignment_format::{MAFBlock, MAFItem, MAFBlockEntry, MAFBlockAlignedEntry};
use std::collections::HashMap;
use std::io::{BufRead, Read, Write, BufWriter, Seek, SeekFrom};
use std::fs::File;

use tempfile::tempfile;

const DEFAULT_FASTA_WIDTH: u32 = 120;

/// Stores the aligned sequence for each genome, storing them in
/// temporary files when necessary, then spits out an aligned FASTA
/// when done.
struct GenomeFastaAccumulator {
    files: HashMap<String, BufWriter<File>>,
    /// Current position within the reference sequence.
    cur_pos: u64,
    /// Max width of the aligned FASTA entries.
    fasta_width: u32,
    ref_name: String,
}

impl GenomeFastaAccumulator {
    fn new(ref_name: String) -> Self {
        Self {
            files: HashMap::new(),
            cur_pos: 0,
            fasta_width: DEFAULT_FASTA_WIDTH,
            ref_name,
        }
    }

    fn push(&mut self, seq: &str, chars: &[u8], ref_pos: u64) {
        if seq == self.ref_name {
            if ref_pos != self.cur_pos + 1 {
                panic!("Ref pos skipped from {} to {}", self.cur_pos, ref_pos);
            }
            self.cur_pos == ref_pos;
        }
        let file_opt = self.files.get_mut(seq);
        let file = match file_opt {
            None => {
                self.files.insert(seq.to_string(), BufWriter::new(tempfile().expect("Couldn't open temporary file")));
                let file = self.files.get_mut(seq).unwrap();
                for _ in 0..ref_pos {
                    file.write(&[b'-']);
                }
                file
            },
            Some(file) => file,
        };
        file.write(chars);
    }

    fn write_fasta(&mut self, output: &mut dyn Write) {
        let mut files: Vec<_> = self.files.values().collect();
        for file in files.iter_mut() {
            file.flush();
            let inner = file.get_mut();
            inner.seek(SeekFrom::Start(0));
            let line_buf = String::new();
            inner.read_to_string(&mut line_buf);
        }
    }
}

pub fn maf_to_fasta(input: &mut dyn BufRead, output: &mut dyn Write) {
    while let Ok(item) = next_maf_item(input) {
        match item {
            MAFItem::Block(block) => {
                
            },
            _ => {},
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}
