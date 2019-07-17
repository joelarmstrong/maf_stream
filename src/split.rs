use multiple_alignment_format::parser::next_maf_item;
use multiple_alignment_format::{MAFBlock, MAFItem};
use std::io::{BufRead, Write, BufWriter};
use std::path::PathBuf;
use std::fs::File;
use itertools::Itertools;

struct MAFSplit {
    /// Chromosome of reference within current file.
    cur_chrom: Option<String>,
    /// Total length of blocks (in reference coordinates) within current file.
    cur_length: Option<u64>,
    cur_file: Option<BufWriter<File>>,
    output_dir: PathBuf,
    /// Maximum aligned length (in reference) per file.
    max_length: u64,
}

impl MAFSplit {
    fn new(output_dir: &str, max_length: u64) -> MAFSplit {
        Self {
            cur_chrom: None,
            cur_length: None,
            cur_file: None,
            output_dir: PathBuf::from(output_dir),
            max_length
        }
    }

    /// Outputs this block to the correct file, opening a new one if
    /// needed.
    fn output_block(&mut self, block: &MAFBlock) {
        let ref_line = block.aligned_entries().next();
        if let Some(ref_aln) = ref_line {
            let chr = ref_aln.seq.split(".").skip(1).join(".");
            // On any new reference chromosome, or if the file would grow too
            // large, we switch to a new file.
            if self.cur_chrom.is_none() ||
                self.cur_length.is_none() ||
                &chr != self.cur_chrom.as_ref().unwrap() ||
                self.cur_length.unwrap() + ref_aln.aligned_length > self.max_length {
                    self.new_file(&chr, ref_aln.start);
            }
            self.cur_length = self.cur_length.map(|l| l + ref_aln.aligned_length);
        }
        write!(self.cur_file.as_mut().unwrap(), "{}", block).expect("failed to write");
    }

    /// Starts a new file and flushes the old one.
    fn new_file(&mut self, chrom: &str, start_pos: u64) {
        let f = File::create(self.output_dir.join(format!("{}.{}.maf", chrom, start_pos))).expect("Couldn't create file");
        self.cur_file = Some(BufWriter::new(f));
        self.cur_length = Some(0);
        self.cur_chrom = Some(chrom.to_string());
        writeln!(self.cur_file.as_mut().unwrap(), "##maf version=1").expect("failed to write");
    }
}

pub fn split_maf(input: &mut dyn BufRead, max_length: u64, output_dir: &str) {
    let mut splitter = MAFSplit::new(output_dir, max_length);

    while let Ok(item) = next_maf_item(input) {
        match item {
            MAFItem::Block(block) => {
                splitter.output_block(&block);
            },
            _ => {},
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;
    use std::path::Path;
    use std::fs::read_to_string;

    #[test]
    fn test_simple_split() {
        let input_maf = "##maf version=1
a
s       Rhesus.chr21_chr20      0       54      +       19571763        AATTCTGTGAAGCTTCTTTGAGAGGCTTGGATTTATTTCACACATTCGAACATT
s       Human.chr21     0       54      +       9688985 AGTTCTGAGAAGCTTCTTTGTGAGGCTTGGATTCATTTCACACATTTGAACAtt

a
s       Rhesus.chr21_chr20      54      28      +       19571763        TGATTGAAGATTTGGAAACAGTCTTTTT
s       Human.chr21     58      27      +       9688985 tgattgtagatctggaaacagtctt-tt

a
s       Rhesus.chr21_chr20      82      16      +       19571763        TGTAAAATCTATAAAG
s       Human.chr21     85      16      +       9688985 tgtgaaatctataaag

a
s       Rhesus.chr22      193     32      +       19571763        aacctttcctttgctagagcactttggaaata
s       Human.chr21     217     32      +       9688985 aacctttcctttgctagagcactttgaaaata
";
        let tempdir = TempDir::new().unwrap();
        let output_dir = tempdir.path().to_str().unwrap();
        split_maf(&mut input_maf.as_bytes(), 84, &output_dir);

        // The first two blocks should fit in one file, the third
        // should spill over into another file, and the fourth should
        // go into yet another file (since it switches reference
        // chromosome).
        assert!(Path::exists(&tempdir.path().join("chr21_chr20.0.maf")));
        assert!(!Path::exists(&tempdir.path().join("chr21_chr20.54.maf")));
        assert!(Path::exists(&tempdir.path().join("chr21_chr20.82.maf")));
        assert!(Path::exists(&tempdir.path().join("chr22.193.maf")));

        assert_eq!(read_to_string(&tempdir.path().join("chr21_chr20.0.maf")).unwrap(),
                                  "##maf version=1
a
s Rhesus.chr21_chr20 0 54 + 19571763 AATTCTGTGAAGCTTCTTTGAGAGGCTTGGATTTATTTCACACATTCGAACATT
s Human.chr21 0 54 + 9688985 AGTTCTGAGAAGCTTCTTTGTGAGGCTTGGATTCATTTCACACATTTGAACAtt

a
s Rhesus.chr21_chr20 54 28 + 19571763 TGATTGAAGATTTGGAAACAGTCTTTTT
s Human.chr21 58 27 + 9688985 tgattgtagatctggaaacagtctt-tt

");

        assert_eq!(read_to_string(&tempdir.path().join("chr21_chr20.82.maf")).unwrap(),
                                  "##maf version=1
a
s Rhesus.chr21_chr20 82 16 + 19571763 TGTAAAATCTATAAAG
s Human.chr21 85 16 + 9688985 tgtgaaatctataaag

");

        assert_eq!(read_to_string(&tempdir.path().join("chr22.193.maf")).unwrap(),
                                  "##maf version=1
a
s Rhesus.chr22 193 32 + 19571763 aacctttcctttgctagagcactttggaaata
s Human.chr21 217 32 + 9688985 aacctttcctttgctagagcactttgaaaata

");
    }
}
