use multiple_alignment_format::parser::next_maf_item;
use multiple_alignment_format::{MAFBlock, MAFItem, MAFBlockEntry, MAFBlockAlignedEntry};
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, Write};

#[derive(PartialEq, Eq, Debug)]
struct Range {
    seq: String,
    start: u64,
    end: u64,
}

impl Ord for Range {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.seq, &self.start, Reverse(&self.end)).cmp(&(&other.seq, &other.start, Reverse(&other.end)))
    }
}

impl PartialOrd for Range {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct MAFCoverage {
    /// Coverage by genome.
    coverage: HashMap<String, u64>,
    ranges: Option<BTreeSet<Range>>,
    ref_genome: String,
    ref_length: u64,
}

fn aligned_base(base: u8) -> bool {
    match base {
        b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n' => true,
        _ => false,
    }
}

impl MAFCoverage {
    fn new(ref_genome: &str, ranges: Option<BTreeSet<Range>>) -> Self {
        MAFCoverage {
            coverage: HashMap::new(),
            ref_genome: ref_genome.to_string(),
            ranges,
            ref_length: 0,
        }
    }

    fn add_block(&mut self, block: MAFBlock) {
        let entries = block.entries_as_hash();
        let ref_entries_opt = entries.get::<str>(&self.ref_genome);
        if let Some(ref_entries) = ref_entries_opt {
            for ref_entry in ref_entries {
                self.add_block_with_ref_entry(ref_entry, &entries);
            }
        }
    }

    fn add_block_with_ref_entry(&mut self, ref_entry: &MAFBlockAlignedEntry, entries: &HashMap<&str, Vec<&MAFBlockAlignedEntry>>) {
        for i in 0..ref_entry.alignment.len() {
            // Within each column, add a base of coverage to a genome
            // if at least one entry in the genome is aligned (not a
            // gap) and the reference is aligned (not a gap).
            if !aligned_base(ref_entry.alignment[i]) {
                continue;
            }
            for (genome, genome_entries) in entries {
                let mut found_alignment = false;
                for genome_entry in genome_entries {
                    if aligned_base(genome_entry.alignment[i]) {
                        found_alignment = true;
                        break;
                    }
                }
                if found_alignment {
                    if !self.coverage.contains_key(*genome) {
                        self.coverage.insert(genome.to_string(), 0);
                    }
                    let coverage = self.coverage.get_mut(*genome).unwrap();
                    *coverage += 1;
                }
            }
        }
        self.ref_length += ref_entry.aligned_length;
    }

    fn print(&self, output: &mut dyn Write) {
        writeln!(output, "# referenceSpecies/Chr\tquerySpecies/Chr\tlengthOfReference\tpercentCoverage\tbasesCoverage");
        let total = match &self.ranges {
            None => self.ref_length,
            Some(set) => set.iter().map(|p| p.end - p.start).sum(),
        };
        for (genome, coverage) in self.coverage.iter() {
            writeln!(output, "{}\t{}\t{}\t{}\t{}", self.ref_genome, genome,
                     total, (*coverage as f64) / (total as f64), coverage);
        }
    }

    fn in_range(&self, chrom: &str, position: u64) -> bool {
        match &self.ranges {
            None => true,
            Some(set) => {
                let pos = Range {
                    seq: chrom.to_string(),
                    start: position,
                    end: position + 1,
                };
                dbg!(set.range(..=&pos).next_back());
                match set.range(..=pos).next_back() {
                    None => false,
                    Some(range) => range.seq == chrom && range.start <= position && range.end > position,
                }
            },
        }
    }
}

fn parse_bed(bed: impl BufRead) -> BTreeSet<Range> {
    bed.lines()
       .map(|line_res| {
           let line = line_res.expect("Can't read line");
           let fields: Vec<_> = line.split_whitespace().collect();
           if fields.len() > 9 {
               panic!("BED12 input not supported");
           }
           let seq = fields[0].to_string();
           let start: u64 = fields[1].parse().expect("Can't parse start position");
           let end: u64 = fields[2].parse().expect("Can't parse start position");
           Range { seq, start, end }
       }).collect()
}

pub fn coverage(input: &mut dyn BufRead, output: &mut dyn Write,
                ref_genome: &str, bed: Option<impl BufRead>) {
    let ranges = bed.map(parse_bed);

    let mut maf_coverage = MAFCoverage::new(ref_genome, ranges);

    while let Ok(item) = next_maf_item(input) {
        match item {
            MAFItem::Block(block) => {
                maf_coverage.add_block(block);
            },
            _ => {},
        }
    }

    maf_coverage.print(output);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_in_range() {
        let ranges: BTreeSet<_> = vec![
            Range {
                seq: "chr1".to_string(),
                start: 20,
                end: 30,
            },
            Range {
                seq: "chr1".to_string(),
                start: 30,
                end: 31,
            },
            Range {
                seq: "chr1".to_string(),
                start: 32,
                end: 40,
            },
            Range {
                seq: "chrW".to_string(),
                start: 0,
                end: 300
            }
        ].into_iter()
         .collect();
        let maf_coverage = MAFCoverage::new("none", Some(ranges));
        assert!(!maf_coverage.in_range("chr0", 0));
        assert!(maf_coverage.in_range("chr1", 20));
        assert!(maf_coverage.in_range("chr1", 21));
        assert!(maf_coverage.in_range("chr1", 30));
        assert!(!maf_coverage.in_range("chr1", 31));
        assert!(maf_coverage.in_range("chr1", 35));
        assert!(!maf_coverage.in_range("chr2", 12));
        assert!(maf_coverage.in_range("chrW", 35));
        assert!(!maf_coverage.in_range("chrZ", 36));
    }
}
