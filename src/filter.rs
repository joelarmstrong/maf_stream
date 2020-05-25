use crate::lib::{chrom_part, overlapping_ranges, parse_bed, Range};
use multiple_alignment_format::parser::next_maf_item;
use multiple_alignment_format::{MAFBlock, MAFBlockAlignedEntry, MAFBlockEntry, MAFItem, Strand};
use std::collections::BTreeSet;
use std::io::{BufRead, Write};

/// Run of columns.
#[derive(Debug, PartialEq)]
struct Run {
    start: usize,
    length: usize,
}

/// Get the columns within the block which should be kept.
fn get_filtered_columns(ref_entry: &MAFBlockAlignedEntry, ranges: &BTreeSet<Range>) -> Vec<Run> {
    assert!(ref_entry.strand == Strand::Positive);
    let chrom = chrom_part(&ref_entry.seq);
    let mut runs: Vec<Run> = vec![];
    let mut relevant_ranges = overlapping_ranges(
        ranges,
        &Range {
            seq: chrom_part(&ref_entry.seq),
            start: ref_entry.start,
            end: ref_entry.start + ref_entry.aligned_length,
        },
    );
    let mut current_range = relevant_ranges.next();
    let mut current_pos = ref_entry.start;
    let mut was_within_run = false;
    for (i, c) in ref_entry.alignment.iter().enumerate() {
        while current_range.map_or(false, |r| r.precedes(&chrom, current_pos)) {
            current_range = relevant_ranges.next();
        }
        if current_range.is_none()
            || current_range.map_or(false, |r| {
                r.succeeds(&chrom, ref_entry.start + ref_entry.aligned_length)
            })
        {
            break;
        }
        let mut within_run = false;
        if *c != b'-' {
            if current_range.unwrap().overlaps(&chrom, current_pos) {
                if was_within_run {
                    runs.last_mut().unwrap().length += 1;
                } else {
                    runs.push(Run {
                        start: i,
                        length: 1,
                    });
                }
                within_run = true;
            }
            current_pos += 1;
        }
        was_within_run = within_run;
    }
    runs
}

fn filter_entry_columns(entry: &MAFBlockAlignedEntry, run: &Run) -> MAFBlockAlignedEntry {
    let before_range_offset = entry.alignment[..run.start]
        .iter()
        .filter(|c| **c != b'-')
        .count() as u64;
    let inside_range_offset = entry.alignment[run.start..run.start + run.length]
        .iter()
        .take_while(|c| **c == b'-')
        .count() as u64;

    MAFBlockAlignedEntry {
        seq: entry.seq.clone(),
        sequence_size: entry.sequence_size,
        strand: entry.strand,
        start: entry.start + before_range_offset + inside_range_offset,
        alignment: entry.alignment[run.start..run.start + run.length].to_vec(),
        aligned_length: entry.alignment[run.start..run.start + run.length]
            .iter()
            .filter(|c| **c != b'-')
            .count() as u64,
        // TODO. But no one uses/cares about these anyway
        context: None,
        qualities: None,
    }
}

fn filter_block_columns(block: &MAFBlock, run: &Run) -> MAFBlock {
    MAFBlock {
        entries: block
            .aligned_entries()
            .map(|e| MAFBlockEntry::AlignedEntry(filter_entry_columns(e, run)))
            .collect(),
        metadata: block.metadata.clone(),
    }
}

fn filter_block(block: &MAFBlock, ranges: &BTreeSet<Range>) -> Vec<MAFBlock> {
    if block.aligned_entries().next().is_none() {
        return vec![];
    }
    let ref_entry = block.aligned_entries().next().unwrap();
    get_filtered_columns(ref_entry, ranges)
        .iter()
        .map(|run| filter_block_columns(block, run))
        .collect()
}

pub fn filter(input: &mut dyn BufRead, output: &mut dyn Write, bed: impl BufRead) {
    let ranges = parse_bed(bed);

    while let Ok(item) = next_maf_item(input) {
        match item {
            MAFItem::Comment(comment) => {
                write!(output, "#{}\n", comment).ok();
            }
            MAFItem::Block(block) => {
                for filtered_block in filter_block(&block, &ranges) {
                    write!(output, "{}", filtered_block).ok();
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_block() {
        let block = "a
s       Gallus_gallus.chr1 4432333   5       +       157682039  CAGT-A
s       Alca_torda.scaffold4709 42333   6       -       157682  TAGTAA
s       Alca_torda.scaffold4709 41641   3       -       157682  G-AA--
";
        let item = next_maf_item(&mut block.as_bytes()).expect("Couldn't parse MAF block");
        if let MAFItem::Block(block) = item {
            let block = filter_block_columns(
                &block,
                &Run {
                    start: 2,
                    length: 3,
                },
            );
            assert_eq!(
                format!("{}", block),
                "a
s Gallus_gallus.chr1 4432335 2 + 157682039 GT-
s Alca_torda.scaffold4709 42335 3 - 157682 GTA
s Alca_torda.scaffold4709 41642 2 - 157682 AA-

"
            );
        } else {
            assert!(false, "Got unexpected maf item {:?}", item);
        }
    }

    #[test]
    fn test_get_filtered_columns() {
        let block = "a
s       Gallus_gallus.chr1 4432333   5       +       157682039  CAGT-A
s       Alca_torda.scaffold4709 42333   6       -       157682  TAGTAA
s       Alca_torda.scaffold4709 41641   3       -       157682  G-AA--
";
        let regions: BTreeSet<_> = vec![
            Range {
                seq: "chr1".to_string(),
                start: 4432333,
                end: 4432334,
            },
            Range {
                seq: "chr1".to_string(),
                start: 4432336,
                end: 4432338,
            },
        ]
        .into_iter()
        .collect();

        let item = next_maf_item(&mut block.as_bytes()).expect("Couldn't parse MAF block");
        if let MAFItem::Block(block) = item {
            let ref_entry = block.aligned_entries().next().unwrap();
            assert_eq!(
                get_filtered_columns(ref_entry, &regions),
                vec![
                    Run {
                        start: 0,
                        length: 1
                    },
                    Run {
                        start: 3,
                        length: 1
                    },
                    Run {
                        start: 5,
                        length: 1
                    }
                ]
            );
        } else {
            assert!(false, "Got unexpected maf item {:?}", item);
        }
    }
}
