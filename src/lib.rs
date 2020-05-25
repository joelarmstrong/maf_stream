use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::io::BufRead;

#[derive(PartialEq, Eq, Debug)]
pub struct Range {
    pub seq: String,
    pub start: u64,
    pub end: u64,
}

impl Ord for Range {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.seq, &self.start, &self.end).cmp(&(&other.seq, &other.start, &other.end))
    }
}

impl PartialOrd for Range {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Range {
    pub fn overlaps(&self, chrom: &str, position: u64) -> bool {
        self.seq == chrom && self.start <= position && self.end > position
    }

    pub fn precedes(&self, chrom: &str, position: u64) -> bool {
        self.seq.as_str() < chrom || self.end <= position
    }

    pub fn succeeds(&self, chrom: &str, position: u64) -> bool {
        self.seq.as_str() > chrom || self.start > position
    }
}

pub fn parse_bed(bed: impl BufRead) -> BTreeSet<Range> {
    bed.lines()
        .filter_map(|line_res| {
            let line = line_res.expect("Can't read line");
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() > 9 {
                panic!("BED12 input not supported");
            } else if fields.len() > 0 {
                let seq = fields[0].to_string();
                let start: u64 = fields[1].parse().expect("Can't parse start position");
                let end: u64 = fields[2].parse().expect("Can't parse start position");
                Some(Range { seq, start, end })
            } else {
                // Blank line
                None
            }
        })
        .collect()
}

pub fn range_contains_pos(set: &BTreeSet<Range>, chrom: &str, position: u64) -> bool {
    let pos = Range {
        seq: chrom.to_string(),
        start: position + 1,
        end: position + 1,
    };
    match set.range(..=pos).next_back() {
        None => false,
        Some(range) => range.overlaps(chrom, position),
    }
}

/// Gives (potentially) overlapping ranges
pub fn overlapping_ranges<'a>(
    set: &'a BTreeSet<Range>,
    range: &Range,
) -> impl Iterator<Item = &'a Range> {
    let end = Range {
        seq: range.seq.clone(),
        start: range.end,
        end: range.end,
    };
    set.range(..=range)
        .rev()
        .take(1)
        .chain(set.range(range..=&end))
}

/// Get "chr.name" from "genome.chr.name".
pub fn chrom_part(seq: &str) -> String {
    seq.split(".").skip(1).join(".")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_range_contains_pos() {
        let regions: BTreeSet<_> = vec![
            Range {
                seq: "chr1".to_string(),
                start: 4432333,
                end: 4432334,
            },
            Range {
                seq: "chr1".to_string(),
                start: 4432333,
                end: 4432335,
            },
            Range {
                seq: "chr1".to_string(),
                start: 4432336,
                end: 4432338,
            },
            Range {
                seq: "chr2".to_string(),
                start: 4,
                end: 5,
            },
        ]
        .into_iter()
        .collect();
        assert!(range_contains_pos(&regions, "chr1", 4432333));
        assert!(range_contains_pos(&regions, "chr1", 4432334));
        assert!(!range_contains_pos(&regions, "chr1", 4432335));
        assert!(!range_contains_pos(&regions, "chr1", 4));
        assert!(range_contains_pos(&regions, "chr2", 4));
    }

    #[test]
    fn test_overlapping_ranges() {
        let regions: BTreeSet<_> = vec![
            Range {
                seq: "chr1".to_string(),
                start: 4432333,
                end: 4432334,
            },
            Range {
                seq: "chr1".to_string(),
                start: 4432333,
                end: 4432335,
            },
            Range {
                seq: "chr1".to_string(),
                start: 4432336,
                end: 4432338,
            },
            Range {
                seq: "chr2".to_string(),
                start: 4,
                end: 5,
            },
        ]
        .into_iter()
        .collect();

        assert_eq!(
            overlapping_ranges(
                &regions,
                &Range {
                    seq: "chr1".to_string(),
                    start: 4430000,
                    end: 8000000
                }
            )
            .collect::<Vec<_>>(),
            vec![
                &Range {
                    seq: "chr1".to_string(),
                    start: 4432333,
                    end: 4432334,
                },
                &Range {
                    seq: "chr1".to_string(),
                    start: 4432333,
                    end: 4432335,
                },
                &Range {
                    seq: "chr1".to_string(),
                    start: 4432336,
                    end: 4432338,
                }
            ]
        );
        assert_eq!(
            overlapping_ranges(
                &regions,
                &Range {
                    seq: "chr1".to_string(),
                    start: 4430000,
                    end: 4432335
                }
            )
            .collect::<Vec<_>>(),
            vec![
                &Range {
                    seq: "chr1".to_string(),
                    start: 4432333,
                    end: 4432334,
                },
                &Range {
                    seq: "chr1".to_string(),
                    start: 4432333,
                    end: 4432335,
                }
            ]
        );
        assert_eq!(
            overlapping_ranges(
                &regions,
                &Range {
                    seq: "chr1".to_string(),
                    start: 4432335,
                    end: 4432337
                }
            )
            .collect::<Vec<_>>(),
            vec![
                &Range {
                    seq: "chr1".to_string(),
                    start: 4432333,
                    end: 4432335,
                },
                &Range {
                    seq: "chr1".to_string(),
                    start: 4432336,
                    end: 4432338,
                }
            ]
        );
    }
}
