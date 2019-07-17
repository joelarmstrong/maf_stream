use multiple_alignment_format::parser::next_maf_item;
use multiple_alignment_format::{MAFBlock, MAFItem, MAFBlockEntry, MAFBlockAlignedEntry};
use std::collections::HashMap;
use std::io::{BufRead, Write};

fn dup_entries_from_block(block: &MAFBlock) -> HashMap<&str, Vec<&MAFBlockAlignedEntry>> {
    let mut hash = block.entries_as_hash();
    hash.retain(|_, v| v.len() > 1);
    hash
}

fn block_contains_dups(block: &MAFBlock) -> bool {
    !dup_entries_from_block(block).is_empty()
}

/// When merging 2+ alignment entries from the same species within a
/// single block, this describes what the base call will be for the
/// merged entry.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum ConsensusMode {
    /// The base is kept as a nucleotide if all dups have the same
    /// base in this position, and N otherwise.
    Unanimity,
    /// The base is set to the consensus of the dups within the
    /// species, using the most frequent base within the entire column
    /// to break ties within the species. If there is still a tie even
    /// after considering the entire column, the base is set to N.
    Consensus,
    /// All bases which are duplicated are set to N for that species's block entry.
    Mask,
}

fn unanimous_base(base_counts: &BaseCounts) -> u8 {
    if base_counts.a > 0 && base_counts.c == 0 && base_counts.g == 0 && base_counts.t == 0 {
        b'A'
    } else if base_counts.a == 0 && base_counts.c > 0 && base_counts.g == 0 && base_counts.t == 0 {
        b'C'
    } else if base_counts.a == 0 && base_counts.c == 0 && base_counts.g > 0 && base_counts.t == 0 {
        b'G'
    } else if base_counts.a == 0 && base_counts.c == 0 && base_counts.g == 0 && base_counts.t > 0 {
        b'T'
    } else {
        b'N'
    }
}

fn max_among_possibilities(base_counts: &BaseCounts, possibilities: &mut [bool; 4]) {
    let mut max_so_far = 0;
    if possibilities[0] {
        max_so_far = base_counts.a;
    }
    if possibilities[1] && max_so_far < base_counts.c {
        possibilities[0] = false;
        max_so_far = base_counts.c;
    } else if base_counts.c < max_so_far {
        possibilities[1] = false;
    }
    if possibilities[2] && max_so_far < base_counts.g {
        possibilities[0] = false;
        possibilities[1] = false;
        max_so_far = base_counts.g;
    } else if base_counts.g < max_so_far {
        possibilities[2] = false;
    }
    if possibilities[3] && max_so_far < base_counts.t {
        possibilities[0] = false;
        possibilities[1] = false;
        possibilities[2] = false;
    } else if base_counts.t < max_so_far {
        possibilities[3] = false;
    }
}

fn consensus_base(base_counts: &BaseCounts, tie_breaker: &BaseCounts) -> u8 {
    let mut possible_bases = [true, true, true, true];
    max_among_possibilities(base_counts, &mut possible_bases);
    max_among_possibilities(tie_breaker, &mut possible_bases);
    match possible_bases {
        [true, false, false, false] => b'A',
        [false, true, false, false] => b'C',
        [false, false, true, false] => b'G',
        [false, false, false, true] => b'T',
        _ => b'N',
    }
}

fn merge_dup_entries(dup_entries: &HashMap<&str, Vec<&MAFBlockAlignedEntry>>, block_consensus: &[BaseCounts], mode: ConsensusMode) -> Vec<MAFBlockEntry> {
    let mut merged_entries = vec![];
    for (_, alignments) in dup_entries.iter() {
        let mut merged_alignment = alignments[0].clone();
        let base_counts = get_consensus_info(&alignments);
        for i in 0..merged_alignment.alignment.len() {
            merged_alignment.alignment[i] = match mode {
                ConsensusMode::Mask => b'N',
                ConsensusMode::Unanimity => unanimous_base(&base_counts[i]),
                ConsensusMode::Consensus => consensus_base(&base_counts[i], &block_consensus[i]),
            }
        }
        merged_entries.push(MAFBlockEntry::AlignedEntry(merged_alignment));
    }
    merged_entries
}

#[derive(Debug, PartialEq)]
struct BaseCounts {
    a: usize,
    c: usize,
    g: usize,
    t: usize,
}

fn get_consensus_info(entries: &[&MAFBlockAlignedEntry]) -> Vec<BaseCounts> {
    let bases = [b'a', b'c', b'g', b't'];
    let mut counts = vec![];
    let length = entries[0].alignment.len();
    for i in 0..length {
        let mut count_iter = bases.iter()
            .map(|b| entries.iter().filter(|e| e.alignment[i].eq_ignore_ascii_case(b)).count());
        counts.push(BaseCounts {
            a: count_iter.next().unwrap(),
            c: count_iter.next().unwrap(),
            g: count_iter.next().unwrap(),
            t: count_iter.next().unwrap(),
        })
    }
    counts
}

pub fn output_merged_consensus_blocks(input: &mut dyn BufRead, output: &mut dyn Write, mode: ConsensusMode) {
    while let Ok(item) = next_maf_item(input) {
        match item {
            MAFItem::Comment(comment) => { write!(output, "#{}\n", comment).ok(); },
            MAFItem::Block(mut block) => {
                let dup_entries = dup_entries_from_block(&block);
                let aligned_entries: Vec<_> = block.aligned_entries().collect();
                let block_counts = get_consensus_info(&aligned_entries);

                // Clear out duplicated entries within the block.
                let values: Vec<_> = dup_entries.values().flatten().collect();
                let mut new_block_entries = block.entries.clone();
                new_block_entries.retain(|e| match e { MAFBlockEntry::AlignedEntry(a) => !values.contains(&&a), _ => true });
                let dup_entries = merge_dup_entries(&dup_entries, &block_counts, mode);
                block.entries = new_block_entries;
                block.entries.extend(dup_entries);
                write!(output, "{}\n", block).ok();
            },
        }
    }
}

pub fn output_dup_blocks(input: &mut dyn BufRead, output: &mut dyn Write) {
    while let Ok(item) = next_maf_item(input) {
        match item {
            MAFItem::Comment(comment) => { write!(output, "#{}\n", comment).ok(); },
            MAFItem::Block(block) => {
                if block_contains_dups(&block) {
                    write!(output, "{}", block).ok();
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unanimous_base() {
        let counts = BaseCounts {
            a: 0,
            c: 10,
            g: 0,
            t: 0,
        };
        assert_eq!(unanimous_base(&counts), b'C');

        let counts = BaseCounts {
            a: 0,
            c: 0,
            g: 0,
            t: 10,
        };
        assert_eq!(unanimous_base(&counts), b'T');

        // Should give N if any disagreement
        let counts = BaseCounts {
            a: 9,
            c: 0,
            g: 0,
            t: 1,
        };
        assert_eq!(unanimous_base(&counts), b'N');
    }

    #[test]
    fn test_consensus_base() {
        // If a majority points to a base, choose that
        let counts = BaseCounts {
            a: 6,
            c: 5,
            g: 5,
            t: 5,
        };
        let tiebreakers = BaseCounts {
            a: 0,
            c: 0,
            g: 0,
            t: 0,
        };
        assert_eq!(consensus_base(&counts, &tiebreakers), b'A');

        // If a tie exists, go to the tiebreaking counts to resolve it
        let counts = BaseCounts {
            a: 4,
            c: 4,
            g: 5,
            t: 5,
        };
        let tiebreakers = BaseCounts {
            a: 1,
            c: 1,
            g: 2,
            t: 1,
        };
        assert_eq!(consensus_base(&counts, &tiebreakers), b'G');

        // If a tie still exists even after tiebreaking, return an N
        let counts = BaseCounts {
            a: 4,
            c: 5,
            g: 4,
            t: 5,
        };
        let tiebreakers = BaseCounts {
            a: 1,
            c: 2,
            g: 1,
            t: 2,
        };
        assert_eq!(consensus_base(&counts, &tiebreakers), b'N');
    }

    #[test]
    fn test_get_consensus_info() {
        let block = "a
s       Gallus_gallus.chr1 4432333   6       +       157682039  CAG
s       Alca_torda.scaffold4709 42333   6       -       157682  TAG
s       Alca_torda.scaffold4709 41641   6       -       157682  G-A
";
        let item = next_maf_item(&mut block.as_bytes()).expect("Couldn't parse MAF block");
        if let MAFItem::Block(block) = item {
            let alignments: Vec<_> = block.aligned_entries().collect();
            let counts = get_consensus_info(&alignments);
            let expected_counts = vec![
                BaseCounts { a: 0, c: 1, t: 1, g: 1, },
                BaseCounts { a: 2, c: 0, t: 0, g: 0, },
                BaseCounts { a: 1, c: 0, t: 0, g: 2, },
            ];
            assert_eq!(counts, expected_counts);
        } else {
            assert!(false, "Got unexpected maf item {:?}", item);
        }
    }

    #[test]
    fn test_block_contains_dups() {
        let dup_block = "a
s       Gallus_gallus.chr1 4432333   6       +       157682039  CAACAG
s       Alca_torda.scaffold4709 42333   6       -       157682  CAACAG
s       Alca_torda.scaffold4709 41641   6       -       157682  CAACAG
";
        let item = next_maf_item(&mut dup_block.as_bytes()).expect("Couldn't parse MAF block");
        if let MAFItem::Block(block) = item {
            assert!(block_contains_dups(&block));
        } else {
            assert!(false, "Got unexpected maf item {:?}", item);
        }

        let non_dup_block = "a
s       Erythrocercus_mccallii.scaffold_2093    58535   1       +       127396  T
s       Eubucco_bourcierii.scaffold13745        58548   1       +       73788   C
s       Eudromia_elegans.scaffold_5     12617300        1       +       13876364        A
s       Eulacestoma_nigropectus.scaffold148     3647521 1       +       4202789 T
s       Eurypyga_helias.scaffold13804   27799   1       +       30980   G
s       Eurystomus_gularis.scaffold6487 121546  1       +       203918  T
s       Formicarius_rufipectus.scaffold473      3224110 1       -       3420713 T
s       Fregata_magnificens.C5769372__2.0       142     1       +       150     T
s       Fregetta_grallaria.scaffold_297 174414  1       -       673556  T
s       Fulmarus_glacialis.scaffold7044 80945   1       +       82154   T
s       Furnarius_figulus.scaffold_634  115343  1       +       392412  T
s       Galbula_dea.scaffold1422        3938    1       -       1348798 C
s       Gavia_stellata.scaffold9486     35556   1       +       49599   T
s       Geococcyx_californianus.scaffold6221    68277   1       +       96248   T
s       Geospiza_fortis.scaffold54      15705654        1       -       19033121        T
s       Glareola_pratincola.scaffold_8  396272  1       -       2357087 C
s       Glaucidium_brasilianum.scaffold_161     1648450 1       -       1875072 T
";
        let item = next_maf_item(&mut non_dup_block.as_bytes()).expect("Couldn't parse MAF block");
        if let MAFItem::Block(block) = item {
            assert!(!block_contains_dups(&block));
        } else {
            assert!(false, "Got unexpected maf item {:?}", item);
        }
    }
}
