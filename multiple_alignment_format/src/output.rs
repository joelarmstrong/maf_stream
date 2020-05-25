use MAFBlock;
use MAFBlockEntry;
use Strand;
use AlignedContextStatus;
use UnalignedContextStatus;
use std::str;
use std::fmt;

fn aligned_context_status_char(status: &AlignedContextStatus) -> &'static str {
    use AlignedContextStatus::*;
    match status {
        Contiguous => "C",
        Insertion => "I",
        FirstInSequence => "N",
        FirstInSequenceBridged => "n",
        MissingData => "M",
        AlreadyUsed => "T",
    }
}

fn unaligned_context_status_char(status: &UnalignedContextStatus) -> &'static str {
    use UnalignedContextStatus::*;
    match status {
        Deletion => "C",
        Insertion => "I",
        MissingData => "M",
        NewSequence => "n",
        AlreadyUsed => "T",
    }
}

impl fmt::Display for MAFBlock {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "a")?;
        for (key, value) in self.metadata.iter() {
            write!(f, " {}={}", key, value)?;
        }
        writeln!(f)?;
        for entry in self.entries.iter() {
            match entry {
                MAFBlockEntry::AlignedEntry(e) => {
                    writeln!(f, "s {} {} {} {} {} {}",
                           e.seq, e.start, e.aligned_length,
                           if e.strand == Strand::Positive { "+" } else { "-" },
                           e.sequence_size,
                           str::from_utf8(&e.alignment).expect("alignment not utf8 compatible"))?;
                    if let Some(ref context) = e.context {
                        writeln!(f, "i {} {} {} {} {}",
                               e.seq,
                               aligned_context_status_char(&context.left_status),
                               context.left_count,
                               aligned_context_status_char(&context.right_status),
                               context.right_count)?;
                    }
                },
                MAFBlockEntry::UnalignedEntry(e) => {
                    writeln!(f, "e {} {} {} {} {} {}",
                    e.seq, e.start, e.size,
                    if e.strand == Strand::Positive { "+" } else { "-" },
                    e.sequence_size,
                    unaligned_context_status_char(&e.status))?;
                },
            }
        }
        writeln!(f)
    }
}

#[cfg(test)]
mod tests {
    use MAFBlock;
    use MAFBlockEntry;
    use MAFBlockAlignedEntry;
    use MAFBlockUnalignedEntry;
    use Strand;
    use AlignedContext;
    use AlignedContextStatus;
    use UnalignedContextStatus;
    #[test]
    fn display_block() {
        let block = MAFBlock {
            metadata: btreemap!{"meta1".to_owned() => "val1".to_owned(),
                                "meta2".to_owned() => "val2".to_owned()},
            entries: vec![
                MAFBlockEntry::AlignedEntry(MAFBlockAlignedEntry {
                    seq: "panTro1.chr6".to_owned(),
                    start: 28869787,
                    aligned_length: 13,
                    sequence_size: 161576975,
                    strand: Strand::Positive,
                    alignment: "gcagctgaaaaca".as_bytes().to_vec(),
                    context: Some(AlignedContext {
                        left_status: AlignedContextStatus::FirstInSequence,
                        left_count: 0,
                        right_status: AlignedContextStatus::Contiguous,
                        right_count: 0,
                    }),
                    qualities: None,
                }),
                MAFBlockEntry::AlignedEntry(MAFBlockAlignedEntry {
                    seq: "baboon".to_owned(),
                    start: 249182,
                    aligned_length: 13,
                    sequence_size: 4622798,
                    strand: Strand::Positive,
                    alignment: "gcagctgaaaaca".as_bytes().to_vec(),
                    context: Some(AlignedContext {
                        left_status: AlignedContextStatus::Insertion,
                        left_count: 234,
                        right_status: AlignedContextStatus::FirstInSequenceBridged,
                        right_count: 19,
                    }),
                    qualities: None,
                }),
               MAFBlockEntry::UnalignedEntry(MAFBlockUnalignedEntry {
                    seq: "mm4.chr6".to_owned(),
                    start: 53310102,
                    size: 13,
                    sequence_size: 151104725,
                    strand: Strand::Positive,
                    status: UnalignedContextStatus::Insertion,
                }),
            ],
        };
        assert_eq!(block.to_string(), "a meta1=val1 meta2=val2
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
i panTro1.chr6 N 0 C 0
s baboon 249182 13 + 4622798 gcagctgaaaaca
i baboon I 234 n 19
e mm4.chr6 53310102 13 + 151104725 I

");
    }

}
