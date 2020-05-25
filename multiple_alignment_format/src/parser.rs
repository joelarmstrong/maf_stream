use MAFItem;
use MAFBlock;
use MAFBlockEntry;
use MAFBlockAlignedEntry;
use MAFBlockUnalignedEntry;
use AlignedContext;
use AlignedContextStatus;
use UnalignedContextStatus;
use Strand;
use std::collections::BTreeMap;
use std::io;


pub struct LinesRef<'a, B: 'a> {
    buf: &'a mut B,
}

impl<'a, B: io::BufRead> Iterator for LinesRef<'a, B> {
    type Item = io::Result<String>;

    fn next(&mut self) -> Option<io::Result<String>> {
        let mut buf = String::new();
        match self.buf.read_line(&mut buf) {
            Ok(0) => None,
            Ok(_n) => {
                if buf.ends_with('\n') {
                    buf.pop();
                    if buf.ends_with('\r') {
                        buf.pop();
                    }
                }
                Some(Ok(buf))
            }
            Err(e) => Some(Err(e))
        }
    }
}

#[derive(Debug)]
pub enum MAFParseError {
    IOError(io::Error),
    UnexpectedLine(String),
    BadMetadata,
    BadLineType(String),
    Misc(&'static str),
    EOF,
}

impl From<io::Error> for MAFParseError {
    fn from(err: io::Error) -> Self {
        MAFParseError::IOError(err)
    }
}

/// Get the next MAFItem out of the input.
pub fn next_maf_item<T: io::BufRead + ?Sized>(mut input: &mut T) -> Result<MAFItem, MAFParseError> {
    let mut header: Option<String> = None;
    {
        let lines = LinesRef { buf: &mut input };
        for line_res in lines {
            let line: String = line_res?;
            if line.trim().is_empty() {
                // Blank line
                continue;
            }
            if line.starts_with('#') {
                // MAF comment
                return Ok(MAFItem::Comment(line.chars().skip(1).collect()));
            } else if line.starts_with('a') {
                // Start of a block
                header = Some(line);
                break;
            } else {
                // Shouldn't see this.
                return Err(MAFParseError::UnexpectedLine(line))
            }
        };
    }
    let block = parse_block(header.ok_or(MAFParseError::EOF)?, LinesRef { buf: &mut input })?;
    Ok(MAFItem::Block(block))
}

// Go from "key=value" to "(key, value)".
fn split_metadata_pairs(pair: &str) -> Result<(String, String), MAFParseError> {
    let mut iter = pair.split('=');
    let first = iter.next().ok_or(MAFParseError::BadMetadata)?;
    let second = iter.next().ok_or(MAFParseError::BadMetadata)?;
    Ok((first.to_string(), second.to_string()))
}

// Parse block metadata (the header looks like "a key1=value1 key2=value2").
fn metadata_from_header(header: &str) -> Result<BTreeMap<String, String>, MAFParseError> {
    header.split_whitespace().skip(1).map(split_metadata_pairs).collect()
}

fn parse_strand(strand: &str) -> Result<Strand, MAFParseError> {
    match strand {
        "+" => Ok(Strand::Positive),
        "-" => Ok(Strand::Negative),
        _ => Err(MAFParseError::Misc("Strand not valid")),
    }
}

fn update_from_s_line(fields: &mut Vec<&str>, block_entries: &mut Vec<MAFBlockEntry>) -> Result<(), MAFParseError> {
    let alignment = fields.pop()
        .ok_or(MAFParseError::Misc("s line incomplete"))?;
    let sequence_size = fields.pop()
        .ok_or(MAFParseError::Misc("s line incomplete"))
        .and_then(|s| s.parse::<u64>().map_err(|_| MAFParseError::Misc("invalid sequence size")))?;
    let strand = fields.pop()
        .ok_or(MAFParseError::Misc("s line incomplete"))
        .and_then(parse_strand)?;
    let aligned_length = fields.pop()
        .ok_or(MAFParseError::Misc("s line incomplete"))
        .and_then(|s| s.parse::<u64>().map_err(|_| MAFParseError::Misc("invalid aligned length")))?;
    let start = fields.pop()
        .ok_or(MAFParseError::Misc("s line incomplete"))
        .and_then(|s| s.parse::<u64>().map_err(|_| MAFParseError::Misc("invalid start")))?;
    let seq = fields.pop()
        .ok_or(MAFParseError::Misc("s line incomplete"))?;
    block_entries.push(MAFBlockEntry::AlignedEntry(MAFBlockAlignedEntry {
        alignment: alignment.as_bytes().to_vec(),
        seq: seq.to_string(),
        start,
        aligned_length,
        sequence_size,
        strand,
        context: None,
        qualities: None,
    }));
    Ok(())
}

fn parse_aligned_context_status(status: &str) -> Result<AlignedContextStatus, MAFParseError> {
    use AlignedContextStatus::*;
    match status {
        "C" => Ok(Contiguous),
        "I" => Ok(Insertion),
        "N" => Ok(FirstInSequence),
        "n" => Ok(FirstInSequenceBridged),
        "M" => Ok(MissingData),
        "T" => Ok(AlreadyUsed),
        _   => Err(MAFParseError::Misc("invalid aligned context status"))
    }
}

fn update_from_i_line(fields: &mut Vec<&str>, block_entries: &mut Vec<MAFBlockEntry>) -> Result<(), MAFParseError> {
    let right_count = fields.pop()
        .ok_or(MAFParseError::Misc("i line incomplete"))
        .and_then(|s| s.parse::<u64>().map_err(|_| MAFParseError::Misc("invalid right count")))?;
    let right_status = fields.pop()
        .ok_or(MAFParseError::Misc("i line incomplete"))
        .and_then(parse_aligned_context_status)?;
    let left_count = fields.pop()
        .ok_or(MAFParseError::Misc("i line incomplete"))
        .and_then(|s| s.parse::<u64>().map_err(|_| MAFParseError::Misc("invalid left count")))?;
    let left_status = fields.pop()
        .ok_or(MAFParseError::Misc("i line incomplete"))
        .and_then(parse_aligned_context_status)?;
    let seq = fields.pop()
        .ok_or(MAFParseError::Misc("i line incomplete"))?;

    let context = AlignedContext {
        left_status,
        left_count,
        right_status,
        right_count,
    };

    let last_entry = block_entries.pop()
        .ok_or_else(|| MAFParseError::UnexpectedLine("i line cannot be first in block".to_owned()))?;
    match last_entry {
        MAFBlockEntry::AlignedEntry(mut e) => {
            if e.seq != seq {
                return Err(MAFParseError::UnexpectedLine("i line must follow a corresponding s line".to_owned()))
            }
            e.context = Some(context);
            block_entries.push(MAFBlockEntry::AlignedEntry(e));
            Ok(())
        },
        MAFBlockEntry::UnalignedEntry(_) => Err(MAFParseError::UnexpectedLine("i line must follow a corresponding s line".to_owned())),
    }
}

fn update_from_e_line(fields: &mut Vec<&str>, block_entries: &mut Vec<MAFBlockEntry>) -> Result<(), MAFParseError> {
    let status_char = fields.pop()
        .ok_or(MAFParseError::Misc("e line incomplete"))?;
    let sequence_size = fields.pop()
        .ok_or(MAFParseError::Misc("e line incomplete"))
        .and_then(|s| s.parse::<u64>().map_err(|_| MAFParseError::Misc("invalid sequence size")))?;
    let strand = fields.pop()
        .ok_or(MAFParseError::Misc("e line incomplete"))
        .and_then(parse_strand)?;
    let unaligned_length = fields.pop()
        .ok_or(MAFParseError::Misc("e line incomplete"))
        .and_then(|s| s.parse::<u64>().map_err(|_| MAFParseError::Misc("invalid unaligned length")))?;
    let start = fields.pop()
        .ok_or(MAFParseError::Misc("e line incomplete"))
        .and_then(|s| s.parse::<u64>().map_err(|_| MAFParseError::Misc("invalid start")))?;
    let seq = fields.pop()
        .ok_or(MAFParseError::Misc("e line incomplete"))?;
    let status = match status_char {
        "C" => UnalignedContextStatus::Deletion,
        "I" => UnalignedContextStatus::Insertion,
        "M" => UnalignedContextStatus::MissingData,
        "n" => UnalignedContextStatus::NewSequence,
        "T" => UnalignedContextStatus::AlreadyUsed,
        _   => return Err(MAFParseError::Misc("invalid unaligned context status character")),
    };
    block_entries.push(MAFBlockEntry::UnalignedEntry(MAFBlockUnalignedEntry {
        status,
        seq: seq.to_string(),
        start,
        sequence_size,
        size: unaligned_length,
        strand,
    }));
    Ok(())
}


pub fn parse_block(header: String, iter: impl Iterator<Item = Result<String, io::Error>>) -> Result<MAFBlock, MAFParseError> {
    let mut block_entries: Vec<MAFBlockEntry> = vec![];
    let block_metadata = metadata_from_header(&header)?;
 
    for line_res in iter {
        let line: String = line_res?;
        if line.is_empty() {
            // Blank lines terminate the "paragraph".
            break;
        }
        let mut fields: Vec<_> = line.split_whitespace().collect();
        match fields[0] {
            "s" => update_from_s_line(&mut fields, &mut block_entries)?,
            "i" => update_from_i_line(&mut fields, &mut block_entries)?,
            "e" => update_from_e_line(&mut fields, &mut block_entries)?,
//            "q" => update_from_q_line(&mut fields, &mut block_entries)?,
            _ => return Err(MAFParseError::BadLineType(fields[0].to_string())),
        };
    }
    Ok(MAFBlock {
        metadata: block_metadata,
        entries: block_entries,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufRead, BufReader};

    #[test]
    fn metadata_from_header_filled() {
        let header = "a score=23262.0 pass=2";
        match metadata_from_header(header) {
            Err(e) => assert!(false, "got error {:?}", e),
            Ok(val) => assert_eq!(val, btreemap!{"score".to_owned() => "23262.0".to_owned(),
                                                 "pass".to_owned() => "2".to_owned()}),
        }
    }

    #[test]
    fn metadata_from_header_blank() {
        let header = "a";
        match metadata_from_header(header) {
            Err(e) => assert!(false, "got error {:?}", e),
            Ok(val) => assert_eq!(val, btreemap!{}),
        }
    }

    #[test]
    fn parse_block_only_s_lines() {
        let block_str = "a meta1=val1 meta2=val2
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
s baboon         249182 12 -   4622798 gcagctgaa-aca
s mm4.chr6     53310102 12 + 151104725 ACAGCTGA-AATA

this line is a canary to ensure it stops after a 'paragraph'";
        let mut lines = BufReader::new(block_str.as_bytes()).lines();
        let header = lines.next().unwrap().unwrap();
        match parse_block(header, lines) {
            Err(e) => assert!(false, "got error {:?}", e),
            Ok(val) => assert_eq!(val, MAFBlock {
                metadata: btreemap!{"meta1".to_owned() => "val1".to_owned(),
                                    "meta2".to_owned() => "val2".to_owned()},
                entries: vec![
                    MAFBlockEntry::AlignedEntry(MAFBlockAlignedEntry {
                        seq: "hg16.chr7".to_owned(),
                        start: 27707221,
                        aligned_length: 13,
                        sequence_size: 158545518,
                        strand: Strand::Positive,
                        alignment: "gcagctgaaaaca".as_bytes().to_vec(),
                        context: None,
                        qualities: None,
                    }),
                    MAFBlockEntry::AlignedEntry(MAFBlockAlignedEntry {
                        seq: "baboon".to_owned(),
                        start: 249182,
                        aligned_length: 12,
                        sequence_size: 4622798,
                        strand: Strand::Negative,
                        alignment: "gcagctgaa-aca".as_bytes().to_vec(),
                        context: None,
                        qualities: None,
                    }),
                    MAFBlockEntry::AlignedEntry(MAFBlockAlignedEntry {
                        seq: "mm4.chr6".to_owned(),
                        start: 53310102,
                        aligned_length: 12,
                        sequence_size: 151104725,
                        strand: Strand::Positive,
                        alignment: "ACAGCTGA-AATA".as_bytes().to_vec(),
                        context: None,
                        qualities: None,
                    }),
                ],
            }),
        }
    }

    #[test]
    fn parse_block_i_lines() {
        let block_str = "a
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
i panTro1.chr6 N 0 C 0
s baboon         249182 13 +   4622798 gcagctgaaaaca
i baboon       I 234 n 19";
        let mut lines = BufReader::new(block_str.as_bytes()).lines();
        let header = lines.next().unwrap().unwrap();
        match parse_block(header, lines) {
            Err(e) => assert!(false, "got error {:?}", e),
            Ok(val) => assert_eq!(val, MAFBlock {
                metadata: btreemap!{},
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
                ],
            })
        };
    }

    #[test]
    fn parse_block_e_lines() {
        let block_str = "a
s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
e mm4.chr6     53310102 13 + 151104725 I";
        let mut lines = BufReader::new(block_str.as_bytes()).lines();
        let header = lines.next().unwrap().unwrap();
        match parse_block(header, lines) {
            Err(e) => assert!(false, "got error {:?}", e),
            Ok(val) => assert_eq!(val, MAFBlock {
                metadata: btreemap!{},
                entries: vec![
                    MAFBlockEntry::AlignedEntry(MAFBlockAlignedEntry {
                        seq: "hg16.chr7".to_owned(),
                        start: 27707221,
                        aligned_length: 13,
                        sequence_size: 158545518,
                        strand: Strand::Positive,
                        alignment: "gcagctgaaaaca".as_bytes().to_vec(),
                        context: None,
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
            })
        };
    }


    #[test]
    fn parse_comment() {
        let comment_str = "##maf version=1";
        let mut buf_reader = BufReader::new(comment_str.as_bytes());
        match next_maf_item(&mut buf_reader) {
            Err(e) => assert!(false, "got error {:?}", e),
            Ok(val) => assert_eq!(val, MAFItem::Comment("#maf version=1".to_owned())),
        }
    }

    #[test]
    fn parse_blank_comment() {
        let comment_str = "#";
        let mut buf_reader = BufReader::new(comment_str.as_bytes());
        match next_maf_item(&mut buf_reader) {
            Err(e) => assert!(false, "got error {:?}", e),
            Ok(val) => assert_eq!(val, MAFItem::Comment("".to_owned())),
        }
    }
}
