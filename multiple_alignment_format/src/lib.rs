#[cfg(test)]
#[macro_use]
extern crate maplit;

pub mod parser;
pub mod output;
use std::collections::{BTreeMap, HashMap};

/// Structure representing a MAF item (comment or block).
#[derive(Debug, PartialEq, Eq)]
pub enum MAFItem {
    Block(MAFBlock),
    Comment(String),
}

/// A MAF alignment block.
#[derive(Debug, PartialEq, Eq)]
pub struct MAFBlock {
    pub entries: Vec<MAFBlockEntry>,
    pub metadata: BTreeMap<String, String>,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum MAFBlockEntry {
    AlignedEntry(MAFBlockAlignedEntry),
    UnalignedEntry(MAFBlockUnalignedEntry),
}

/// An alignment entry within a MAF block. Corresponds to the "s"
/// line, as well as the "i" and "q" lines if they are present.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct MAFBlockAlignedEntry {
    /// Actual sequence of bases/amino acids, including gaps.
    pub alignment: Vec<u8>,
    /// The sequence name.
    pub seq: String,
    /// Start of the aligned region within this sequence.
    pub start: u64,
    /// Length of the aligned region (not including gaps).
    pub aligned_length: u64,
    /// The total length of this sequence (including regions outside
    /// this alignment).
    pub sequence_size: u64,
    /// Which strand the aligned sequence is on.
    pub strand: Strand,
    /// Context about what's happening before or after the alignment
    /// within this sequence.
    pub context: Option<AlignedContext>,
    /// Optional scores indicating the alignment quality for each
    /// base, ranging from 0-100. No one uses this as far as I can
    /// tell.
    pub qualities: Option<Vec<u8>>,
}

/// Indicates one of the two strands.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum Strand {
    Positive,
    Negative,
}

/// Corresponds to the "i" line.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct AlignedContext {
    pub left_status: AlignedContextStatus,
    pub left_count: u64,
    pub right_status: AlignedContextStatus,
    pub right_count: u64,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum AlignedContextStatus {
    /// Corresponds to "C" in the original MAF line. "The sequence
    /// before or after is contiguous with this block."
    Contiguous,
    /// Corresponds to "I" in the original MAF line. "There are bases
    /// between the bases in this block and the one before or after
    /// it."
    Insertion,
    /// Corresponds to "N" in the original MAF line. "This is the
    /// first sequence from this src chrom or scaffold."
    FirstInSequence,
    /// Corresponds to "n" in the original MAF line. "This is the
    /// first sequence from this src chrom or scaffold but it is
    /// bridged by another alignment from a different chrom or
    /// scaffold."
    FirstInSequenceBridged,
    /// Corresponds to "M" in the original MAF line. "There is missing
    /// data before or after this block (Ns in the sequence)."
    MissingData,
    /// Corresponds to "T" in the original MAF line. "The sequence in
    /// this block has been used before in a previous block (likely a
    /// tandem duplication)."
    AlreadyUsed,
}

/// Indicates that the entry is unaligned, but there is a chain
/// "bridging" two alignment blocks on either side. Corresponds to the
/// "e" line.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct MAFBlockUnalignedEntry {
    pub seq: String,
    /// Start of the unaligned region.
    pub start: u64,
    /// Size of the unaligned region.
    pub size: u64,
    /// Strand aligned to by the chain bridging this unaligned region.
    pub strand: Strand,
    /// Size of the entire sequence.
    pub sequence_size: u64,
    /// The relationship between this unaligned region and the
    /// bridging regions.
    pub status: UnalignedContextStatus,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum UnalignedContextStatus {
    /// "C" in the original MAF line. "The sequence before and after
    /// is contiguous implying that this region was either deleted in
    /// the source or inserted in the reference sequence. The browser
    /// draws a single line or a "-" in base mode in these blocks."
    Deletion,
    /// "I" in the original MAF line. "There are non-aligning bases in
    /// the source species between chained alignment blocks before and
    /// after this block. The browser shows a double line or "=" in
    /// base mode."
    Insertion,
    /// "M" in the original MAF line. "There are non-aligning bases in
    /// the source and more than 90% of them are Ns in the source. The
    /// browser shows a pale yellow bar."
    MissingData,
    /// "n" in the original MAF line. "There are non-aligning bases in
    /// the source and the next aligning block starts in a new
    /// chromosome or scaffold that is bridged by a chain between
    /// still other blocks. The browser shows either a single line or
    /// a double line based on how many bases are in the gap between
    /// the bridging alignments."
    NewSequence,
    /// "T" in the original MAF line. This isn't supposed to happen
    /// according to the docs, but MULTIZ seems to output it. Whatevs.
    AlreadyUsed,
}

impl MAFBlock {
    pub fn aligned_entries(&self) -> impl Iterator<Item=&MAFBlockAlignedEntry> {
        self.entries.iter()
            .filter_map(|e| match e { MAFBlockEntry::AlignedEntry(a) => Some(a), _ => None })
    }

    pub fn entries_as_hash(&self) -> HashMap<&str, Vec<&MAFBlockAlignedEntry>> {
        self.aligned_entries()
            .map(|a| (a.seq.split('.').next().unwrap(), a))
            .fold(HashMap::new(), |mut acc: HashMap<&str, Vec<&MAFBlockAlignedEntry>>, (species, a)| { acc.entry(species).or_insert_with(Vec::new).push(a); acc })
    }
}
