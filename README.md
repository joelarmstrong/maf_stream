## maf_stream
`maf_stream` is a collection of utilities to manipulate multiple alignments in the Multiple Alignment Format (MAF).
### Installation
First, install [Rust](https://rustup.rs/). Then, run:
```
cargo install --path .
```
which will install the program into your Cargo bin directory, or use `cargo build --release` and use the binary in `target/release/maf_stream`.
### Usage
By default, `maf_stream` sets the input MAF to stdin and the output to stdout.
#### Finding duplicated blocks
`maf_stream dup_blocks <input maf> <output maf>`
#### Resolving duplicated entries
`maf_stream merge_dups <merging mode> <input maf> <output maf>`
The resulting blocks always contain at most one entry per species; where there were previously duplicated entries only one entry will remain.

Available merging modes:
- `consensus`: Replace duplicated entries with a single entry (containing a consensus of the dups, with ties broken by consensus with the rest of the column).
- `unanimity`: Replace duplicated entries with a single entry (containing N if there are different bases within the duplicates, containing the unanimous base if the duplicate entries all agree).
- `mask`: Replace duplicated entries with single masked entry (containing all Ns).
#### Splitting a MAF (by reference sequence and maximum length)
`maf_stream <output dir> --max_length <max length per file> <input maf>`
#### Calculating coverage
`maf_stream <reference genome> [--bed BED_FILE] <input maf> <output file>`

If `--bed BED_FILE` is provided, coverage is restricted to be of bases within the regions within the BED file. Note that the BED file should not contain overlaps, i.e. it should be run through `bedtools merge` before being used. BED12 input is also currently disallowed, but will work if split up into BED3.

The output is similar to [mafCoverage](https://github.com/dentearl/mafTools/tree/master/mafCoverage).
