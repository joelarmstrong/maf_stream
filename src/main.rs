use std::fs::File;
use std::io;
use std::io::{stdout, BufRead, BufReader, Write};
use clap::{App, Arg, SubCommand};

mod dup_blocks;
use dup_blocks::{output_dup_blocks, output_merged_consensus_blocks, ConsensusMode};
mod fasta;
use fasta::maf_to_fasta;

fn main() {
    let matches = App::new("maf_junk")
        .arg(Arg::with_name("input_maf")
             .global(true))
        .arg(Arg::with_name("output")
             .global(true))
        .subcommand(SubCommand::with_name("dup_blocks"))
        .subcommand(SubCommand::with_name("merge_dups")
                    .arg(Arg::with_name("mode")
                         .required(true)
                         .possible_values(&["unanimity",
                                            "consensus",
                                            "mask"])))
        .subcommand(SubCommand::with_name("to_fasta"))
        .get_matches();

    let stdin = io::stdin();
    let mut input = matches.value_of("input_maf").map(|p| Box::new(BufReader::new(File::open(p).expect("Couldn't open input file"))) as Box<BufRead>).unwrap_or(Box::new(stdin.lock()));
    let mut output = matches.value_of("output").map(|p| Box::new(File::create(p).expect("Couldn't create output file")) as Box<Write>).unwrap_or(Box::new(stdout()));

    if let Some(_) = matches.subcommand_matches("dup_blocks") {
        output_dup_blocks(&mut input, &mut output);
    } else if let Some(matches) = matches.subcommand_matches("merge_dups") {
        let mode = match matches.value_of("mode").unwrap() {
            "unanimity" => ConsensusMode::Unanimity,
            "consensus" => ConsensusMode::Consensus,
            "mask" => ConsensusMode::Mask,
            _ => panic!("Unknown consensus mode"),
        };
        output_merged_consensus_blocks(&mut input, &mut output, mode);
    } else if let Some(matches) = matches.subcommand_matches("to_fasta") {
        maf_to_fasta(&mut input, &mut output);
    }
}
