use std::fs::File;
use std::io;
use std::io::{stdout, BufRead, BufReader, Write};
use clap::{App, Arg, SubCommand, value_t};

mod dup_blocks;
use dup_blocks::{output_dup_blocks, output_merged_consensus_blocks, ConsensusMode};
mod split;
use split::split_maf;
mod coverage;
use coverage::coverage;

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
        .subcommand(SubCommand::with_name("split")
                    .arg(Arg::with_name("output_dir")
                         .required(true))
                    .arg(Arg::with_name("max_length")))
        .subcommand(SubCommand::with_name("coverage")
                    .arg(Arg::with_name("ref_genome")
                         .required(true))
                    .arg(Arg::with_name("bed")
                         .long("bed")
                         .takes_value(true)))
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
    } else if let Some(_matches) = matches.subcommand_matches("to_fasta") {
        unimplemented!();
        //maf_to_fasta(&mut input, &mut output);
    } else if let Some(matches) = matches.subcommand_matches("split") {
        let max_length = value_t!(matches, "max_length", u64).unwrap_or(100000);
        split_maf(&mut input, max_length, matches.value_of("output_dir").unwrap());
    } else if let Some(matches) = matches.subcommand_matches("coverage") {
        let bed_file = matches.value_of("bed").map(|path| BufReader::new(File::open(path).expect("Couldn't open bed file")));
        let ref_genome = matches.value_of("ref_genome").unwrap();
        coverage(&mut input, &mut output, ref_genome, bed_file);
    }
}
