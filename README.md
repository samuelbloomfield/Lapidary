# Lapidary
Identifying amino acid sequences using sequenced reads

## Introduction
Lapidary uses Diamond to identify reads that align to an amino acid sequence database and calculates the amino acid coverage, identity and mean read depth using the  translated sections of the read alignment.

## Quick start
Bring up the full list of options:
```
%perl Lapidary.pl -version

Lapidary: a software for identifying amino acid sequences using sequenced reads

	Options:
	read_1			Location of first read file (required)
	read_2			Location of second read file if read files are paired
	db			Full location to fasta file containing amino acid sequences (required)
	threads			Number of threads to use for Diamond (default: 1)
	identity		Diamond identity percentage cut-off to use (default: 80)
	coverage		Diamond coverage percentage cut-off to use (default: 20)
	read_type		Types of reads used (required): single or paired
	sequence_identification	Method for calling most likely sequence: identity (default) or consensus
	help			Display help screen
	version			Return version of Lapidary
```

Run Lapidary on single read files:
```
perl Lapidary.pl -read_1 Examples/Reads/Single_reads.fq.gz --read_type single -db Examples/Amino_acid_database.fasta
```

Run Lapidary on paired read files:
```
perl Lapidary.pl -read_1 Examples/Reads/Paired_reads_1.fq.gz -read_2 Examples/Reads/Paired_reads_2.fq.gz -read_type paired -db Examples/Amino_acid_database.fasta
```
## Output
Lapidary produces a tab-delimited file for each read file or read pair with the following columns:

| Column | Description | Example |
| --- | --- | --- |
| Protein | Reference amino acid sequence | BAC0273 |
| Coverage | Proportion of reference amino acid sequence that reads alined to | 0.7563248764 |
| Identity | Identity proportion between reads and reference amino acid sequence and reads | 0.9515154354 |
| Mean_read_depth | Mean read depth for the reference amino acid sequence | 1.2111562156 |
| Alignment_start | Position of first amino acid that reads aligned to | 23 |
| Alignment_end | Position of last amino acid that reads aligned to | 324 |
| Most_likely_sequence | Most likely sequence predicted | VGVAVYELDLFGRLRNL |

## Requirements
 - perl >=5.32
 - diamond >=2.1.7
 
## Installation
Lapidary can be installed using the following cpanm command:
 - cpanm App::lapidary
