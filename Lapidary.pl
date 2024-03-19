###
#	Lapidary.pl
#	Script to determine the read coverage and identity to a protein database using diamond
#	Last modified: 19/03/2024
#	Author: Samuel Bloomfield
###

use warnings;
use strict;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use FindBin qw($RealBin);
use File::Basename;
use File::Spec::Functions qw(catfile);
use LWP::Simple;
use Archive::Extract;


# here check if no arguments were given, then print help message
if (@ARGV == 0) {
	# Print the help message
	print "This script requires some arguments to run.\n";
	print "use --help to see the options\n";
	exit;
}

my $script_location = $RealBin;
my $path_variable = $ENV{'PATH'};
my $read_1;
my $read_2="";
my $db;
my $threads = 1; #Default to 1 thread
my $identity = 70; #Default to 70% identity
my $coverage = 50; #Default to 50% coverage
my $read_type;
my $help;
my $version;
my $sequence_identification = "identity";

# Check if the script location is already in the PATH
if ($path_variable !~ m/$script_location/) {
	# Add the script location to the PATH
	$ENV{'PATH'} = "$path_variable:$script_location";
}

GetOptions (	'read_1:s'		=> \$read_1,
		'read_2:s'		=> \$read_2,
		'db:s'			=> \$db,
		'threads:i'		=> \$threads,
		'identity:i'	=> \$identity,
		'coverage:i'	=> \$coverage,
		'read_type:s'	=> \$read_type,
		'sequence_identification:s'	=> \$sequence_identification,
		'help'			=> \$help,
		'version'		=> \$version);

#Print out help message
if(defined($help)){
	die "\n\nLapidary: a software for identifying amino acid sequences using sequenced reads\n\n
	Options:\n
	read_1\tLocation of first read file (required)\n
	read_2\tLocation of second read file if read files are paired\n
	db\tFull location to fasta file containing amino acid sequences (required)\n
	threads\tNumber of threads to use for Diamond (default: 1)\n
	identity\tDiamond identity percentage cut-off to use (default: 70)\n
	coverage\tDiamond coverage percentage cut-off to use (default: 50)\n
	read_type\tTypes of reads used (required): single or paired\n
	sequence_identification\tMethod for calling most likely sequence: identity (default) or consensus\n
	help\tDisplay help screen\n
	version\tReturn version of Lapidary\n\n";
} 

#Print out version
if(defined($version)){
	die "\n\nLapidary version 0.5.0\n\n";
} 


# Check if we have diamond
# Define the URL to download diamond
my $diamond_url = 'http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz';

# Check if diamond is in the PATH
unless (system("which diamond >/dev/null 2>&1") == 0) {
	print "diamond not found in PATH. Downloading and installing...\n";
	
	# Get the filename from the URL
	my $filename = basename($diamond_url);
	
	# Set the destination path for download and extraction
	my $download_path = catfile($RealBin, $filename);
	my $extraction_path = $RealBin;
	
	# Download the diamond tar.gz file
	my $status = getstore($diamond_url, $download_path);
	
	# Check if download was successful
	unless ($status == 200) {
		die "Failed to download diamond: $status";
	}
	
	# Extract the tar.gz file
	my $extractor = Archive::Extract->new(archive => $download_path);
	my $result = $extractor->extract(to => $extraction_path);
	
	# Check if extraction was successful
	unless ($result) {
		die "Failed to extract diamond: " . $extractor->error;
	}
	
	# Print success message
	print "diamond has been downloaded and installed in $extraction_path\n";
}


#Create diamond database of protein sequences
my $db_name;
if ($db=~m/^.+\/(.*?)\.f.*?$/) {
	$db_name = $1;
	
	system "diamond makedb --in $db --db $db_name"

} else {
	die "File $db not in fasta format $!\n";
}

#Open fasta file and extract amino acid sequences
my @protein_names=();
my @protein_sequences=();

my @sequence_partial=();

open DB, $db or die "File $db NOT FOUND - $!\n";
foreach $_ (<DB>) {
	if ($_=~m/^>(.*?)$/) {
	push @protein_names, $1;
	if((scalar(@sequence_partial))>0){
		push @protein_sequences, join("", @sequence_partial);
		@sequence_partial = ();
	}
	} elsif ($_=~m/^(.*?)$/) { 
	push @sequence_partial, $1;
	}
}

push @protein_sequences, join("", @sequence_partial);
close DB;

my $concatenated;

my $diamond_output;
my $protein_cut_off;

my $read_protein;
my $read_name;
my $read_translate;
my $read_position;
my $read_start;
my $read_end;
my $read_length;
my $protein_length;

my $start_temp;
my $end_temp;

my $read_alignment_QC;

my @sample_proteins;
my @sample_reads;
my @sample_positions;
my @sample_amino_acids;

my @protein_reads;
my @protein_positions;
my @protein_amino_acids;

my @unique_proteins;
my @unique_sequences;

my @amino_acids_split;

my @protein_split;

my $depth;
my $match;
my $cover;

my @indexes;
my @index_sorted;

my $position_depth;

my $identity_proportion;
my $coverage_proportion;
my $mean_read_depth;

my @protein_indexes;
my @protein_index_sorted;

my @position_depth_array;
my @position_cover_array;
my @position_match_array;
my @position_start_array;
my @position_end_array;
my @index_sorted_size;

my $sample_position_current;

my $alignment_start;
my $alignment_end;

my @consensus;
my @uniq_consensus;
my @uniq_consensus_count;
my $variant_count;

my @consensus_positions;
my @consensus_chunks;

my $read_identity;
my @sample_identity;
my @protein_identity;
my @consensus_identity;
my @sample_length;

my @identity_index;

my $adjusted_identity;
my $adjusted_coverage;
my @adjusted_amino_acids;

my $partial_match;
my $partial_count;
my $partial_position_current;  

my $temp;

#Use first read file to name outputs
if ($read_1=~m/^.+\/(.*?)\.f.*?$/) { 
	$concatenated = ($1 . "_concatenated.fastq.gz");
	$diamond_output = ($1 . "_diamond_read_alignments.txt");
	$protein_cut_off = ($1 . "_protein_cut_offs.txt");

	if($read_type eq "paired") {
		#Concaternate read files
		system "cat $read_1 $read_2 > $concatenated";
		
		#Run diamond on concatenated reads
		system "diamond blastx -k 1000000000000 -e 1.0 -k0 --matrix BLOSUM45 -q $concatenated --db $db_name --threads $threads --id $identity --query-cover $coverage -o $diamond_output -f 6 stitle qseqid qseq_translated qstart qend sstart send qcovhsp pident qlen slen";
		
		#Remove concatenated reads
		system "rm $concatenated";
	} elsif ($read_type eq "single") {
		#Run diamond on single reads
		system "diamond blastx -k 1000000000000 -e 1.0 -k0 --matrix BLOSUM45 -q $read_1 --db $db_name --threads $threads --id $identity --query-cover $coverage -o $diamond_output -f 6 stitle qseqid qseq_translated qstart qend sstart send qcovhsp pident qlen slen";
	}

	#Read in diamond_output and extract read information
	@sample_proteins = ();	
#	@sample_reads = ();
	@sample_positions = ();
	@sample_amino_acids = ();
	@sample_identity = ();
	@sample_length = ();

	open DIAMOND, $diamond_output or die "File $diamond_output NOT FOUND - $!\n";
	foreach $_ (<DIAMOND>) {
		if ($_=~m/^(.*?)\t(.*?)\t(.*?)\t(\d+\.*\d*)\t(\d+\.*\d*)\t(\d+\.*\d*)\t\d+\.*\d*\t\d+\.*\d*\t(\d+\.*\d*)\t(\d+)\t(\d+)$/) {
		$read_protein = $1;
#		$read_name = $2;
		$read_translate = $3;
		$start_temp = $4;
		$end_temp = $5;
		$read_position = ($6 - 1);
		$read_identity = $7;
		$read_length = $8;
		$protein_length = $9;
			
		#Make read_start smallest and read_end largest
		if ($start_temp < $end_temp) {
			$read_start = $start_temp;
			$read_end = $end_temp;
		} else {
			$read_start = $end_temp;
			$read_end = $start_temp;
		}
			
			
		#Ignore read if alignment is in the middle of the read - account for two nucleotides due to frame shifts
		$read_alignment_QC = 0;
			
		if($read_start > 3 ){
			$read_alignment_QC++;
		}
			
		if($read_end < ($read_length - 2)){
			$read_alignment_QC++;
		}					
			
		if($read_alignment_QC < 2) {
			push @sample_proteins, $read_protein;
#			push @sample_reads, $read_name;
			push @sample_positions, $read_position;
			push @sample_amino_acids, $read_translate;	
			push @sample_identity, $read_identity;
			push @sample_length, $read_length;

			}
		}
	}

	close DIAMOND;
		
	#Identify protein sequences with read matches
	@unique_proteins = uniq @sample_proteins;
	@unique_sequences = ();

	for (my $n=0; $n < scalar(@unique_proteins); $n++) {
		for (my $o=0; $o < scalar(@protein_names); $o++) {
			if($unique_proteins[$n] eq $protein_names[$o]) {
				push @unique_sequences, $protein_sequences[$o];
			}
		}
	}
		
	#Open cut-off output file
	open(OUT, ">$protein_cut_off") or die "Couldn't open OUT $protein_cut_off $!\n";

	print OUT "Protein\tCoverage\tIdentity\tMean_read_depth\tAlignment_start\tAlignment_end\tMost_likely_sequence\n";
		
	#Loop through each protein sequence
	for (my $j=0; $j < scalar(@unique_proteins); $j++) {
		
#		@protein_reads=();
		@protein_positions=();
		@protein_amino_acids=();
		@protein_identity=();

		@protein_indexes=();

		#Loop through each read base match
		for (my $l=0; $l < scalar(@sample_proteins); $l++) {

			#Check if protein names match
			if($sample_proteins[$l] eq $unique_proteins[$j]){
				
				push @protein_indexes, $l;

				@amino_acids_split = split(//, $sample_amino_acids[$l]);
				
				@adjusted_amino_acids=();
								
				#Loop through alignment and adjust if interupted
				$partial_match = 0;
				$partial_count = 0;
				$partial_position_current =  $sample_positions[$l];
				
				#Trim read alignments if there is stop codon
				for (my $z=0; $z < scalar(@amino_acids_split); $z++) {
					#Check if read amino acid position is within reference amino acid alignment
					if(length($unique_sequences[$j]) < $partial_position_current ) {
					} else {
						if($amino_acids_split[$z] eq "*"){
							push(@adjusted_amino_acids, $amino_acids_split[$z]);
							$partial_count++;
							
							if($amino_acids_split[$z] eq (substr($unique_sequences[$j],$partial_position_current, 1))){
								$partial_match++;
							}
							$partial_position_current++;
							last;
						}elsif($amino_acids_split[$z] eq "X"){
							push(@adjusted_amino_acids, $amino_acids_split[$z]);
							$partial_count++;
							
							if($amino_acids_split[$z] eq (substr($unique_sequences[$j],$partial_position_current, 1))){
								$partial_match++;
							}
							$partial_position_current++;
							last;
						} else {
							push(@adjusted_amino_acids, $amino_acids_split[$z]);
							$partial_count++;

							if($amino_acids_split[$z] eq (substr($unique_sequences[$j],$partial_position_current, 1))){
								$partial_match++;
							}
							$partial_position_current++;
						}
					}
				}
				
				#Adjust read coverage and identity
				$adjusted_identity = ($partial_match / $partial_count) * 100;
				
				$adjusted_coverage = (($partial_count * 3)/ $sample_length[$l]) * 100;
				
				#Extract amino acids from reads that pass QC
				if($adjusted_identity >= $identity) {
					if($adjusted_coverage >= $coverage) {
						$sample_position_current = $sample_positions[$l];
						
						for (my $i=0; $i < scalar(@adjusted_amino_acids); $i++) {
#							push @protein_reads, $sample_reads[$l];
							push @protein_positions, $sample_position_current;
							push @protein_amino_acids, $adjusted_amino_acids[$i];
							push @protein_identity, $adjusted_identity;
						
							$sample_position_current++;
						}
					}
				}
			}
		}

		#Remove protein matched positions
		if(scalar(@protein_indexes) > 0) {
			my @protein_index_sorted = sort { $b <=> $a } @protein_indexes;

			for (my $r=0; $r < (scalar(@protein_index_sorted)); $r++) {
				splice @sample_proteins, $protein_index_sorted[$r], 1;
#				splice @sample_reads, $protein_index_sorted[$r], 1;
				splice @sample_positions, $protein_index_sorted[$r], 1;
				splice @sample_amino_acids, $protein_index_sorted[$r], 1;
				splice @sample_identity, $protein_index_sorted[$r], 1;
			}
		}

		#Loop through each amino acid of protein sequence
		@protein_split = split(//, $unique_sequences[$j]);

		@position_depth_array = ();
		@position_match_array = ();
		@position_cover_array = ();
		@position_start_array = ();
		@position_end_array = ();
		@consensus_chunks = ();
		@consensus_positions = ();
		
		$depth = 0;
		$match = 0;
		$cover = 0;

		#Set alignment_start to a very large value
		$alignment_start = 1000000000000000000;
		
		for (my $k=0; $k < scalar(@protein_split); $k++) {

			$position_depth = 0;

			@indexes = ();					
		
			$alignment_end = $k;

			@consensus = ();
			@consensus_identity = ();
		
			#Loop through each protein position
			for (my $p=0; $p < scalar(@protein_positions); $p++) {	
				
				#Check if positions match
				if($protein_positions[$p] == $k){
				push @indexes, $p;

				$depth++;
				$position_depth++;
				
				#Store amino acid at position for consensus
				push(@consensus, $protein_amino_acids[$p]);
				push(@consensus_identity, $protein_identity[$p]);
				
				}
			}

			if($sequence_identification eq "consensus") {
				#Determine most common amino acid for position
				if (scalar(@consensus) > 0) {
					@uniq_consensus_count = ();
				
					@uniq_consensus = uniq(@consensus);
				
					for (my $p=0; $p < (scalar(@uniq_consensus)); $p++) {

						$variant_count = 0;

						for (my $q=0; $q < (scalar(@consensus)); $q++) {
							if($uniq_consensus[$p] eq $consensus[$q]) {
								$variant_count++;
							}
						}

						push(@uniq_consensus_count, $variant_count);
					}

					#Sort count from largest to smallest
					my @sorted_uniq_consensus_count = sort { $uniq_consensus_count[$b] <=> $uniq_consensus_count[$a] } 0..$#uniq_consensus_count;

					@uniq_consensus = @uniq_consensus[@sorted_uniq_consensus_count];
			
					push (@consensus_positions, $uniq_consensus[0]);

					#Determine if consensus and reference position match
					if($protein_split[$k] eq $uniq_consensus[0]){
						$match++;
					}	
				}
			} elsif ($sequence_identification eq "identity") {

				#Determine amino acid for position based on the best identity for the read
				if (scalar(@consensus) > 0) {

					#Sort read positions by read identity
					@identity_index = sort { $consensus_identity[$b] <=> $consensus_identity[$a] } 0 .. $#consensus_identity;

					@consensus_identity = @consensus_identity[@identity_index];
					@consensus = @consensus[@identity_index];

					push (@consensus_positions, $consensus[0]);
				
					#Determine if consensus and reference position match
					if($protein_split[$k] eq $consensus[0]){
						$match++;
					}	
				}
			} else {
				if (scalar(@consensus) > 0) {
					#If unknown consensus method is given, leave "-" as consensus position
					push (@consensus_positions, "-");
				}
			}
				
			#Remove matched positions
			if(scalar(@indexes) > 0) {
				my @index_sorted = sort { $b <=> $a } @indexes;

				for (my $m=0; $m < (scalar(@index_sorted)); $m++) {
#					splice @protein_reads, $index_sorted[$m], 1;
					splice @protein_positions, $index_sorted[$m], 1;
					splice @protein_amino_acids, $index_sorted[$m], 1;
					splice @protein_identity, $index_sorted[$m], 1;
				}
			}
		
			
		
			if($position_depth == 0) {
				push @position_match_array, $match;
				push @position_cover_array, $cover;
				push @position_depth_array, $depth;
				push @position_start_array, ($alignment_start + 1);
				push @position_end_array, $alignment_end;
				push @consensus_chunks, join("", @consensus_positions);
					
				$match = 0;
				$cover = 0;
				$depth = 0;
				
				@consensus_positions=();
				
				$alignment_start = 1000000000000000000;
			} elsif($consensus_positions[-1] eq "*") {
				push @position_match_array, $match;
				push @position_cover_array, $cover;
				push @position_depth_array, $depth;
				push @position_start_array, ($alignment_start + 1);
				push @position_end_array, $alignment_end;
				push @consensus_chunks, join("", @consensus_positions);
					
				$match = 0;
				$cover = 0;
				$depth = 0;
				
				@consensus_positions=();
				
				$alignment_start = 1000000000000000000;
			} elsif($consensus_positions[-1] eq "X") {
				push @position_match_array, $match;
				push @position_cover_array, $cover;
				push @position_depth_array, $depth;
				push @position_start_array, ($alignment_start + 1);
				push @position_end_array, $alignment_end;
				push @consensus_chunks, join("", @consensus_positions);
					
				$match = 0;
				$cover = 0;
				$depth = 0;
				
				@consensus_positions=();
				
				$alignment_start = 1000000000000000000;
			} else {
				$cover++;
				if($k < $alignment_start) {
					$alignment_start = $k;
				}
			}
		}

		push @position_match_array, $match;
		push @position_cover_array, $cover;
		push @position_depth_array, $depth;
		push @position_start_array, $alignment_start;
		push @position_end_array, ($alignment_end - 1);
		push @consensus_chunks, join("", @consensus_positions);

		#Sort position arrrays from largest to smallest fragment
		@index_sorted_size = sort { $position_cover_array[$b] <=> $position_cover_array[$a] } 0 .. $#position_cover_array;

		@position_cover_array = @position_cover_array[@index_sorted_size];
		@position_match_array = @position_match_array[@index_sorted_size];
		@position_depth_array = @position_depth_array[@index_sorted_size];
		@position_start_array = @position_start_array[@index_sorted_size];
		@position_end_array = @position_end_array[@index_sorted_size];
		@consensus_chunks = @consensus_chunks[@index_sorted_size];

		#Calcular coverage, identity and mean read depth from largest fragment
		if($position_cover_array[0] == 0) {
			$identity_proportion = 0;
		} else {
			$identity_proportion = $position_match_array[0] / $position_cover_array[0];
		}

		$coverage_proportion = $position_cover_array[0] / scalar(@protein_split);
		$mean_read_depth = $position_depth_array[0] / scalar(@protein_split);
		
		if($coverage_proportion > 0){
			print OUT "$unique_proteins[$j]\t$coverage_proportion\t$identity_proportion\t$mean_read_depth\t$position_start_array[0]\t$position_end_array[0]\t$consensus_chunks[0]\n";
		}
	}
		
	close OUT;
} else {
	die "Please provide reads in fastq format\n";
}
