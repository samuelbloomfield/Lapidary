###
#	Random_amino_acid_alterations
#	Script to alter amino acid alterations
#	The script can be run by: perl Random_amino_acid_alterations.pl -input <Location of amino acid sequences in fasta format> -coverage <coverage to reduce amino acid sequences to> -identity <identity to change amino acid sequences to>
#	Last modified: 11/07/2024
#	Author: Samuel Bloomfield
###


use warnings;
use strict;
use Getopt::Long;
use List::MoreUtils qw[uniq];
use POSIX;


# here check if no arguments were given, then print help message
if (@ARGV == 0) {
	# Print the help message
	print "\n\tThis script requires some arguments to run.\n\n";
	print "\tUse --help to see the options\n";
	exit;
}


my $input;
my $identity;
my $coverage;
my $help;


GetOptions (	'input:s'			=> \$input,
				'identity:s'		=> \$identity,
				'coverage:s'		=> \$coverage,
				'help'				=> \$help);

#Print out help message
if(defined($help)){
	die "\n\tRandom_amino_acid_alterations: software for changing the coverage and identity of amino acid sequences\n\n
	Options:\n
	input\t\tLocation of fasta file\n
	coverage\tCoverage to change amino acid sequences to as a proportion\n
	identity\tIdentity to change amino acid sequences to as a proportion\n
	help\t\tDisplay help screen\n";
} 

#Array of potential aminoc acid codes
my @amino_acids = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y");


#Read in fasta file and place sequences and sequence names into arrays
my @substring=();
my @contigs=();
my @sequences=();

my $line;

open FILE, $input or die "File $input NOT FOUND - $!\n";

foreach $_ (<FILE>) {
	$line = $_;
	if ($line=~m/^>(.*?)$/) {
		push @contigs,$1;
		if (scalar(@substring)>0){
			push @sequences, join("", @substring);
			@substring=();
		}
	} elsif ($line=~m/^(.*?)$/) {
		push @substring, $1;
	}
}

push @sequences, join("", @substring);

close FILE;


my $output;

my $protein_length;
my $protein_reduced_length;

my $sequence_reduced;
my $amino_acids_to_alter;
my $sequence_altered;

my $contig_altered;

my @sequence_array;
my @positions_to_change;

my @unique_array;

my $amino_acid_position;
my $amino_acid_replacement;

my $different_amino_acid;


#Output results to output file sepecific to coverage and identity
$output = ($input . "_coverage_" . $coverage . "_identity_" . $identity . ".fasta");
	
open(OUT, ">$output") or die "Couldn't open OUT $output $!\n";

#Loop through sequences
for(my $i=0; $i < scalar(@contigs); $i++) {
	
	#Rename sequence names to include coverage and identity values
	$contig_altered = ($contigs[$i] . "_coverage_" . $coverage . "_identity_" . $identity);
	
	#Calculate protein length using coverage cut-off - round up to integer
	$protein_length = length($sequences[$i]);
	$protein_reduced_length = ceil($protein_length * $coverage);
	
	#Extract first section of protein sequence
	$sequence_reduced = substr $sequences[$i], 0, $protein_reduced_length;
	
	#Create array of each sequence
	@sequence_array = split(//, $sequence_reduced);
	
	#Calculate amino acids to alter using identity - round down
	$amino_acids_to_alter = floor($protein_reduced_length * (1 - $identity));
	
	#Make array of positions to change
	@positions_to_change = ();
	
	#Determine if any should be change
	if($amino_acids_to_alter < 1) {
		$sequence_altered = $sequence_reduced;
	} else {
	
		do {
			#Randomly choose position
			push @positions_to_change, floor(rand($protein_reduced_length));
							
			#Remove duplicates
			@unique_array = uniq(@positions_to_change);
			
			@positions_to_change = @unique_array;
			
		} until(scalar(@positions_to_change)>($amino_acids_to_alter - 1));
		
		
		#Loop through positions to change
		for(my $m=0; $m < scalar(@positions_to_change); $m++) {
			$amino_acid_position = $sequence_array[$positions_to_change[$m]];
			
			#Find replacement that is different
			$different_amino_acid = 0;
			
			do{
				$amino_acid_replacement = $amino_acids[rand(@amino_acids)];
				
				if($amino_acid_replacement eq $amino_acid_position) {
				} else {
					$different_amino_acid++;
				}
				
			} until ($different_amino_acid > 0);
			
			$sequence_array[$positions_to_change[$m]] = $amino_acid_replacement;
		}
		
		#Convert altered amino acids into an array
		$sequence_altered = join("", @sequence_array);
		
	}
	
	#Print altered protein sequences to file
	print OUT ">$contig_altered\n$sequence_altered\n"
}