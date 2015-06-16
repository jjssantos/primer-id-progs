#!/usr/bin/env perl

use strict;
use warnings;
$| = 1;
use aomisc;
use Data::Dumper;
use File::Basename;

my $exe = basename($0);

my $usage = "$exe <input aa/nuc.tally.xls>
Prints to STDOUT.  Takes aa or nuc tally.xls file as input (not codon).
";

die $usage unless $ARGV[0];

my $readfh = open_to_read($ARGV[0]);

my @header;
my $residue_to_header_index; # hashref, keys residues (nuc/aa), values, index in the line to get the count
LINE: while(<$readfh>){
	chomp;
	my @line = split(/\t/);
	if ((m/^#/)|( m/name/i && m/consensus/i && m/ref/i && m/coverage/i) ){		# Sometimes the # is gone, in case of compare_variant_frequency.pl where the file is going into R, which doesn't like # in the header.
		@header = @line;
		$residue_to_header_index = map_residues_to_header_indexes(\@header);
		print join "\t", @line, "numConsensus", "numNonConsensus", "numMajorAltAllele", "majorAltAllele"; 
		print "\n";
		next LINE;
	} 
	else {
		my $consensus_residue = $line[4];
		my ($consensus_count,$nonconsensus_count,$major_alt_allele,$major_alt_allele_count) = (0,0,"",0);
		
		# If it's a normal residue, get consensus residue count and tally up nonconsensus residue counts and 
		if (exists ($residue_to_header_index->{$consensus_residue})){		# Only if it is one of the regular aa/nuc residues (not X, -, etc.)
			my $index = $residue_to_header_index->{$consensus_residue};
			$consensus_count = $line[$index];
			foreach my $residue (keys %$residue_to_header_index){
				unless ($residue eq $consensus_residue){		# Don't add the consensus base/aa 
					my $this_count = $line[$residue_to_header_index->{$residue}];
					$nonconsensus_count += $this_count; 	# Add all other bases/aa's
					($major_alt_allele,$major_alt_allele_count) = ($residue, $this_count) if ($this_count > $major_alt_allele_count);
				}
			}
		}
		
		# Print out the data.
		print join "\t", @line, $consensus_count, $nonconsensus_count, $major_alt_allele_count, $major_alt_allele; 
		print "\n";		
	}
}




##-----------------------------------
sub map_residues_to_header_indexes {
	my $header = shift;	# arrayref
	my $map;	# hashref where keys are residues (aa or nuc) and value is the index in the header to find this residue
	for (my $i = 0; $i < @$header; $i++){
		if ( (@header < 15 && $header->[$i] =~ m/^num([A-MO-WY])$/) || (@header > 15 && $header->[$i] =~ m/^num([A-WY])$/) ){		# All amino acids and nucleotides; not including "X", "*", "-", or "Other".  "N" is included for aa but not for nuc.  
			$map->{$1} = $i;
		} 
	}	
	return $map;
# Example
#$VAR1 = {
#          'A' => 8,
#          'T' => 11,
#          'C' => 9,
#          'G' => 10
#        };

}