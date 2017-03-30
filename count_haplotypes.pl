#!/usr/bin/env perl 
use warnings;
$| = 1;

#Add use lib statement to assume there is a directory at the same level as bin in which the script is run, called 'lib'
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin";

use strict;
use FileHandle;
use aomisc;
use Cwd;
use diagnostics; 
use Getopt::Long;
use Data::Dumper;
use File::Basename;


# Andrew J. Oler, PhD
# Computational Biology Section
# Bioinformatics and Computational Biosciences Branch (BCBB)
# OCICB/OSMO/OD/NIAID/NIH
# Bethesda, MD 20892
#
# andrew.oler@gmail.com
# 
#This package is free software; you can redistribute it and/or modify
#it under the terms of the GPL (either version 1, or at your option,
#any later version) or the Artistic License 2.0.


#my %options;	#Hash in which to store option arguments
#use vars qw($opt_s $opt_f $opt_n $opt);
#getopts ('s:f:n:bda:g:om:i:r:ecu:kh:j:l:ypx:',\%options);	

#Print out the options
print STDERR "Arguments: ", join " ", @ARGV, "\n" if (@ARGV); 

my $save;
my $files;
my $verbose;
my $output;
my $gzip;
my $pos;
my $hap_min = 2;
my $skip = 'X';
GetOptions(
	'save=s' => \$save, 
	'output=s' => \$output, 
	'verbose' => \$verbose, 
	'files=s' => \$files, 
	'gzip' => \$gzip,
	'pos=s' => \$pos,
	'hap_min=s' => \$hap_min,
	'skip=s' => \$skip,
);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $exe = basename($0);

my $usage = "$exe [options] --pos <positions> --files <files> 
This script takes translated reads (output of convert_reads_to_amino_acid.pl), 
extracts amino acids at a particular position and tallies the motifs up across
one or multiple files.

Example input:
#readID|gene|positionForFirstAminoAcid	cleanPeptide
4-866|HA|143 HNTTGVTAACSHEGKNSFYRNLLWLTKKESSYPELKNSYVNKKRKEVLVLWGIHHPPNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVKGQAGXMNYYWTLLKPGDTIIFEANGNLIAPMYAFALSRGFGSGIITSNASMHECNTKCQTPLG
1-21583|HA|143 HNTTGVTAACSHEGKNSFYRNLLWLTKKESSYPELKNSYVNKKRKEVLVLWGIHHPPNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVKGQAGRMNYYWTLLKPGDTIIFEANGNLIAPMYAFALSRGFGSGIITSNASMHECNTKCQTPLG
2-2631|HA|143 HNTNGVTAACSHEGKNSFYRNLLWLTKKESSYPELKNSYVNKKRKEVLVLWGIHHPPNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVKDQAGRMNYYWTLLKPGDTIIFEANGNLIAPMYAFALSRGFGSGIITSNASMHECNTKCQTPLG
...

Example output:
motif	141_73_S5	151_71_S3	152_70_S3	151_78_S2
NTKNG	103	28	33	81257
NTENG	7706	3828	3205	14697
TTKNG	50587	27490	37243	2553
...

OPTIONS:
--files 	Files, comma-delimited, or a directory of files.  Takes *.majority.cons.uniq.cleanpeptides.txt
		file(s) as input.  Files can also be provided as arguments on the command-line.
--pos 		1-based amino acid positions, comma-delimited.  Required.
--hap_min	minimum count in at least one file for a particular haplotype to include in the output. 
		Default = 2.  
--skip	Skip motifs containing this string. Default = 'X' (unknown amino acids)
";

# Change log
# 2017-03-30
# Created script.


# To do

unless ($pos && ($files||$ARGV[0])){	print STDERR "$usage\n";	exit;	}
unless($files){	
	if (scalar(@ARGV)>1){	#This will allow you to use wildcard to pass files to the script from the command line. e.g, script.pl *.txt, and it will run through each.
		$files = join ",", @ARGV;
	}
	else{
		$files = $ARGV[0];	
	}
}

my $start_time = time;
my @files = &aomisc::get_files($files);		#If allowing a directory, specify extension of the files in second argument, e.g., my @files = &get_files($files, 'bed');

my $hash; 

foreach my $file (@files){
	(my $base = $file ) =~ s/\.contigs.pid.btrim.\d+.majority.cons.uniq.cleanpeptides.txt//;
	$hash->{$base} = read_and_write($file, $pos, $skip);
}


print_table($hash, $hap_min);






&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub read_and_write {
	my ($file, $pos, $skip) = @_;

	# Get the index position within the sequence
	my @pos = split(/,/, $pos);
	my @index;
	foreach my $position (@pos){
		push @index, $position - 1;
	}

	my $readfh = open_to_read($file);
	
	my $data;	# Hashref of the data for this file.
	my $haplotypes;	# Hashref of all haplotypes

	while (<$readfh>){
		# head 152_70_S3.contigs.pid.btrim.1.majority.cons.uniq.cleanpeptides.txt
		#	  9682 158
		#	#readID|gene|positionForFirstAminoAcid	cleanPeptide
		#	4-866|HA|143 HNTTGVTAACSHEGKNSFYRNLLWLTKKESSYPELKNSYVNKKRKEVLVLWGIHHPPNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVKGQAGXMNYYWTLLKPGDTIIFEANGNLIAPMYAFALSRGFGSGIITSNASMHECNTKCQTPLG
		#	1-21583|HA|143 HNTTGVTAACSHEGKNSFYRNLLWLTKKESSYPELKNSYVNKKRKEVLVLWGIHHPPNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVKGQAGRMNYYWTLLKPGDTIIFEANGNLIAPMYAFALSRGFGSGIITSNASMHECNTKCQTPLG
		#	2-2631|HA|143 HNTNGVTAACSHEGKNSFYRNLLWLTKKESSYPELKNSYVNKKRKEVLVLWGIHHPPNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVKDQAGRMNYYWTLLKPGDTIIFEANGNLIAPMYAFALSRGFGSGIITSNASMHECNTKCQTPLG
		chomp;
		my @F = split(/\s+/);
		next unless ($F[0] =~ m/HA/);

		my @id = split(/\|/, $F[0]);
		my ($rank,$count) = split(/-/, $id[0]);	# Get the multiplier

		my @seq = split(/|/, $F[1]);	 #amino acid sequence
		my $hap = "";
		foreach my $index (@index){
			$hap .= $seq[$index]; 
		}
		$data->{$hap} += $count unless $hap =~ m/$skip/;
	}
	
	#	close ($writefh);
	close($readfh);
	
	return $data;
}
#-----------------------------------------------------------------------------
sub print_table {
	my ($hash,$hap_min) = @_;

	my @files = sort keys %$hash; 

	# Get haplotypes passing $hap_min threshold
	my %hap;
	foreach my $file (@files){
		foreach my $hap (keys %{$hash->{$file}}){
			if ($hash->{$file}->{$hap} >= $hap_min){
				$hap{$hap}++;
			}
		}
	}
	my @hap_pass = keys %hap;

	# print header
	print join "\t", "motif", @files; 
	print "\n"; 	
	
	foreach my $hap (@hap_pass){
		print "$hap"; 
		foreach my $file (@files){
			my $count = exists($hash->{$file}->{$hap}) ? $hash->{$file}->{$hap} : 0; 
			print "\t$count"; 
		}
		print "\n";
	}
}
#-----------------------------------------------------------------------------
# for i in *_7*btrim.1.*pep*; do echo $i; cat $i | awk '$1 ~ /HA/{print $2}' | cut -c 4,7,15,16,96 | grep -v X | sort | uniq -c | sort -k2,2 | egrep "(NMKNG|NTEDG|NTENG|NTKDG|NTKND|NTKNG|NTKSG|NTKTG|NTTNG|STENG|STKNG|TMKNG|TSKNG|TTENG|TTKNG|TTKNS|TTKTG)"; done | perl -e'my $hash; my $file; my %hap; while(<>){chomp; if (m/contigs/){$file = $_; next;}else { s/^\s+//; my @F = split(/\s+/); $hash->{$file}->{$F[1]} = $F[0]; $hap{$F[1]}++; } } my @files = keys %$hash; print join "\t", "motif", @files; print "\n"; foreach my $hap (keys %hap){ print "$hap"; foreach my $file (@files){ my $count = exists($hash->{$file}->{$hap}) ? $hash->{$file}->{$hap} : 0; print "\t$count"; }print "\n"; }' | sed 's/.contigs.pid.btrim.1.majority.cons.uniq.cleanpeptides.txt//g'
#-----------------------------------------------------------------------------

