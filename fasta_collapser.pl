#!/usr/bin/env perl
#	#!/usr/local/bio_apps/perl-5.16.2/bin/perl

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
use Bio::SeqIO;
#use Getopt::Std;
#use PostData;
#use Fasta_utils;
#use feature ":5.10";	#allows say (same as print, but adds "\n"), given/when switch, etc.


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
if (@ARGV){		print STDERR "Arguments: ", join " ", @ARGV, "\n";	}

my $save;
my $files;
my $verbose;
my $output;
my $gzip;
GetOptions('save=s' => \$save, 'output=s' => \$output, 'o=s' => \$output, 'verbose' => \$verbose, 'files=s' => \$files, 'i=s' =>\$files, 'gzip' => \$gzip);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $exe = basename($0);

my $usage = "$exe is my Perl implementation of fastx_collapser from 
the FASTX toolkit.  This version is just for fasta files and allows ambiguous
bases.  RAM usage is approximately equivalent to the input file size.

$exe [OPTIONS] -i <fasta> -o <output fasta>
OPTIONS:
-i/--files 	Files, comma-delimited, or a directory of files. Required.
-o/--output	Output file name.   Required (for now).
";




# Change log
# 2014-06-12
# Created script.  Purpose of this script is to do the same thing that fastx_collapser (FASTX Toolkit) does for fasta files, but also permits ambiguous bases.  merge_primerid_read_groups.pl sometimes produces sequences with ambiguous bases, e.g., when there is a non-N ambiguous base in the middle gap in merge_primerid_read_groups.pl, it is not treated as ambiguous so the sequence is retained.  fastx_collapser dies with an error when an ambiguous base in encountered.  
# From my notes: "I wrote a script called fasta_collapser.pl.  It's pretty fast.  For the file Sample_6.contigs.pid10bp.btrim.1.majority.fasta (342562 sequences; 162MB), it takes 2.666 seconds for fastx_collapser and 17.362 sec for fasta_collapser.pl, so about 6.5-fold slower, but that's okay.  Gave the same results.  My Perl script uses memory relatively efficiently; at the peak, it is about 150MB, so that's good."
# This could be a subroutine in convert_reads_to_amino_acids.pl


unless ($files||$ARGV[0]){	print STDERR "$usage\n";	exit;	}
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
my @suffixes = (qw(.bed .bed.gz .bed12 .bed12.gz .txt .txt.gz .BED .BED.gz .BED12 .BED12.gz .fasta .fa .FA .FASTA .FAS .fas), @SUFFIXES);	#for fileparse.  Feel free to add more accepted extensions.  @SUFFIXES comes from aomisc.pm.  

foreach my $file (@files){
	read_and_write($file);
}




&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub read_and_write {
	my $file = shift;

	# Initialize hashref to count the sequences.  
	my $sequences;  
	
	# Read the fasta file, store the sequences. 
	my $in  = Bio::SeqIO->new(-file => "$file" ,		
                           -format => 'Fasta');			# http://search.cpan.org/~cjfields/BioPerl-1.6.922/Bio/SeqIO.pm
	
	my $total = 0;
	FASTA: while ( my $seqobj = $in->next_seq() ) {
		# e.g., 
		#	>MISEQ:50:000000000-A4142:1:1101:10000:4300:ACGACTTGCA
		#	ATTCGAAAGATTCAAAATATTTCCCAAAGAAAGCTCATGGCCCGACCACAACACAAACGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCNNNNNNNTAACAGTAAGGAACAACAGGATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
		$total++;
		my $seq = $seqobj->seq();
		$sequences->{$seq}++;
	}
	print STDERR "Read $total sequences\n";
	
	# Print sequences in order from highest to lowest count and print them out.  Give each an id based on order and count e.g., >1-15809 for the first, which has 15809 sequences in the original fasta file.
	my $writefh = open_to_write($output,0,0,1);
	my $count = 0;
	foreach my $seq (sort { $sequences->{$b} <=> $sequences->{$a} } keys %$sequences){
		$count++;
		print $writefh ">$count-".$sequences->{$seq}."\n$seq\n";
#		print ">$count-".$sequences->{$seq}."\n$seq\n";
	}
	print STDERR "Printed $count unique sequences.\n";
	
}
#-----------------------------------------------------------------------------

