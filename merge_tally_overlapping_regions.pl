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
#use Bio::SeqIO;


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
my $label = "Sample_1";
GetOptions(
	'save=s' => \$save, 
	'output=s' => \$output, 
	'verbose' => \$verbose, 
	'files=s' => \$files, 
	'gzip' => \$gzip,
	'label=s' => \$label,
);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $exe = basename($0);

my $usage = "$exe [options] --files <region1.tally.xls> <region2.tally.xls>
$exe takes output from convert_reads_to_amino_acid.pl and merges frequency counts for 
positions where the amplicons overlap.  Right now, only includes a subset of the columns 
in the output, necessary for making a scatterplot, i.e., name, Sample, Position, numNonConsensus,
coverageDepth.
OPTIONS:
--files 	Files, comma-delimited, or a directory of files.  File names can also be 
			received as arguments
--output	Output file.  Prints to STDOUT by default.  
--label		Sample ID.  All files should be for the same sample. Default is 'Sample_1'.
";

# Change log
# 2015-06-16
# Created script.
# 2016-09-12
# Modified column header lookups to match the new column headers in convert_reads_to_amino_acid.pl (instead of the columns output by add_consensus_columns_to_frequency_tables.pl)


# To do

unless ($files||$ARGV[0]){	print STDERR "$usage\n";	exit;	}
unless($files){	
	if (scalar(@ARGV)>1){	#This will allow you to use wildcard to pass files to the script from the command line. e.g, script.pl *.txt, and it will run through each.
		$files = join ",", @ARGV;
	}
	else{
		$files = $ARGV[0];	
	}
}

my $save_dir = $save || Cwd::cwd();
unless (-d $save_dir){	mkdir($save_dir) or warn "$!\n"	}
warn "save directory: $save_dir\n";

my $start_time = time;
my @files = &aomisc::get_files($files);		#If allowing a directory, specify extension of the files in second argument, e.g., my @files = &get_files($files, 'bed');
my @suffixes = (qw(.bed .bed.gz .bed12 .bed12.gz .txt .txt.gz .BED .BED.gz .BED12 .BED12.gz .fasta .fa .FA .FASTA .FAS .fas), @SUFFIXES);	#for fileparse.  Feel free to add more accepted extensions.  @SUFFIXES comes from aomisc.pm.  

my $data;	# HoHoH where first key is the position, second keys are name, coverageDepth, numNonConsensus, merged.  As files are being read, coverageDepth and numNonConsensus will be incremented.  

foreach my $file (@files){
	read_tally($file);
}

# Print out the table.
my $writefh = 
	($output) ? open_to_write($output) :
		*STDOUT;

print $writefh join "\t", "sample", "name", "position", "coverageDepth", "unambigCoverageDepth", "numUnambigNonConsensus", "merged";
print $writefh "\n";
foreach my $pos (sort {$a <=> $b} keys %$data){
	print $writefh join "\t", $label, $data->{$pos}->{name}, $pos, $data->{$pos}->{coverageDepth}, $data->{$pos}->{unambigCoverageDepth}, $data->{$pos}->{numUnambigNonConsensus}, $data->{$pos}->{merged};
	print $writefh "\n";
}

close($writefh) if ($output);	# Don't close STDOUT

&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub read_tally {
	my $file = shift;
	
	my @header = get_header($file);
	$header[0] =~ s/^#//;
	my $h_lookup = column_header_lookup_hash(\@header);
	my $name_index 			= $h_lookup->{'name'};
	my $cov_index 			= $h_lookup->{'coverageDepth'};
	my $noncon_index 		= $h_lookup->{'numUnambigNonConsensus'};
	my $unambig_cov_index 	= $h_lookup->{'unambigCoverageDepth'};
	my $pos_index = 															# aminoAcidPosition OR nucleotidePosition OR codonPosition
		(exists($h_lookup->{aminoAcidPosition}))		? $h_lookup->{aminoAcidPosition} :
		(exists($h_lookup->{nucleotidePosition}))	? $h_lookup->{nucleotidePosition} :
		(exists($h_lookup->{codonPosition}))		? $h_lookup->{codonPosition} :   
			2;	# default
	
	
	my $readfh = open_to_read($file);
	while (<$readfh>){
		chomp;
		next if (m/coverageDepth/);		 #skip header column
		my @F = split(/\t/);
		my ($name, $pos, $cov, $noncon, $unambig_cov) = ($F[$name_index], $F[$pos_index], $F[$cov_index], $F[$noncon_index], $F[$unambig_cov_index]);

		$data->{$pos}->{name} = $pos unless exists($data->{$pos}->{name});
		if (exists($data->{$pos}->{coverageDepth})){
			$data->{$pos}->{merged} = "T";
		}
		else {
			$data->{$pos}->{merged} = "F";
		}
		
		$data->{$pos}->{coverageDepth} += $cov;
		$data->{$pos}->{numUnambigNonConsensus} += $noncon;
		$data->{$pos}->{unambigCoverageDepth} += $unambig_cov;
	}
	
	close($readfh);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#Useful things

#my $total_lines = 0;
#	$total_lines++;		#somewhere in subroutine, if I run through the file first.
#	my $lines_done = 0;
#	$lines_done++;		#in while loop, increment as each new line is processed
#	if (($lines_done % 25) == 0){	#put this in the while loop too.
#		print "Done processing $lines_done of $total_lines sequences.";
#		&elapsed($start_time, ' Elapsed', $verbose);
#	}
