#!/usr/bin/env perl
use strict;
use Data::Dumper; 
use File::Basename; 
use aomisc; 

my $usage="This script takes the output of Btrim64 and checks for the number of reads to decide which files to continue with in the pipeline
For example, if Btrim64 is supplied with a primers file with 3 sets of primers, and finds 100, 500000, and 1000000 matches to the three primer sets, respectively, it is most likely that the matches to the first primer set are false positives since they represent such a small fraction.
The output is a list of the file extensions for the files that passed, e.g., if the btrim output files are file.fastq.0, file.fastq.1 and file.fastq.2 then it would output this:
\"1 2 \"
Use in a for loop such as this:
base=file.fastq; for i in \$(check_btrimmed_read_count.pl \${base}*); do f=\${base}.\${i}; bwa mem genome.fa \$f > \${base/.fastq/}.sam; done
";

unless (@ARGV){
	print STDERR "$usage\n";
	exit;
}
my $min_fraction = 0.05;		# Required setting.  Largest file (most reads) will be multiplied by this value to get the minimum # of reads required to consider the file.  Default = 0.05, or 5% of the largest file

my $hash; 

# Get possible suffixes
# Example file names: Sample_5.concat.primerid10bp.btrim.fastq.0, Sample_5.concat.primerid10bp.btrim.fastq.1, Sample_5.concat.primerid10bp.btrim.fastq.2
my @a = (0..scalar(@ARGV) + 10); 
my @suffixes = map {$_ = "." . $_} @a; 

# Save the number of reads in each file with key as the suffix.
foreach my $file (@ARGV){ 
#	print STDERR "Reading $file\n"; 
	my ($fn,$dir,$ext) = fileparse($file,@suffixes); 
	$ext =~ s/\.//; 
	unless($ext =~ m/^\d+$/){
		print STDERR "Problem with extension: $ext\nIs this file from Btrim64?\n$file\n"; 
		exit; 
	}
	my $wc = `wc -l $file`; 
	$wc = $wc / 4; 
	$hash->{$ext} = $wc; 
} 

my $max_value = find_key_with_biggest_value($hash,1); 
my $min_value = $max_value * $min_fraction; 	

# Walk through each and print out if it passes the minimum threshold
foreach (sort {$a <=> $b} keys %$hash){ 
	if ($hash->{$_} >= $min_value){
		print "$_ "; 
	} 
} 
