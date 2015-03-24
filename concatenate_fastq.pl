#!/usr/local/bio_apps/perl-5.16.2/bin/perl

use warnings;
# select STDOUT;		# When printing to STDOUT from multiple threads/forks/processes at a time, you need to make STDOUT "hot" instead of the default buffered.  When buffered, the output can be garbled and intermingled, especially when redirecting STDOUT to a file.  See http://www.perlmonks.org/?node_id=619092 and http://perl.plover.com/FAQs/Buffering.html .   # This by itself slows things down by about 2-fold
$| = 1;

#Add use lib statement to assume there is a directory at the same level as bin in which the script is run, called 'lib'
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin";

use strict;
use FileHandle;
use aomisc;
#use fastq_clustal;
use Cwd;
use diagnostics; 
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use File::Temp;


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

# Some code and ideas for reading two fastq files at the same time borrowed from http://www.dylanstorey.com/node/50

#Print out the options
if (@ARGV){		print STDERR "Arguments: ", join " ", @ARGV, "\n";	}

my $save = Cwd::cwd();
my $files;
my $verbose;
my $output;
my $gzip;
my $number = 0;
my $ref; 
my $buffer = 20;
my $cpu	= 1;
my $char = "N"; 
my $clustalw;
GetOptions('save=s' => \$save, 'output=s' => \$output, 'o=s' => \$output, 'verbose' => \$verbose, 'files=s' => \$files, 'gzip' => \$gzip, 'number=s' => \$number, 'n=s' => \$number, 'ref=s' => \$ref, 'r=s' => \$ref, 'buffer=s' => \$buffer, 'b=s' => \$buffer, 'cpu=s' => \$cpu, 'p=s' => \$cpu, 'char=s' => \$char, 'c=s' => \$char, 'clustalw' => \$clustalw, 'w' => \$clustalw);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = basename($0)." [options] <fastq1> <fastq2>
Takes paired-end fastq read files and concatenates the R2 read onto the R1 read and the R2 
quality scores onto the R1 quality scores.  R2 sequence is reverse complemented before 
joining to R1; R2 quality scores are reversed before joining to R1.  Assumes there is no 
overlap between the R1 and R2 reads.  (If you expect there could be overlap, even a few 
bp, use something like pandaseq or pear to get concatenated/joined reads.)  
fastq1 is the first read in a pair.
fastq2 is the second read in a pair.
Processes about 1000 pairs/second with 12 processors (if no other instances running) using 
a reference.  5x faster without using a reference.
OPTIONS:
-n/--number	Number of Ns to use in concatenation. Default = 0.
-r/--ref	Single fasta file containing reference sequence to use to decide a custom 
	number of Ns to insert in concatenation.  R1 and R2 reads will be aligned to the 
	reference with clustalw	and the internal gap size will be determined.   This is 
	recommended for long read lengths, e.g., MiSeq with 250bp reads. 
-b/--buffer	Number of bases up and down of the end of the first read to look for the 
	gap between reads with respect to reference.  Only required if --ref specified.  
	Default = 20.  This should be at least twice the expected gap size.  
-w/--clustalw	Use ClustalW alignment with reference sequence to determine gap size. 
	Default is to use bwa mem.  Using ClustalW is about 3-4x slower.  Requires --ref.
-p/--cpu	Number of processors to use.  Default = 1. Determines the number of threads 
	to use for bwa if --ref is supplied.  
-c/--char	Character to use for concatenation.  Default = N.  
-o/--output	Output file name.  Default is [file1 basename, minus \"_1\"].concat.fastq . 
Example:
concatenate_fastq.pl -o Sample_5.concat.fastq -r Seq12_093009_HA_cds.fa --cpu 14 Sample_5_1.fastq Sample_5_2.trim.fastq 
";

# Change log
# 2013-09-23
# Created script.
# 2013-11-14
# Added option --number to put Ns between R1 and R2.  
# 2013-11-15
# Added option --ref to align to reference and add a custom number of Ns between the reads during concatenation.  
# 2013-11-18
# Parallelized concatentation. Optimal cpu number is about 12-14.  ~17x speed-up with 14 cpus.  Processes ~800,000 pairs of reads in 3 hours about 4400 pairs/minute.
# 2013-11-27
# Added option --char/-c.  Also changed short option for --cpu to -p instead of -c.  
# Removed some alternate code in sub get_custom_number_Ns that was commented out (running clustalw manually instead of using Bioperl Run module)
# 2013-12-03
# Added select STDOUT; line before $| = 1;  to try to fix the problem of lines in the output being intermingled.  New problem is that there are many blank lines in the output now. 
# 2013-12-04
# Removed select STDOUT as it was printing blank lines.  Instead, open the output file in append mode in the child processes, lock the file handle, write to the file, then close it.
# Added -o/--output option to give output file a name
# Added the function to align ungapped concatenated reads to the reference sequence with bwa mem to look for gap size to increase speed.  Made this the default mode.   Clustalw is optional, via -w/--clustalw.  Renamed sub get_custom_number_Ns to get_custom_number_Ns_clustalw
# Commented out the Parallel Loop so the BWA method could be fast (3850 pairs/second).  Copying a large hash to all child processes was slowing it down to 10 minutes per 1000 pairs.  
# 2013-12-17
# Changed the base quality to 3 '$' instead of 2 '#' for the inserted Ns.  This way, I can use fastq_quality_filter -q 3 -p 100 to filter out any concatenated reads where either R1 or R2 had an N in the original read if I want to.  
# 2014-01-09
# Removed the Parallel Loop lines that were commented out.
# 2015-01-16
# Changed Processed read counter to every 100K instead of every 10K. 
# 2015-03-16
# Added table of gap lengths
# 2015-03-18
# Added basename to temporary output files to avoid collisions in case no output directory is specified

# To do
# Make a distribution graph of insert sizes when using --ref. 
# Test MAFFT to see if we can speed it up even further compared to Clustalw 
# Make clustalw-specific modules (e.g., Bioperl Run) optional using eval.
# Put the Parallel Loop for clustalw in a separate subroutine.  Right now it's not running in parallel. 
# Maybe make an option to remove any sequences that don't align to the reference...
# Add a script before this one the check whether reads of a pair overlap or not; if so, send to pandaseq, if not, run in concatenate_fastq.pl
# Check first whether bwa index exists before running bwa index
# Maybe change default output directory to something other than Cwd?

unless ($ARGV[1]){	print STDERR "$usage\n";	exit;	}
my $start_time = time;

unless (-d $save){	mkdir($save) or warn "$!\n"	}

if ($clustalw){
	unless ($ref){
		die "Please enter reference sequence for alignment\n";
	}
}
my $refseq = "";	# reference sequence itself
if ($ref && $clustalw){
	my $in = Bio::SeqIO->new(-file => $ref , -format => 'Fasta' );		# There should only be one sequence in this file, so it should be the first.   # This conflicts with the Bioperl module
	my $seq = $in->next_seq();
	$refseq = $seq->seq();
					if ($verbose){	print "refseq: $refseq\n";	}
}

# Run BWA to get gap sizes if reference sequence file is supplied.
my $id_to_gap_hash;		# hashref with keys as fastq id and value is the gap size to insert between R1 and R2.  Determined from CIGAR string in alignment with bwa mem.
if ($ref){
	unless($clustalw){
		if (check_for_bwa() ){		# check_for_bwa() sub will return 1 if bwa mem is not found; will also print an appropriate error message.
			exit 1;
		}	
		else {
			print STDERR "Found bwa mem.\n";
			$id_to_gap_hash = get_gaps_with_bwa($ref, $ARGV[0], $ARGV[1], $buffer);
		}	
	}
}

						# print Dumper($id_to_gap_hash);	exit;

# Set up output file
unless ($output){
	my ($filename,$dir,$ext) = fileparse($ARGV[0],@SUFFIXES);
	$filename =~ s/_\d$/\.concat/;	# e.g.,  Sample_5_1.fastq changed to Sample_5.concat (.fastq added next line)
	$output = $save . "/" . $filename . $ext;  # Don't open here, open in each child fork
}
	
if (-e $output){
	print STDERR "File already exists: $output -- will overwrite\n";
	# sleep 5; 
	unlink($output);
}



open (my $FH1 , '<' , $ARGV[0]) || die $!;
open (my $FH2 , '<' , $ARGV[1]) || die $!;

my $fastq_iterator = iterator($FH1);
my $fastq_iterator_2 = iterator($FH2);


my $count = 0;		# Keep a count of progress.  Need to increment this in the first sub.  

if ($ref){
	print STDERR "Walking through reads to insert gaps between reads. ";
	&elapsed($start_time, 'Elapsed', $verbose);
}

my $entry 	= "";
my $entry2	= "";
my $writefh = open_to_write($output); 
# The BWA method is very fast if you don't use a Parallel Loop.  
while( $entry = $fastq_iterator->()){
	$entry2 = $fastq_iterator_2->(); 
	$count++; 
	
	if ($count % 10000 == 0){
		print STDERR "Processed $count reads. "; 
		&elapsed($start_time, 'Elapsed', $verbose);
	}
			
#	print "record $count\nread1:\n$entry\nread2:\n$entry2\n";
	my @fastq = split(/\n/, $entry);
	my @fastq2 = split(/\n/, $entry2);		# About 225,000 records per second before doing this step.  # Adding this conversion to an array brings speed to about 130,000 records per second
	
	# Check that the ids match.  
	my @ids = split(/\s+/, $fastq[0]);		# e.g., split this:  @MISEQ:50:000000000-A4142:1:1101:17205:1539 1:N:0:ACCACCTT
	my @ids2 = split(/\s+/, $fastq2[0]);
# 	die "Read ids do not match! Exiting...\n" unless ($ids[0] eq $ids2[0]);		# Adding this check to all reads brings speed to about 100,000 records per second.  Just checking every 1000th read to save time (to see if the files became out of sync at a certain point), it runs at 110,000 records per second.  Just check all reads in that case -- it doesn't take much more time.

	my $first 	= Bio::Seq->new( -seq => $fastq[1], 	-id => $fastq[0], );
	my $second 	= Bio::Seq->new( -seq => $fastq2[1], 	-id => $fastq[0], );
	$second 	= $second->revcom();		# Seq object, Reverse complement 
	
	my $first_quals		= $fastq[3];
	my $second_quals 	= reverse($fastq2[3]);  # Reverse quality scores.  
	
	if ($ref){
		if ($clustalw){
			$number = get_custom_number_Ns_clustalw($refseq, $first->seq() , $second->seq(), $buffer);
		}
		else {
			my $id = $ids[0];		# Need just the first part of the id, because that is how the hash keys were saved, e.g., @MISEQ:50:000000000-A4142:1:1101:17205:1539 . Also remove the @ 
			$id =~ s/^@//;
			$number = $id_to_gap_hash->{$id}; 
						# print "id: $id number: $number\n\n-->$id_to_gap_hash->{$id}\n"; exit;
		}
	}

	my $new_fastq = $first->id() . "\n" . $first->seq() . "$char"x$number . $second->seq() . "\n" . '+' . "\n" . $first_quals . '$'x$number .  $second_quals . "\n";		# Use '$' as low quality for gap between R1 and R2, (quality 3.  Was '#', quality 2.)

	print $writefh "$new_fastq";
	
	
}

close($FH1);
close($FH2);
close($writefh);

&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub iterator{
	my $file_handle = shift // die $!;
	return sub{
		my $record .= readline($file_handle) // return;
		$record .= readline($file_handle);
		$record .= readline($file_handle);
		$record .= readline($file_handle);
		return $record;
	};
}
#-----------------------------------------------------------------------------
sub get_custom_number_Ns_clustalw {
	# Arguments
		# Reference sequence, First sequence, second sequence, buffer
	# Concatenates reads with no Ns, makes a temporary fasta file, aligns concatenated read to reference with clustalw, then gets the region that should correspond to the gap between the reads and looks for the largest gap.   
		
	my ($refseq,$first,$second,$buffer) = @_;
	
	my $read_length = length($first);
	
	my $concat = $first . $second;
	
					if ($verbose){	print STDERR "second in sub: $second\n";	}
					if ($verbose){	print STDERR "concat: $concat\n";			}
	
	# Create temporary fasta file
	my $cwd = Cwd::cwd();
	my $tempdir = File::Temp->newdir(  "$cwd/fasta_file_tempXXXXX" );
	my $temp_fasta = $tempdir."/temp.fa";
	my $fh = open_to_write($temp_fasta, 0, 0, 1);
	print $fh ">ref\n$refseq\n>read\n$concat\n";
	close($fh);
	
	# Align with clustalw 
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-OUTORDER => "INPUT", -QUIET => 1, -QUICKTREE => 1, );		# -OUTORDER => "INPUT", -QUIET => 1,
	my $aln = $factory->align("$temp_fasta"); 	# $aln is a SimpleAlign object.  http://search.cpan.org/~cjfields/BioPerl-1.6.901/Bio/SimpleAlign.pm.  
		
	
	# Get a slice of the alignment near the gap
	my $primerID_length = 10;
	my $col = $aln->column_from_residue_number("read", $primerID_length) - $primerID_length;		 # Look for the position of the 10th read residue within the alignment and then subtract 10.  Just in case some of the first bases align to an early region of the reference, especially if the first 10 bases are a random primerID, so they could align to anything.
						if ($verbose){	print STDERR "col: $col\n";	}
	my $slice = $aln->slice($col + $read_length - $buffer, $col + $read_length + $buffer + 10);		# If buffer is 20, then go 20 bp upstream and 30 bp downstream of the end of the R1 read.  The extra 10 is based on the expected gap size.  Maybe it should be command-line option....
						if ($verbose){	
							my $seqIO = Bio::AlignIO->new(-fh => \*STDERR, -format => 'clustalw');	# I had this outside of the "if ($verbose)" block and it slowed things down a LOT.  
							$seqIO->write_aln($slice);	
						}

	# Get the aligned read from the alignment
	my @seqs = $slice->each_seq();
	my $aln_read_seq = $seqs[-1]->seq();		# Works.  Needs to have reference first and read second.  
					if ($verbose){	print STDERR "read: $aln_read_seq\n";	}

	# Look for largest gap.  Gap could be '-' or '.' character.
	my $gap_length = 0;
	while($aln_read_seq =~ m/([\-\.]+)/g){		# Looking for consecutive '-' or '.' characters.
		if (length($1)>$gap_length){
			$gap_length = length($1);				
		}
	}
					if ($verbose){	print STDERR "gap length: $gap_length\n";	}
	return $gap_length;	
}
#-----------------------------------------------------------------------------
 sub check_for_bwa {
	# Check for bwa mem on the PATH.
	# Returns 0 if successful, 1 if an error occurs.
     my $pwd = pwd_for_hpc();
	my $output=`$pwd/bwa mem 2>&1`;  # Captures stdout and stderr.  
	$output = trim($output);	
	
	# Test on command-line: perl -e'use aomisc; my $output=`bwa mem 2>&1`; $output = trim($output); print "-$output-\n"'
	
	# If bwa mem is present, value will start: "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]".  
	# If an older version of bwa is present, value will be: "[main] unrecognized command 'mem'"
	# If no bwa is present, the value will be undef.
	
	if ($output =~ m/Usage: bwa mem/){
		return 0;	#successfully found bwa mem
	}
	elsif($output =~ m/unrecognized command/){
		print STDERR "Please update bwa to version 0.7.5a-r405 or newer.\n";
		return 1;
	}
	else{
		print STDERR "No bwa installation detected.  Please install bwa (version 0.7.5a-r405 or newer) or add to PATH environment variable to run this program.\n"; 
		return 1;
	}
	
}
#-----------------------------------------------------------------------------
sub get_gaps_with_bwa{
	my ($ref, $fastq_file1, $fastq_file2, $buffer) = @_;
#	print join "\t", ($ref, $fastq_file1, $fastq_file2, $buffer),"\n"; exit;
	my $id_to_gap_hash; 	# Hash to store 
	
	
	#### First, go through the input fastq to make fasta files with the sequences concatenated with no Ns.   ###
		# This section borrows a lot of code from above.  It might be nice to put it all in a subroutine at some point to make everything more concise.
	
	# Get file base name to add to temp files.
	my ($base,$path,$ext) = fileparse($fastq_file1,@SUFFIXES);
	$base =~ s/_R?1//g;

	# Create temporary fastq file
	print STDERR "Making temporary fastq file by concatenating with no gap.";  &elapsed($start_time, ' Elapsed', $verbose);
#	my $cwd = Cwd::cwd();

#	my $tempdir = File::Temp->newdir(  "$cwd/concat_fastq_tempXXXXXXX" );
	my $tempdir = $save;
	my $temp_fastq = "$tempdir/$base" . ".temp.fastq";
	my $fh = open_to_write($temp_fastq, 0, 0, 1);

	open (my $first , '<' , $fastq_file1) || die $!;
	open (my $second, '<' , $fastq_file2) || die $!;
	my $fastq_iterator = iterator($first);
	my $fastq_iterator_2 = iterator($second);
	
	my $read_length = 0; 
	while (my $entry = $fastq_iterator->() ){
		my $entry2 = $fastq_iterator_2->();		# About 225,000 records per second
		my @fastq = split(/\n/, $entry);
		my @fastq2 = split(/\n/, $entry2);		# Adding this conversion to an array brings speed to about 130,000 records per second
		my @ids = split(/\s+/, $fastq[0]);
		my @ids2 = split(/\s+/, $fastq2[0]);
		die "Read ids do not match! Exiting...\n" unless ($ids[0] eq $ids2[0]);		# Adding this check to all reads brings speed to about 100,000 records per second.  Just checking every 1000th read to save time (to see if the files became out of sync at a certain point), it runs at 110,000 records per second.  Just check all reads in that case -- it doesn't take much more time.
		$read_length = length($fastq[1]) unless ($read_length);
		my $first 	= Bio::Seq->new( -seq => $fastq[1], 	-id => $fastq[0], );
		my $second 	= Bio::Seq->new( -seq => $fastq2[1], 	-id => $fastq[0], );
		$second 	= $second->revcom();		# Seq object, Reverse complement 
		my $first_quals		= $fastq[3];
		my $second_quals 	= reverse($fastq2[3]);  # Reverse quality scores.  		
		my $number = 0;	
		my $new_fastq = $first->id() . "\n" . $first->seq() . "$char"x$number . $second->seq() . "\n" . '+' . "\n" . $first_quals . "#"x$number .  $second_quals . "\n";		# Use '#' as low quality for gap between R1 and R2, 2
		print $fh "$new_fastq";
	}
	close($first);
	close($second);
	close($fh);
	
#						system("wc $temp_fastq");
	
	####  Now align to the reference sequence with bwa mem  ###
	print STDERR "Aligning temporary fastq file to reference with bwa mem.";  &elapsed($start_time, ' Elapsed', $verbose);
	my $refbase = basename($ref);
	my $temp_sam = "$tempdir/$base" . ".temp.sam";
#	my $temp_sam = "temp.sam";
	my $pwd = pwd_for_hpc();
	my $bwa_index_out 	= "$tempdir/$base" . ".bwa_index_out.txt";
	my $bwa_mem_out 	= "$tempdir/$base" . ".bwa_mem_out.txt";
	my $samtools_view_out 	= "$tempdir/$base" . ".samtools_view_out";

	# [Need to check first whether the bam index exists before running bwa index...]	
	# [Also check whether $ref exists in $tempdir before copying...]

#	my $cmd = "cp $ref $tempdir/; $pwd/bwa index $tempdir/$refbase 2> $bwa_index_out && $pwd/bwa mem -t $cpu -M $tempdir/$refbase $temp_fastq 2> $bwa_mem_out | $pwd/samtools view -S -F 256 - 2> $samtools_view_out > $temp_sam";
	my $cmd = "$pwd/bwa mem -t $cpu -M $tempdir/$refbase $temp_fastq 2> $bwa_mem_out | $pwd/samtools view -S -F 256 - 2> $samtools_view_out > $temp_sam";
#	print $cmd; exit;
	print STDERR "Executing: $cmd\n";
	system($cmd);  # This will copy the reference sequence to the temporary directory, index it, align the reads, filter to get only primary alignments and output a SAM file.
					#	system("cat $tempdir/bwa_mem_out.txt");
					#	system("wc $temp_sam");
	
	#### Now parse the sam file to get id and cigar string and save gaps between reads
	print STDERR "Alignment done, parsing output.";  &elapsed($start_time, ' Elapsed', $verbose);
	my $samfh = open_to_read($temp_sam);
	my $unmapped = 0;
	my $mapped = 0; 
	my $gap_sizes; 	#Hashref of gap size counts
	while(<$samfh>){ 
		chomp; 
		my @F = split(/\t/); 
		my ($id,$cigar_string) = ($F[0],$F[5]);
		unless($id && $cigar_string =~ m/\w+/){		# Check to see that the reads were mapped.  If not mapped, cigar string will be *.  If not mapped, assign gap as zero.
			$id_to_gap_hash->{$id} = 0;	
			$unmapped++;
			next;
		}
		$mapped++;			
		my @numbers = split(/[A-Z]/, $cigar_string); 
		my @letters = split(/[0-9]+/, $cigar_string); 
#					 print "cigar: $cigar_string\nnumbers: ", join " ", @numbers, "letters:", @letters; print "\n"; 
		shift(@letters); 
		my $cigar; 		# AoA of Cigar elements
		for (my $i=0; $i < @numbers; $i++){
			$cigar->[$i] = [$letters[$i], $numbers[$i]]; 
		} 
					#  print Dumper($cigar); 
		my $sum = 0;
		my $stored = 0;
		CIGAR: for (my $i=0; $i < @$cigar; $i++){ 
			if ($sum >= ($read_length - $buffer) && $sum <= ($read_length + $buffer * 2)){
					# position of the deletion/insertion in the read with respect to reference is within one buffer upstream and 2 buffer distances downstream of the end of the R1 read
			 	if($cigar->[$i]->[0] eq 'D'){
						#	print "id: $F[0]\tsum: $sum\tmatch: $cigar->[$i]->[1]\n";
					$id_to_gap_hash->{$id} = $cigar->[$i]->[1]; # store the size of the gap.  Positive value for gap.  
					$gap_sizes->{$cigar->[$i]->[1]}++;
					$stored++;
					last CIGAR;		#?  Is there ever insertion and deletion??
				}
				elsif($cigar->[$i]->[0] eq 'I'){
					# somehow store the position to delete the overlap region.  
					# Negative value for insertion.  
					$id_to_gap_hash->{$id} = -1 * $cigar->[$i]->[1]; # store the size of the gap.  Negative value for an insertion.   (Does this work???)
					$stored++;
					last CIGAR;		#? Unless sometimes insertion and deletion??
				}
			}
			$sum+= $cigar->[$i]->[1];
#			unless ($cigar->[$i]->[1] =~ m/^\d+$/){
#				print STDERR "id: $id\n";
#				die Dumper($cigar->[$i]);
#			}
		}
		unless ($stored){
			$id_to_gap_hash->{$id} = 0;		# If there was no D cigar element in the expected location, then there is no gap in the read, so assign a zero gap.  (actually, there might actually be overlap between these reads... I've found that sometimes there is a gap in the reference with these reads)
			$gap_sizes->{0}++;
		}
	}
	close($samfh);
	print STDERR "$mapped\treads mapped to the reference sequence.\n" if ($mapped);
	print STDERR "$unmapped\treads not mapped to the reference sequence.\n" if ($unmapped);
	&elapsed($start_time, ' Elapsed', $verbose);
	
	# Now print out counts for gap sizes.  
	print STDERR "Gap_size\tCount\n";
	for my $gap (sort {$gap_sizes->{$b} <=> $gap_sizes->{$a}} keys %$gap_sizes){
		print STDERR "$gap\t$gap_sizes->{$gap}\n";
	}

	return $id_to_gap_hash;
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
#Useful things

#my $total_lines = 0;
#	$total_lines++;		#somewhere in subroutine, if I run through the file first.
#	my $lines_done = 0;
#	$lines_done++;		#in while loop, increment as each new line is processed
#	if (($lines_done % 25) == 0){	#put this in the while loop too.
#		print "Done processing $lines_done of $total_lines sequences.";
#		&elapsed($start_time, ' Elapsed', $verbose);
#	}
