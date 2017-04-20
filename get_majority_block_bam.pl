#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use aomisc;
use Cwd;

my $PWD  = pwd_for_hpc();

# Use tools packaged along with this script (same directory)
my $SORTSAM = $PWD . "/picard.jar SortSam";		
my $BWA=$PWD."/bwa";
my $INTERSECTBED=$PWD."/intersectBed";
my $BAMTOBED=$PWD."/bamToBed";
my $GET_MAJORITY_START_STOP_PL=$PWD."/get_majority_start_stop.pl";
my $SAMTOOLS=$PWD."/samtools";
my $JAVA='java';	# Use system java 

my $cpu = 8;		# Number of threads to use for BWA. 
my $ref = '';		# eg data_default/Seq12_093009_HA_cds.fa
my $fastq = '';
my $output_dir = Cwd::cwd();	# Default is the current directory

GetOptions(
   'ref|r=s' => \$ref,
   'fastq|i=s' => \$fastq,
   'cpu|p=s' => \$cpu, 
   'output_dir=s' => \$output_dir,
);

my $usage="
OPTIONS:
-r/--ref	Reference sequence.  Required.  There should be a corresponding BWA index
	for this reference. (Run 'bwa index <reference>' prior to running this script.)
-i/--fastq	FASTQ input file.  Required.
--output_dir	Output directory.  Default = current directory.  Will be created if it
	does not exist.
-p/--cpu	Number of threads for BWA.  Default = 8.
";

die "Please provide a reference file (--ref)\n$usage". $! if (($ref eq '') || (! -e $ref));
die "Please provide a FASTQ file as input (--fastq)\n$usage". $! if (($fastq eq '') || (! -e $fastq));

$output_dir =~ s/\/$//;

system ("mkdir $output_dir") unless -e $output_dir;
die "Failed to make output dir $output_dir\n".$! unless -d $output_dir;

#system ("cp $options{ref} $options{fastq} $options{output_dir}");		# commented out by Andrew O. for now.  2016-01-09.
#my ($ref,$fastq) = map{$options{output_dir}.'/'.extract_file_name($_)}($options{ref},$options{fastq});
#die "failed to copy files properly $ref,$fastq\n".$! unless ((-e $ref) && (-e $fastq));

# # Check to see that the reference is indexed.  If not, attempt to index it
# system ("$BWA index $ref") unless -e $ref.'.bwt';
die "No BWA index found for the reference file! Please run 'bwa index $ref'\n".$! unless  -e $ref.'.bwt';

# Get the base name of FASTQ file (minus extension and directory)
my ($name,$dir,$ext) = fileparse($fastq,@SUFFIXES);

# Align with bwa mem
my $bwaopts="-t $cpu -M -B 1";
my $tmp_sam = $output_dir . '/' . $name.'.sam.gz';
#my $bwa_cmd = "$BWA mem $bwaopts $ref $fastq | gzip > $tmp_sam";
my $bwa_cmd = "$BWA mem $bwaopts $ref $fastq | gzip > $tmp_sam";
print_and_execute($bwa_cmd);

# Sort and index the alignment bam file
my $tmp_bam = $output_dir . '/' . $name.'.bam';
my $sort_cmd = "$JAVA -Xmx3G -jar $SORTSAM I=$tmp_sam O=$tmp_bam CREATE_INDEX=true SO=coordinate";
print_and_execute($sort_cmd);

# Get the majority start stop bed file (single bed region)
my $tmp_bed = $output_dir . '/' . $name.'.majority.bed';
my $majority_bed_cmd = "$BAMTOBED -cigar -i $tmp_bam | awk '\$NF !~ /H|S/{print}' | $GET_MAJORITY_START_STOP_PL > $tmp_bed";
print_and_execute($majority_bed_cmd);

# Intersect bam file with majority bed region (-f 1 -r so that all read alignments start and end at the majority positions)
my $finalbam = $output_dir . '/' . $name.'.majority.bam';
my $intersect_cmd = "$INTERSECTBED -abam $tmp_bam -b $tmp_bed -f 1 -r -u > $finalbam";
print_and_execute($intersect_cmd);

#-----------------------------------------------------------------------------
#----------------------------------- SUBS ------------------------------------
#-----------------------------------------------------------------------------
sub print_and_execute {
	my $cmd = shift;
	print STDERR "Executing command: $cmd\n";
	system($cmd);
}

#sub extract_file_name{
#    my $file = shift;
#    my ($filename,$dir,$ext) = fileparse($file,@SUFFIXES)
#    return $filename . $ext;
#}
