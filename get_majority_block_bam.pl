#!/usr/local/bio_apps/perl-5.16.2/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use aomisc;

my $PWD  = pwd_for_hpc();

my $SORTSAMJAR="/usr/local/bio_apps/picard-tools-1.75/SortSam.jar";
#my $BWA=$PWD."/specific_progs/bwa";
my $BWA='/usr/local/bio_apps/bwa/bwa';
my $INTERSECTBED=$PWD."/intersectBed";
my $BAMTOBED=$PWD."/bamToBed";
my $GET_MAJORITY_START_STOP_PL=$PWD."/get_majority_start_stop.pl";
my $SAMTOOLS=$PWD."/samtools";
my $JAVA="/usr/local/bio_apps/java/bin/java";
my $p=8;						# Number of threads to use for BWA. 

my %options;
$options{ref} = '';		# eg data_default/Seq12_093009_HA_cds.fa
$options{fastq} = '';

GetOptions(\%options,
           'ref=s',
	   'fastq=s',
	   'output_dir=s'
	   );

die "I need a ref file to index \n". $! if (($options{ref} eq '') || (! -e $options{ref}));
die "I need a fastq file to index \n". $! if (($options{fastq} eq '') || (! -e $options{fastq}));
die "I need an output_dir \n". $! if ($options{output_dir} eq '');
$options{output_dir} =~ s/\/$//;

system ("mkdir $options{output_dir}") unless -e $options{output_dir};
die "failed to make output dir $options{output_dir}\n".$! unless -d $options{output_dir};

system ("cp $options{ref} $options{fastq} $options{output_dir}");
my ($ref,$fastq) = map{$options{output_dir}.'/'.extract_file_name($_)}($options{ref},$options{fastq});
die "failed to copy files properly $ref,$fastq\n".$! unless ((-e $ref) && (-e $fastq));

# # Check to see that the reference is indexed.  If not, attempt to index it
# system ("$BWA index $ref") unless -e $ref.'.bwt';
# die "I tried to index $ref but failed.".$! unless  -e $ref.'.bwt';

# Check to see if SortSam.jar is present
die unless -e $SORTSAMJAR;

# Align with bwa mem
my $bwaopts="-t $p -M";
my $tmp_sam = $fastq.'.sam.gz';
system "$BWA mem $bwaopts $ref $fastq | gzip > $tmp_sam";

# Sort and index the alignment bam file
my $tmp_bam = $fastq.'.bam';
system "$JAVA -Xmx3G -jar $SORTSAMJAR I=$tmp_sam O=$tmp_bam CREATE_INDEX=true SO=coordinate\n";

# Get the majority start stop bed file (single bed region)
my $tmp_bed=$fastq.'.temp_majority.bed';
system "$BAMTOBED -cigar -i $tmp_bam | awk '\$NF !~ /H|S/{print}' | $GET_MAJORITY_START_STOP_PL > $tmp_bed";

# Intersect bam file with majority bed region (-f 1 -r so that all read alignments start and end at the majority positions)
my $finalbam= $fastq.'.majority.bam';
$finalbam =~ s/\.fastq//;

system("$INTERSECTBED -abam $tmp_bam -b $tmp_bed -f 1 -r -u > $finalbam");

sub extract_file_name{
    my $file = shift;
    my @a = split /\//,$file;
    return $a[$#a];
}
