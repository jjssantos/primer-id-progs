#!/bin/bash

# Script to align reads to a reference with bwa mem, convert to BAM and filter the reads so they all have the same start and stop coordinates.
# Requirements:
#   bwa mem
#   BEDTools (intersectBed bamToBed)
#   picard (SortSam)
# 	get_majority_start_stop.pl (custom script)
# Designed for amplicon sequencing, hence all of the output reads will be in a block
# 
# Change log
# 2015-01-16
# Add filter for step to get the majority start/stop to ignore reads with H or S in the cigar string.  



SORTSAMJAR="/usr/local/bio_apps/picard-tools-1.75/SortSam.jar"		#****Required***  Or you can specify with -s
p=10						# Number of threads to use for BWA. 

if [ "$1" == "-p" ]
then
	shift
	p=$1
	shift
fi

if [ "$1" == "-s" ]
then
	shift
	SORTSAMJAR=$1
	shift
fi


ref=$1		#Reference fasta file for 
shift
fastq=$1	#Input file
shift

if [ -z "$ref" ]
then
	echo "Usage: `basename $0` [-p threads] [-s /path/to/SortSam.jar] <reference> <fastq>" >&2
	echo; echo "fastq file has forward and reverse reads concatenated with concatenate_fastq.pl (if no overlap) or pandaseq (if reads are overlapping)" >&2
	echo "Default value for -p is 10" >&2
	echo "Default value for -s is /usr/local/bio_apps/picard-tools-1.75/SortSam.jar" >&2
	echo "reference is the full path to the reference sequence to align to (should be indexed with bwa index; will attempt to index it if it is not already indexed)" >&2
	echo "Example: `basename $0` Seq12_093009_HA_cds.fa Sample_9.contigs.pid10bp.btrim.fastq.2" >&2
	exit 1
fi

# Check to see that the reference is indexed.  If not, attempt to index it
if [ ! -s ${ref}.bwt ]
then 
	echo "No BWA index found for reference, e.g., ${ref}.bwt . Will attempt to index using default settings for bwa index." >&2
	bwa index $ref
	if [ ! -s ${ref}.bwt ]
	then 
		echo "Unable to index the fasta sequence.  Please index the reference before running." >&2
		exit 1
	fi
fi

# Check to see if SortSam.jar is present
if [ ! -s $SORTSAMJAR ] 
then 
	echo "Picard SortSam.jar not found: $SORTSAMJAR.  Please specify location in -s option." >&2
	exit 1
fi

# Start
starttime=$(date +%s)
echo -e "Fastq file: $fastq\nreference file: $ref"

# Initialize temporary directory and set variable file names
#TMPDIR=$(mktemp -d /tmp/output.XXXXXXXXXX)
TMPDIR=temp${RANDOM}; mkdir $TMPDIR
sam=${fastq/.fastq/}.sam.gz
bam=${sam/.sam.gz/}.bam
bed=temp_majority.bed
finalbam=${bam/.bam/}.majority.bam

# Align with bwa mem
bwaopts="-t $p -M"
echo "Running bwa with options $bwaopts"
bwa mem $bwaopts $ref $fastq | gzip > $TMPDIR/$sam				# Was like this: bwa mem -t 10 -M $ref $fastq | samtools view -F 256 -hS - | gzip > $sam;   Don't need to remove secondary alignments with samtools view because those will be outside the majority region anyway.  

# Sort and index the alignment bam file
echo; echo "Sorting BAM file"
java -Xmx3G -jar $SORTSAMJAR I=$TMPDIR/$sam O=$TMPDIR/$bam CREATE_INDEX=true SO=coordinate

# Get the majority start stop bed file (single bed region)
echo; echo "Getting majority start, stop coordinates"
bamToBed -cigar -i $TMPDIR/$bam | awk '$NF !~ /H|S/{print}' | get_majority_start_stop.pl > $TMPDIR/$bed
echo "Majority bed region:" >&2
cat $TMPDIR/$bed

# Intersect bam file with majority bed region (-f 1 -r so that all read alignments start and end at the majority positions)
echo; echo "Filtering reads from BAM file for those matching the majority region."
intersectBed -abam $TMPDIR/$bam -b $TMPDIR/$bed -f 1 -r -u > $finalbam

# Calculate number of reads in the unfiltered and filtered bam files.
echo -e "Before filter:\t$(samtools view $TMPDIR/$bam | wc -l)"
echo -e "After filter:\t$(samtools view $finalbam | wc -l)"

# Clean up temporary files
#rm -r $TMPDIR

endtime=$(date +%s)
totaltime=$(($endtime - $starttime))
echo; echo "${totaltime} seconds total."

