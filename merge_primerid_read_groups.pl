#!/usr/local/bio_apps/perl-5.16.2/bin/perl
use warnings;
select STDOUT;		 # Turn off buffering for STDOUT for printing from multiple processes at the same time.
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
use Bio::DB::Sam;		#http://search.cpan.org/~lds/Bio-SamTools-1.37/lib/Bio/DB/Sam.pm
use Bio::SeqIO;
use Bio::Tools::Run::Alignment::Clustalw;
#use Bio::Tools::Run::Alignment::MAFFT;
use Bio::AlignIO;
use Bio::SimpleAlign;
use lib dirname (__FILE__);
use Statistics::Distributions;
use Parallel::Loops;
#use Loops;
use File::Temp; 
use Statistics::R;
use File::Copy;

#use Bio::Perl;
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
# pretty grim, but I can't find a better way of doing this right now
my $prog_loc = Cwd::abs_path($0);	  # philip macmenamin
my @a = split /\//,$prog_loc;	  # philip macmenamin
my $PWD = join '/', @a[0..$#a-1]; # philip macmenamin
my $bam2fastx_bin = $PWD.'/bam2fastx'; # philip macmenamin
my $mafft_bin = $PWD.'/mafft';
my $save;
my $files;
my $verbose;
my $output;
my $gzip;
my $baseq = 0;
my $mapq = 0;
my $fasta;
my $debug;
my $min_reads = 3;
my $max_reads;
my $cpu = 1;
my $gap = 'auto';	# Either auto or a number.  If 'auto', will be replaced by a number later for calculations.
my $R1_length = 'auto';  # Either auto or a number.  If 'auto', will be replaced by a number later for calculations.
my $ambig = 0;
my $wide_gap = 0;
my $clustalw;
my $plot_counts;
my $tiebreaker;
my $ref;
GetOptions('save=s' => \$save, 'output=s' => \$output, 'verbose' => \$verbose, 'files=s' => \$files, 'gzip' => \$gzip, 'baseq=s' => \$baseq, 'mapq=s' => \$mapq, 'fasta=s' => \$fasta, 'debug' => \$debug, 'min_reads=s' => \$min_reads, 'm=s' => \$min_reads, 'max_reads=s' => \$max_reads, 'x=s' => \$max_reads, 'cpu=s' => \$cpu, 'p=s' => \$cpu, 'gap=s' => \$gap, 'g=s' => \$gap, 'R1_length=s' => \$R1_length, 'r=s' => \$R1_length, 'ambig=s' => \$ambig, 'n=s' => \$ambig, 'wide_gap' => \$wide_gap, 'clustalw' => \$clustalw, 'plot_counts' => \$plot_counts, 'tiebreaker' => \$tiebreaker, 't' => \$tiebreaker, 'ref=s' => \$ref);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = basename($0)." takes one or more BAM files and outputs a consensus for each 
primerID or barcode (as a BAM file). The primerID or barcode should be in the read name 
(filtering with filter_fastq_by_barcode_length.pl before alignment will do this for you). 
Each primerID/barcode represents a single original library fragment; groups of reads with 
the same primerID/barcode are probably PCR duplicates; a consensus will be called. Each 
group must have at least 2 sequences (this minimum can be increased; default is actually 3); 
reads in a primerID group are aligned with mafft or clustalw and a consensus is decided 
based on majority rule and ties are converted to ambiguous characters.  Requires mafft 
executable on the PATH (or clustalw if --clustalw option is used).  

Prints out a fasta file of consensus reads *.cons.fasta (and a separate fasta file 
*.ambig.fasta for sequences with ambiguous character counts above the threshold set in 
--ambig), one output read per primerID group.  

Processes about 13 primerID groups/second with -p 20.

OPTIONS:
--files 	One or more BAM Files, comma-delimited.  Required.  If multiple samples 
		are input, they will be treated separately.	(The files can also be received as the 
		first argument.)  Alternatively, fasta files can be accepted as input as well (use 
		.fasta or .fa extension).
--save		Directory in which to save files. Default = pwd.  If folder doesn't exist, 
		it will be created.
-m/--min_reads	Minimum number of reads per primerID group.  Default = 3.
-x/--max_reads	Maximum number of reads per primerID group.  No limit by default.
-p/--cpu	Number of processors to use. Default = 1.
-n/--ambig	Number of ambiguous bases in the consensus read (outside of a middle gap) 
		allowed for reporting consensus reads.  Default = 0.  
-t/--tiebreaker	Use intra-sample consensus as a tie-breaker in cases where there is a tie 
		for the majority base at a position in a primerID group.  e.g., 2G, 2A.  
		Recommended if using --min_reads 2. It is recommended to supply reference fasta to 
		--ref option to speed up the step to create intra-sample consensus.  Requires 
		samtools, bcftools, vcfutils.pl and seqtk on the PATH.  
--ref		Reference cds file that the reads in the BAM file were aligned to.  Used with
		--tiebreaker.
--clustalw	Use clustalw to run alignments. Default is to use mafft.  *Executables must be 
		on the PATH.*  Clustalw is about 1.4x faster but mafft seems to do a better job.
-g/--gap	Approximate gap size between the reads.  Default = 'auto', which means the gap
		size will be auto-detected.  The effect of specifying the gap is that Ns are 
		allowed in this region.  Assumes that the reads were not trimmed before 
		concatenating with Ns.  'auto' is recommended for multiple file input as the value 
		will reset for each file.
-r/--R1_length	Average length of R1 portion of the reads in the concatenated reads.  This 
		parameter is used to define the beginning of an internal gap (indicated as stretch 
		of Ns) between two concatenated reads.  Default = 'auto', which means this will
		be auto-detected from the reads.  If --R1_length is specified, --gap will also be
		turned on.  'auto' is recommended for multiple file input as the value 
		will reset for each file.
--wide_gap	When selecting --gap and --R1_length in 'auto' mode, this setting affects 
		whether the widest possible gap is chosen or a more conservative gap is chosen.  
		If this option is selected, the widest possible gap is chosen (gap start position 
		is smallest R1_length and end position is largest R1_length plus largest gap).  
		The default is a more conservative gap with gap start position as the most common 
		R1_length and gap end position as common R1_length plus most common gap size. 
--plot_counts	Make a graph of primerID group counts based on size of group (i.e., number
		of reads with the same primerID).
";

# To Do:

# Add check for samtools on the command line.   
# Add option for maximum per group.  Groups with too many reads might have a repeated barcode.  Maybe this could be chosen automatically by the distribution of the numbers of reads per barcode group (e.g., don't consider the top 5th percentile)  Or maybe have some way to detect barcode groups that are too different/inhomogeneous to automatically reject...
# Make output BAM file conform well to specifications.  Is there any more information from the original BAM files I could include in or recalculate for the output BAM file?
# Add back ability to report average base quality for reads (either root mean square, fisher's method or some other method) "Base quality scores from multiple reads are merged using Fisher's method."
# Make Bio::Tools::Run::Alignment::Clustalw requirement conditional on --clustalw option (using eval)
# Take fasta or fastq input
# Clean out old subs not being used anymore.

# Change log
# 130301
# Added --group.  Minimum number of reads per barcode group.  Default = 2.
# 130403
# Added ability to output a consensus bam file by default.  BAM file that has one read per barcode group so that we can use vprofiler, vphaser, etc.
# 2013-05-24
# Removed the nucleotide and amino acid tally part (that will be in convert_reads_to_amino_acids.pl or potentially another script altogether); now it just outputs a consensus BAM.  As part of this, removed --labels option because that was just related to the output table.
# Removed the part where it first checks if the read overlaps with an exon.  At this stage, it is not important.  As part of this, removed several subroutines and --gff option. 
# Commented out the fasta file input into the BAM object.  It's not being used at the moment since I'm calling the reference sequence with $a->dna, which uses the MD tag.  
#	--fasta		Reference fasta file.  Required.  Used to break ties when reporting the
#			consensus sequence.
# 2013-06-07
# Set $combined_pvalue to 3.162e-26 if calculated $combined_pvalue is 0 from chisq distribution.
# 2013-09-27
# Modified make_consensus sub to change the output file name.  Also, added samtools sort step.
# 2013-09-30
# Modified read_bam sub to increment $read_tally->{$barcode}->{bad} if reference sequence doesn't match for all of the reads of a primerID group (meaning they aligned to different regions).  Then in make_consensus sub, these primerID groups are skipped.
# Added @pos the list of things to check to see if it is the same for all of the reads in a primerID group.  
# 2013-11-14
# Modified make_consensus and read_bam subroutines.  Almost complete rewrite. (saved the old subroutines as .._old)
# 2013-12-10
# Changed --group option to --min_reads and added --max_reads option. 
# Removed check for min_reads in make_consensus subroutine ("if (scalar(@sequences) >= $min_reads){")now since all of those outside of the range set by min_reads and max_reads are deleted from the data structure in the read_bam step. 
# Removed this part from the description: "The program assumes that indels are sequencing errors (e.g., ION Torrent), so insertions to the reference are deleted and deletions are filled with base X; this results in keeping the original reference frame." since indels are considered now.  
# Removing this temporarily, until I add back functionality to get average base quality score of some sort
# Added --gap option.  
# Ties result in ambiguous codes.  If they are outside of the middle gap (if it exists)
# 2013-12-11
# Fixed some bugs regarding checking whether ambiguous bases are found in the internal gap or not
# Added --R1_length option to help with definition of the internal gap
# Added auto detection for --gap and --R1_length parameters.  Set as default.  Uses the first 1000 reads with a stretch of Ns >= 5bp long to determine the widest possible gap region.   
# Changed @ambiguous to save the sequence and then print it out at the end.  
# Added option --skinnier_gap.
# Printing out a table of positions where ambiguous bases were found at the end of the script
# 2013-12-13
# Made mafft the default aligner for making a consensus.  It's about 1.4x slower than clustalw.  Aligning with clustalw is still accessible by using option --clustalw.
# Made the skinny gap the default and changed the wider gap to the option.  
# 2014-01-09 
# Changed the default minimum to 3.  Updated the $usage.
# Added $count to the temporary directory name to further avoid collisions.  Also put the temp directories in /tmp/ instead of Cwd.  
# 2014-04-01
# Changed make_consensus() sub to print out uppercase sequences to be more compatible with fastx_toolkit, which requires upper-case fasta sequences.  
# 2014-06-02
# Changed read_bam to three subroutines: convert_bam_to_fasta, find_gap, and read_fasta.  This way, we can read the file to find the gap and then read the file again to save to a hash.  It is much cleaner than trying to get the gap and save some reads at the same time.  
# Also changed it so that the number of reads to check is not a minimum, but a maximum, so that files with less than 1000 reads can be searched for a gap as well.  In that case, it uses the $fraction of reads with a gap to determine whether to proceed to set $gap and $R1_length or not.  
# Added an output file of counts per PrimerID group size for graphing. Also bridges to R to make a simple barplot of the PrimerID group sizes, saved as PDF.  
# 2014-06-03
# Removed options:
#--mapq		Minimum read mapping quality to consider.  Default = 0.
#--baseq		Minimum base quality to consider for output.  Default = 0. Try 13.
# Added option --plot_counts to make the output file of counts per PrimerID group and the barplot optional, since a similar graph will also be produced by filter_fastq_by_primerid_length.pl, except that the graph by filter_fastq_by_primerid_length.pl will have number with zero as well.  
# Allowed fasta file as input.  (easier for testing/debugging)
# Added variable $ambig_summary_for_fasta_header_def to store information about ambiguities in the fasta header, to explain the reason a read was placed in the ambig fasta file.  
# Added variable $res_num to keep track of residue numbers separately from $pos position in the alignment, in case of gaps inserted that will be removed later.  
# 2014-06-10
# Added subroutine get_sample_consensus and options --tiebreaker and --ref.  Modified make_consensus subroutine appropriately to use the intra-sample consensus as a tie-breaker.  
# 2014-06-11
# Added step to read_fasta to print out fasta files for reads that are above maximum --max_reads value and below minimum --min_reads value.  This is useful for downstream processing to get background frequency rate for unmerged reads.  Added print_out_of_range_fasta sub to do the printing.  
# 2014-07-17
# Fixed a bug in the table of ambiguous positions where I was printing the line returns to STDOUT and the rest of the table to STDERR.  
# 2014-10-14
# Fixed a bug in read_fasta subroutine where the program was appending belowmin and abovemax reads to existing files. Now it checks if the file exists and deletes it first before writing any reads.


unless ($files||$ARGV[0]){	print STDERR "$usage\n";	exit;	}	#fasta and gff are required
unless($files){	
	if (scalar(@ARGV)>1){	#This will allow you to use wildcard to pass files to the script from the command line. e.g, script.pl *.txt, and it will run through each.
		$files = join ",", @ARGV;
	}
	else{
		$files = $ARGV[0];	
	}
}

# Check options
unless ($gap =~ m/^\d+$/ || $gap =~ m/^auto$/i){
	print STDERR "--gap option should either be a number or 'auto'\n";
	exit;
}
unless ($R1_length =~ m/^\d+$/ || $R1_length =~ m/^auto$/i){
	print STDERR "--R1_length option should either be a number or 'auto'\n";
	exit;
}


my $save_dir = $save || Cwd::cwd();
unless (-d $save_dir){	mkdir($save_dir) or warn "$!\n"	}
print STDERR "save directory: $save_dir\n";


my $start_time = time;
my @files = &aomisc::get_files($files);		#If allowing a directory, specify extension of the files in second argument, e.g., my @files = &get_files($files, 'bed');
my @suffixes = (qw(.bed .bed.gz .bed12 .bed12.gz .txt .txt.gz .BED .BED.gz .BED12 .BED12.gz .fasta .fa .FA .FASTA .FAS .fas), @SUFFIXES);	#for fileparse.  Feel free to add more accepted extensions.  @SUFFIXES comes from aomisc.pm.  

my $iupac = iupac_ambiguities();		# Get a hashref of iupac ambiguity codes.  

# Parse the BAM files and output a consensus fasta file for each input file

for (my $i = 0; $i < @files; $i++){
	my $file = $files[$i];

	my $sample_consensus_seq = "";
	if ($tiebreaker){
		$sample_consensus_seq = get_sample_consensus($file, $ref);	# Takes BAM file and optionally a reference, returns sequence with same coordinates as reference.
	}
	my $fasta = convert_bam_to_fasta($file);
	find_gap($fasta);		# Modifies ($gap,$R1_length) global parameters if either is set as 'auto'; otherwise, provides some warnings if user-defined parameters are outside predicted based on auto-detection

	my $read_tally = read_fasta($fasta);

#			print Dumper($read_tally); exit;
	my $consensus_sequences = make_consensus($read_tally, $file, $sample_consensus_seq);	# Need to pass it the file as well, so it can use the name to construct a new name for the consensus bam/fastq/fasta file
	
}




&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub get_sample_consensus {
	# Read a BAM file and create a single fasta consensus sequence for all reads in the file.
	# Requires samtools mpileup, bcftools, vcfutils.pl and seqtk (for BAM input)
	# If a reference fasta is provided (the reference to which the reads were aligned to make the BAM file), this subroutine runs much faster. (e.g., 2 seconds compared to 27 seconds for 150K-read BAM file.)
	# If input file is FASTA format, then align with BWA first
	
	my ($file,$ref) = @_;
	
	# Set up temporary directory 
	my $tempdir = File::Temp->newdir(  "/tmp/fasta_file_sample_consensus_tempXXXXXX" );
	
	# Make sure we have a BAM file for making a consensus.  If not, align to the reference with BWA MEM.
	my $bam = "";	# bam file name  
	if ($file =~ m/.bam$/i){
		$bam = $file;
	}
	else {
		# Input file is FASTA
		# Align the reads to reference first.
		die "Please provide reference file if input file is FASTA format!\n" unless ($ref);
		
		print STDERR "Input is FASTA format -- attempting to align to reference before computing intra-sample consensus.\n\n";
		my ($file_prefix,$dir,$ext) = fileparse($file,@suffixes);
		my $bam_prefix = $tempdir ."/". $file_prefix;
		$bam = $bam_prefix . ".bam";
		my $temp_sam = $tempdir . "/temp.sam.gz"; 
#		my $bwa_mem_cmd = "bwa index $ref; bwa mem -t $cpu -M $ref $file | gzip > $temp_sam; java -Xmx3G -jar $SORTSAMJAR I=$tempdir/temp.sam.gz O=$bam CREATE_INDEX=true SO=coordinate";
		my $bwa_mem_cmd = "bwa index $ref; echo; bwa mem -t $cpu -M $ref $file | samtools view -uS -| samtools sort - $bam_prefix && samtools index $bam";
		print STDERR "Running command:\n$bwa_mem_cmd\n";
		system($bwa_mem_cmd);	 # Run alignment
#		print "Screen Output:$output\n";
		print "\nDone Alignment\n\n";
	}
	
	# Create temporary fasta file of the sequences for this primerID.  
	my $temp_fasta	= $tempdir."/temp.fa";
	
	# Construct the command based on whether a reference fasta is provided or not.  
	# Note that this method can call an ambiguous residue in the consensus if there aren't very many reads
	my $cmd = "";
	if ($ref){
		$cmd = "samtools mpileup -uf $ref $bam | bcftools view -cg - | vcfutils.pl vcf2fq | seqtk seq -A -  > $temp_fasta";
	}
	else {
		$cmd = "samtools mpileup -u $bam | bcftools view -cg - | vcfutils.pl vcf2fq | seqtk seq -A -  > $temp_fasta";	
	}

	# Run the command to create the consensus
	print STDERR "Running command:\n$cmd\n";
	my $output = `$cmd 2>&1`;		# Captures STDERR and STDOUT

	# Read the consensus fasta file.  
	my $in  = Bio::SeqIO->new(-file => "$temp_fasta" ,		
                           -format => 'Fasta');			# http://search.cpan.org/~cjfields/BioPerl-1.6.922/Bio/SeqIO.pm
	my $seqobj = $in->next_seq();
	$in = undef;
	
	# Return consensus sequence
	my $seq = uc( $seqobj->seq() );
	if ($seq){
		# Clean up intra-sample onsensus sequence	
		$seq =~ s/^N+//ig;		# Remove any trailing Ns.  
		$seq =~ s/N+$//ig;
		print "Consensus seq: $seq\n";
		return $seq;
	}
	else {
		print STDERR "Unable to make intrasample consensus. Make sure samtools, bcftools, and seqtk are on your PATH.\nScreen output:\n$output\n";
		exit;
	}
}
#-----------------------------------------------------------------------------
sub convert_bam_to_fasta {
	# This subroutine takes a bam file as input 
	# First step is to convert to fasta using bam2fastx from the tophat package. (or maybe fastq in the future if I want to add quality scores...)
		# Do we still want to have the ability to set a minimum base quality, or should we just ignore this and take the fasta files as they are?  
			# If we want to be able to do this (essentially it converts the base to N so that it will not count towards a majority), then we need to keep quality information, so fastq or BAM.  
			# Otherwise, if we are requiring at least 3 or 5 reads per primerID group then majority rule might be sufficient to correct these things. 
			# But if there is a variant, we might want to look at the average quality of the bases to determine if it is real (like GATK and other variant callers...)
			# Or maybe we validate the accuracy using alternative methods... but how?  
			# Email from Will:
			#	"Indeed, perhaps we should consider the quality score given their data, but you have a better feeling for that. We are already eliminating ambiguous calls – at least in the consensus. (right?)"
	# Also, maybe use the bam file to get an intrapatient consensus to aid in the alignment... (except Will said don't use reference anymore to call consensus...what about ties?  call ambiguous characters)

	# If input file is fasta, no conversion is performed.  
	
	my $bam_file = shift;
	# Convert BAM to Fasta
	# E.g., bam2fastx -a -o Sample_5.concat.20K.primerid10bp.btrim.1.majority.fasta -A Sample_5.concat.20K.primerid10bp.btrim.1.majority.bam 2> /dev/null
	my ($bam_prefix,$dir,$ext) = fileparse($bam_file,@suffixes);
	
	# First check to see whether the input file has been converted to fasta already.  Check based on the extension.  
	#if ($ext =~ m/bam/i){	# < andrew
	if (1){			# < philip macmenamin
		my $fasta = $save_dir . "/" . $bam_prefix . '.fasta'; 
		my $cmd = "$bam2fastx_bin -a -o $fasta -A $bam_file 2> /dev/null";
		print STDERR "Executing command to make fasta file: $cmd\n";
		system($cmd);
		return $fasta;
	}
	elsif($ext =~ m/(fa|fna|fasta|fas)$/i){
		# Already Fasta file
		print STDERR "Input file is FASTA format\n";
		return $bam_file
	}
	else {
		print STDERR "Unrecognized input format: $bam_file\n\tDoes your file have bam or fa/fna/fasta/fas as an extension?\n";
		exit;
	}
}
#-----------------------------------------------------------------------------
sub find_gap {

	my $fasta = shift;
	
	# Now read in the Fasta file and group the sequences by primerID.  
	my $apparent_long_gap_count = 0; 	# Number of reads found with a stretch of Ns >5.  Sometimes the user might forget to set --gap and there are stretches of Ns in the middle.  These will all be counted against the read as ambiguous unless --gap is set.  Help the user out.  
	my $number_to_check = 1000; 	# Maximum number of sequences to read for getting gap size/location information
	my %gap_sizes;	# Save the gap sizes
	my %gap_positions;	# Save the gap positions
	my $total = 0; 		# Total number of reads processed so far.

	my $in  = Bio::SeqIO->new(-file => "$fasta" ,		
                           -format => 'Fasta');			# http://search.cpan.org/~cjfields/BioPerl-1.6.922/Bio/SeqIO.pm
	FASTA: while ( my $seqobj = $in->next_seq() ) {
		# e.g., 
		#	>MISEQ:50:000000000-A4142:1:1101:10000:4300:ACGACTTGCA
		#	ATTCGAAAGATTCAAAATATTTCCCAAAGAAAGCTCATGGCCCGACCACAACACAAACGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCNNNNNNNTAACAGTAAGGAACAACAGGATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
		$total++;
		my $seq = $seqobj->seq();
		
		### Check for an internal stretch of Ns.  Autodetect --gap and --R1_length
		if ($apparent_long_gap_count < $number_to_check){		# Count all reads, or until 1000 reads are found with a gap.  Could do more if desired.  
			if ($seq =~ m/(N+)/i){
				my $length = length($1); 
				if ($length >= 5){		# Long gap set arbitrarily as 5xN.  Could be parameterized...
					$apparent_long_gap_count++;
					$gap_sizes{$length}++;
					my $pos = $-[0];	# Index.  (In 1-based, this would be just before the first N)
					$gap_positions{$pos}++;
#									print STDERR "pos: $pos\tseq: $seq\n";
				}
			}
		}
		last FASTA if ($total >= $number_to_check);
	}
	
	# Now we read enough.  Unset the SeqIO object and process the results of the gap
	$in = undef;
	
	# Process the gap sizes and positions and set/check --gap and --R1_length parameters.
	my $fraction = 0;
	$fraction = $apparent_long_gap_count / $total if ($total); 		# Fraction of the reads with a gap.  
	my $common_size = find_key_with_biggest_value(\%gap_sizes);
	my $common_position = find_key_with_biggest_value(\%gap_positions);
	if ($fraction > 0.2){		# Arbitrarily set to 20%.  It will probably be close to 99% in reality.  
		printf STDERR "Many reads were found with a stretch of Ns >= 5 ($total checked in total, found $apparent_long_gap_count).\nFraction of reads: %.3f\nMost common position: %2d (%.3f of reads with gap)\nMost common size: %2d (%.3f of reads with gap)\n", $fraction, $common_position, $gap_positions{$common_position}/$apparent_long_gap_count, $common_size, $gap_sizes{$common_size}/$apparent_long_gap_count; 		
		my @gap_size_range 	= range([keys %gap_sizes]);
		my @gap_pos_range 	= range([keys %gap_positions]);
#									print join " ", @gap_size_range, "\n";
#									print join " ", @gap_pos_range, "\n";
#									print Dumper(\%gap_positions);
#									print Dumper(\%gap_sizes);
		my $earliest_pos	= $gap_pos_range[0];
		my $biggest_gap		= $gap_size_range[1];
		
		# If the user set 'auto' for $gap, $R1_length, or both, then provide values for those. Otherwise, provide warnings about their parameters as necessary.
		if ($R1_length =~ m/auto/i && $gap =~  m/auto/i){	# Then we'll set both of these parameters for the user
			if ($wide_gap){		# Wider gap
				# If both parameters are using auto, we can set these to produce the widest possible gap
				# Use the earliest position, and the largest stretch of Ns found for the gap size, starting from the latest position.
				$R1_length	= $earliest_pos;
				$gap 		= $biggest_gap + ($gap_pos_range[1] - $gap_pos_range[0]);
			}
			else {				# Default more conservative gap size.  
				$R1_length	= $common_position;
				$gap 		= $common_size;
			}
		}
		elsif($gap =~  m/auto/i){		# User set $R1_length, probably similar to common_position (should we check?).  Then we'll just set gap parameter for the user.  Use something like the common_size
			print STDERR "Setting --gap parameter to $common_size\n";
			$gap = $common_size;  	# Good?
		}
		elsif($R1_length =~ m/auto/i){	# User set $gap but not R1_length.  $gap will probably be a reasonable value.  Set $R1_length for user.
			print STDERR "Setting --R1_length parameter to $common_position\n";
			$R1_length = $common_position;		# Or should it be set based on their gap size?  						
		}
		else {	# User set both parameters.  Maybe suggest better values for the user.
			if ($gap != $common_size){
				print STDERR "WARNING: Most common gap size found was $common_size but $gap was used as input.  You may want to change the gap size or use 'auto'.\n";
			}
			if ($R1_length != $common_position){
				print STDERR "WARNING: Most common R1_length found was $common_position but $R1_length was used as input.  You may want to change the R1_length size or use 'auto'.\n";							
			}
		}
		
		# Additional warning for when user sets --gap parameter that is smaller than auto-detected gap size
		if ($gap < $common_size){
			print STDERR "WARNING: Most common gap size found was $common_size but $gap was used as input.  You may want to increase the gap size or use 'auto'.\n";		
			sleep 5;				
		}

	}
	else {
		# Reads with a stretch of Ns >=5 represent a small fraction of the total reads analyzed so far (< 20%), so they are probably just sequencing errors, not an actual gap from concatenating reads.
		# If the user set some parameters for --gap and --R1_length, then we can tell them that it's not likely that there is a gap actually.  
		# Or if 'auto' was set for either one, we can just set them both to 0.  Don't worry about checking to see if they set a parameter for one but not the other; if either is set to zero, then it's as if both are set to zero.
		if ($R1_length =~ m/auto/i || $gap =~  m/auto/i){		# The user set one or no parameters.
			print STDERR  "No significant gaps >= 5bp were detected in the reads.  $total reads were processed and only $apparent_long_gap_count reads had a gap.\n";
			($gap,$R1_length) = (0,0);
		}
		else{	# The user set both parameters, but we suggest they should reconsider
			if ($gap >= 5){		# See if the user expected gaps >= 5bp long
				 print STDERR  "No significant gaps >= 5bp were detected in the reads.  $total reads were processed and only $apparent_long_gap_count reads had a gap.\n";
			}
		}
	}

	print STDERR "Using $R1_length and $gap as --R1_length and --gap settings, respectively.\n";
#	sleep 5;

}
#-----------------------------------------------------------------------------
sub read_fasta {
	
	# Next, group the reads by primerID (in the header) and then return this as a data structure.  
	# During the process of reading the files, we will also auto-detect the gap sizes and average starting position of the gap
	
	my $fasta = shift;
	
	my $read_tally; 	# HoHoA.  First key is primerID, second key is 'sequence', value is array of sequences for alignment.  

	
	# Now read in the Fasta file and group the sequences by primerID.  
	my $total_reads = 0;
	print STDERR "Grouping sequences by primerID.\n";
	my $in  = Bio::SeqIO->new(-file => "$fasta" ,		
                           -format => 'Fasta');			# http://search.cpan.org/~cjfields/BioPerl-1.6.922/Bio/SeqIO.pm
	while ( my $seqobj = $in->next_seq() ) {
		# e.g., 
		#	>MISEQ:50:000000000-A4142:1:1101:10000:4300:ACGACTTGCA
		#	ATTCGAAAGATTCAAAATATTTCCCAAAGAAAGCTCATGGCCCGACCACAACACAAACGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCNNNNNNNTAACAGTAAGGAACAACAGGATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
		$total_reads++;
		my $seq = $seqobj->seq();
		my $id = $seqobj->display_id();		
		my @id = split(/:/, $id);
		my $primerID = $id[-1];
		warn "Primer ID doesn't look right: $primerID from id: $id\n" unless ($primerID =~ m/^[ACGT]+$/i);	# Should N be allowed?
		if ($primerID =~ m/N/i){
			warn "Primer ID contains an N; skipping sequence\n";
			next;
		}
		push @{$read_tally->{$primerID}->{sequence}}, $seq;
	}	
		
	# Make a plot of counts by primerID group size
	plot_counts($fasta, $read_tally) if ($plot_counts);
	
	## Filter the data structure $read_tally at this point, removing entries for primerIDs that have fewer sequences than $min_reads or more than $max_reads.
	# Adding this step here instead of filtering during the parallel loop in make_consensus seems to make the next step much faster. (about 26x speedup with my test file of 20K reads, albeit that removes most of the reads from the hash...)
	# Also print out the reads that are less than --min_reads and greater than --max_reads to separate files.  "belowmin" "abovemax"
	my $too_few = 0;		# Count for number of primerID groups < $min_reads
	my $too_many = 0;		# Count for number of primerID groups > $max_reads
	my $just_right = 0;		# Count for number of primerID groups within limits set.
	my $total_groups = 0;	# Count for all primerIDs groups. 
	
	# Prepare output files for reads below --min_reads and above --max_reads 
	my ($file_prefix,$dir,$ext) = fileparse($fasta,@suffixes);
	my $belowmin_file 	= $save_dir . "/" . $file_prefix . ".belowmin.fasta";
	my $abovemax_file	= $save_dir . "/" . $file_prefix . ".abovemax.fasta";
	unlink($belowmin_file) if (-e $belowmin_file); 
	unlink($abovemax_file) if (-e $abovemax_file); 
	
	foreach my $primerID (keys %$read_tally){
		if (scalar(@{$read_tally->{$primerID}->{sequence}}) < $min_reads){
			# Save to "belowmin" fasta file
			print_out_of_range_fasta($belowmin_file,$primerID,$read_tally->{$primerID}->{sequence});
			
			# Delete record
			delete($read_tally->{$primerID});
			$too_few++;
		}
		elsif($max_reads && scalar(@{$read_tally->{$primerID}->{sequence}}) > $max_reads){		# Check to see if $max_reads is defined, and if so, check to see whether 
			# Save to "abovemax" fasta file
			print_out_of_range_fasta($abovemax_file,$primerID,$read_tally->{$primerID}->{sequence});
			
			#Delete record
			delete($read_tally->{$primerID});
			$too_many++;
		}
		else {
			$just_right++;
		}
		$total_groups++;
	}
	
	
	printf STDERR "Total primerID groups:\t%2d\n", $total_groups;
	printf STDERR "primerID groups smaller than $min_reads removed:\t\t%i (%.3f)\n", $too_few, $too_few/$total_groups 	if ($too_few);
	printf STDERR "primerID groups larger than $max_reads removed:\t\t%i (%.3f)\n", $too_many, $too_many/$total_groups	if ($too_many);
	printf STDERR "primerID groups retained for analysis:\t\t%i (%.3f)\n", $just_right, $just_right/$total_groups		if ($just_right);
	
	unless($just_right){
		print STDERR "No PrimerID groups retained, exiting.\n";
		&elapsed($start_time, 'Total', $verbose);
		exit;
	}
	
	&elapsed($start_time, 'Elapsed', $verbose);
		
	return $read_tally;
	
}
#-----------------------------------------------------------------------------
sub print_out_of_range_fasta {
	my ($file,$primerID,$seq_array) = @_;
	my $out_of_range_fh		= open_to_write($file,0,">>");		# Open in append mode
	my $j = 0;	# count of seqs printed so far for this primerID group.
	foreach my $seq (@$seq_array){
		$j++;
		print $out_of_range_fh ">$primerID"."-"."$j\n$seq\n";
	}
	close($out_of_range_fh);
}
#-----------------------------------------------------------------------------
sub plot_counts {
	my ($fasta,$read_tally) = @_;
	## Graph Primer ID group size distribution

	# Before filtering the $read_tally data structure, make a graph of the primerID groups by number of reads per group, and then how many groups have that number.  
	# For now just make a simple barplot.  Single sample.  X-axis is "Counts in each PrimerID group" and Y-axis is fraction of Primer IDs
	# It would be nice to make a summary plot with all samples in the workflow as well.  
	# So, I should make a temporary file that can be read into R and parsed later as well.
	# Maybe make this optional in case a similar graph is produced by filter_fastq_by_primerid_length.pl.
	my ($fasta_prefix,$dir,$ext) = fileparse($fasta,@suffixes);
	my $group_count_file = 	$dir . $fasta_prefix . ".group.counts.txt";
	my $group_count_graph = $dir . $fasta_prefix . ".group.counts.barplot.pdf";
	my $group_counts; 	# Hashref to store the counts.  keys, group size; value, number of primerID groups
	foreach my $primer_id_group (keys %$read_tally){
		my $number_in_this_group = scalar(@{$read_tally->{$primer_id_group}->{sequence}});
		$group_counts->{$number_in_this_group}++;
	}
	my $total_primer_id_groups = scalar(keys %$read_tally);
	my $count_fh = open_to_write($group_count_file);
	print $count_fh "Reads_in_PrimerID_group\tNumber_of_groups\tFraction_of_groups\n";
	foreach my $group_size (sort {$a <=> $b} keys %$group_counts){
		printf $count_fh "%2d\t%2d\t%.3f\n", $group_size, $group_counts->{$group_size}, $group_counts->{$group_size} / $total_primer_id_groups;
	}
	close ($count_fh);

	# Now graph this in R
	my $R = Statistics::R->new();		# http://search.cpan.org/~fangly/Statistics-R-0.31/lib/Statistics/R.pm
	my $out1 = $R->run(
		"table <- read.delim(\"$group_count_file\")",
		"pdf(\"$group_count_graph\")",
		"barplot(table\$Fraction_of_groups, names.arg=table\$Reads_in_PrimerID_group, ylim=c(0,1.0), ylab=\"Fraction of PrimerID Groups\", xlab=\"PrimerID Group Size (Number of Reads)\", main=\"PrimerID Group Size Distribution for $fasta_prefix\", cex.main=0.9)",
#		"library(ggplot2)",
#		"ggplot(table, aes(x = Reads_in_PrimerID_group, y = Fraction_of_groups)) + geom_bar(stat = \"identity\")",		# Testing, so that I can add labels above the bars for the actual numbers
		"dev.off()"
	);
	
}
#-----------------------------------------------------------------------------
sub make_consensus {
	# Maybe rename to reflect that we are aligning sequences...
	# Takes read_tally HoHoA, aligns sequences for each barcode and calls a consensus by 
	# looking at the alignment, column by column.
	my ($read_tally,$bam_file,$sample_consensus_seq) = @_;		# $sample_consensus_seq will be empty if --tiebreaker is not specified.  
	
	# Prepare output files
	my ($bam_prefix,$dir,$ext) = fileparse($bam_file,@suffixes);
	my $out_file_good  = $save_dir . "/" . $bam_prefix . '.cons.fasta'; 
	my $out_file_ambig = $save_dir . "/" . $bam_prefix . '.ambig.' . $ambig . '.fasta';

	my @middle = ($R1_length + 1, $R1_length + $gap); 	# Start and stop of the internal gap, basically.  If no gap, then these should both be zero.  These are 1-based coordinates of the gap positions, inclusive.  (check with a file that has no gap...)
	
#	print Dumper(\@middle);exit;
	
	print STDERR "Looking at groups of sequences to call a consensus for each primerID group.\n";
	
	
	my $pl = Parallel::Loops->new($cpu);

	my $count = 0;		
	# Beginning of parallel loop.  See http://search.cpan.org/~pmorch/Parallel-Loops-0.03/lib/Parallel/Loops.pm
	# Syntax:   $pl->while($conditionSub, $childBodySub)
	# Conceptually similar to 
	#	while($conditionSub->()) {
	#		$childBodySub->();
	#	}
	my $primerID 	= "";
	my @primerIDs	= keys %$read_tally; 
#	my $sequences; 	# AoA[oA]
#	foreach my $key (keys %$read_tally){
#		my $sequence;	# arrayref
#		$sequence->[0] = $key;
#		$sequence->[1] = $read_tally->{$key}->{sequence}; 
#		push @$sequences, $sequence;		
#	}
#	print Dumper($sequences); #exit;
	
	# Set up shared variables to store data.
	my (@ambiguous,@good,@ambiguous_positions);	 # arrays to save or count the number of primerid groups that are completely good or that have one or more ambiguous character.  For @good and @ambiguous, each element is an arrayref with primerID and consensus sequence.  @ambiguous_positions will have an array of all ambiguous positions found in the reads; this array can't be used to count the number of ambiguous sequences since some sequences will have multiple ambiguous sequences.
	$pl->share(\@ambiguous, \@good, \@ambiguous_positions);
	
	$pl->while( sub {$count++; $primerID = $primerIDs[$count - 1]; },	sub {			# was $pl->foreach( \@primerIDs,	sub {	
#		$primerID  = $_; 
#		print "id: $primerID\n";
		my @sequences = @{$read_tally->{$primerID}->{sequence}}; 
#								if($verbose){		print Dumper(\@sequences);	}
	
		if ($count % 1000 == 0){
			print STDERR "Processed $count primerID groups. "; 
			&elapsed($start_time, 'Elapsed', $verbose);
		}

		# Create temporary fasta file of the sequences for this primerID.  
		my $cwd = Cwd::cwd();
		my $tempdir = File::Temp->newdir(  "/tmp/fasta_file_".$count."_tempXXXXXX" );
		my $temp_fasta	= $tempdir."/temp.fa";
		my $temp_aln 	= $tempdir."/temp.aln";
		my $fh = open_to_write($temp_fasta, 0, 0, 1);
		my $seq_count = 1;	# sequence number within the primer id group
		foreach my $sequence (@sequences){
			print $fh ">$primerID-$seq_count\n$sequence\n";
			$seq_count++;
		}
		
		# Add reference if needed for tiebreaker
		if ($tiebreaker){
			print $fh ">ref\n$sample_consensus_seq\n";
		}
		
		close($fh);
#								if ($verbose){		system("cat $temp_fasta");		}
		
		# Align with clustalw of mafft
		my $aln;
		if ($clustalw){
			my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-OUTORDER => "INPUT", -QUIET => 1, -QUICKTREE => 1, );		# -OUTORDER => "INPUT", -QUIET => 1, -QUICKTREE => 1
			$aln = $factory->align("$temp_fasta"); 	# $aln is a SimpleAlign object.  http://search.cpan.org/~cjfields/BioPerl-1.6.901/Bio/SimpleAlign.pm.  
		}
		else {

	#		my $factory = Bio::Tools::Run::Alignment::MAFFT->new("-clustalout" => 1, "-quiet" => 1, );		# For some reason, it's ignoring these parameters...
	#		my $aln = $factory->align("$temp_fasta"); 	# $aln is a SimpleAlign object.  http://search.cpan.org/~cjfields/BioPerl-1.6.901/Bio/SimpleAlign.pm.  

			my $cmd = "$mafft_bin --quiet --thread $cpu --op 1.4 $temp_fasta  > $temp_aln";		# Default gap opening penalty 1.53 is a little too stringent sometimes.  e.g., 
#>Seq12_part
#TTCGAAAGATTCAAAATATTT CCC AAA GAA AGC TCA TGG CCC GACCACAACACAACCGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCACCCGCCTAACAGTAAGGAACAACAGAATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
#>GGAAAGACGG
#TTCGAAAGATTCAAAATATTT CCC AAA AA AGC TCA TGG GCCC GACCACAACACAAACGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCNNNNNNNTAACAGTAAGGAACAACAGAATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
			system($cmd);
			my $aln_in = Bio::AlignIO->new(
				-file 	=> $temp_aln,
				-format	=> 'fasta',
			);
			$aln = $aln_in->next_aln();
		}
		
		my $seqIO = Bio::AlignIO->new(-fh => \*STDERR, -format => 'clustalw');
								if ($verbose){	print "ALIGNMENT: \n"; $seqIO->write_aln($aln);	}
		
		# Save the aligned reference and them remove from the alignment if $tiebreaker is set.  See http://search.cpan.org/~cjfields/BioPerl-1.6.923/Bio/Align/AlignI.pm
		my $aligned_ref_seq = "";
		if ($tiebreaker){
#					$seqIO->write_aln($aln);		# Before removing ref
			my @seqs = $aln->each_seq_with_id("ref");	# Array of seq objects containing ref only
			my $aligned_ref_locatable_seq_obj = shift(@seqs);		# Bio::LocatableSeq object with aligned reference sequence
			$aligned_ref_seq = $aligned_ref_locatable_seq_obj->seq();	# aligned reference sequence.  Same positions as the aligned reads.
#					print STDERR "reference: $aligned_ref_seq\n";
			$aln->remove_seq($aligned_ref_locatable_seq_obj);		# alignment of reads only, without reference.  
#					$seqIO->write_aln($aln);		# After removing ref
		}		
		
		# Walk through the alignment column by column
		my $len = $aln->length();
		my @cons_seq; 	# Consensus sequence for this primerID group
		my $ambig_pos; 	# Hashref of ambiguous positions where key is position and value is an AoA of the tied bases and their counts.
		my $nongap_ambig_count = 0; 
		my $ambig_summary_for_fasta_header_def = "";	# String to add to fasta header after the primerID.  positions and ambiguous characters that were found outside of the gap.  Helpful for understanding why a read was put in the ambig fasta file.  
		my $res_num = 1;		# Residue number.  This is the position of the residue within the final consensus sequence.  
		for (my $pos = 1; $pos <= $len; $pos++){		# $pos is the position in the alignment.  It is usually the same as $res_num, residue number, except when there is an indel (-)
			my %count; 
			foreach my $seq ($aln->each_seq() ) {		# Looking at each base in the column, counting up the bases.
				#http://search.cpan.org/~cjfields/BioPerl-1.6.901/Bio/SimpleAlign.pm
				my $res = $seq->subseq($pos, $pos);
				$count{$res}++;
			}
			foreach my $res (keys %count) {
#								if($verbose){		printf "Pos: %2d  Res: %s  Count: %2d\n", $pos, $res, $count{$res};		}
			}		
			my ($max_key_value) = find_key_with_biggest_value(\%count,2,1);	# returns arrayref of results, including ties
#								if($verbose){	if (scalar(keys %count)>1){	my $top = find_key_with_biggest_value(\%count); printf "Pos: %2d  Res: %s  Count: %2d\n", $pos, $top, $count{$top};	}	}
			my $res = $max_key_value->[0]->[0]; 	# If there is a clear winner, this will work.  Note that this can sometimes be an N, which we need to check.  
			if (scalar(@$max_key_value) > 1){
				# Then there is an ambiguous base because it couldn't be decided on a winner
				$ambig_pos->{$res_num}=$max_key_value;		# Store the AoA of residue counts, to keep track of residue position
				#$VAR1 = [
				#          [
				#            'T',
				#            2
				#          ],
				#          [
				#            'C',
				#            2
				#          ]
				#        ];
				if ($tiebreaker){
					# Assign sample consensus residue 
					$res = substr($aligned_ref_seq, $pos - 1, 1);	# Check is it right?  yes, it's getting the right base, even with a gap.  
					if ($verbose){	print STDERR "assign sample cons: $pos, $res\n"; 	}
				}
				else {
					# Assign ambiguous IUPAC residue or N
					my @ambig_res;		# Array of the nucleotide residues that tied for maximum.
					foreach(@$max_key_value){
						push @ambig_res, $_->[0];
					}
					my $ambig_res = join "", sort @ambig_res;		# Join the two or more residues to get IUPAC consensus		
					if ($ambig_res =~ m/[N\.-]/i){	# ambiguities with Ns or gaps are assigned as "N" (?)
									if ($debug){	print STDERR "An ambiguous base with N, '.' or '-'\tPos: $res_num\n";	}
						$res = "N";
					}
					else {
						# Max residues do not include N.  Get IUPAC ambiguous base.  
									if ($debug){	print STDERR "An ambiguous base with $ambig_res\tPos: $res_num\n";	}
						$res = $iupac->{$ambig_res};				
					}
	#				$ambig_summary_for_fasta_header_def .= " " . $res_num . uc($res);	# add the position and resulting ambiguous character to the fasta header.  Do below instead...
				}
			}
			
			unless($res){	
				# Then something's wrong  
				warn "Residue has no value...\n";
				print STDERR Dumper($max_key_value);	
			}
			
			# If the base is ambiguous or N and OUTSIDE of the gap region (between concatenated reads), flag this column as ambiguous.  
			# If --tiebreaker is set, then the ambiguous base above will be overwritten by the sample consensus base, so it shouldn't be flagged as ambiguous.  
			if ($res =~ m/[NMRWSYKVHDB]/i){	# ambiguous bases OR Ns from gap
				unless ( $gap && (($res_num >= $middle[0]) && ($res_num <= $middle[1])) ){ 		# If the "N" residue is found in the middle gap created by concatenating two reads together, then it is not counted against the read as an ambiguous base.  Only check this if $gap is defined.  Otherwise, all Ns are treated as ambiguous bases.  
#									print STDERR "nongap_ambig at position $res_num\n";
					$nongap_ambig_count++;		# The ambiguous residue is found outside of the gap, so counts against the read.
					$ambig_summary_for_fasta_header_def .= " " . $res_num . uc($res);
#							if ($debug){
#								$seqIO->write_aln($aln);					
#							}
				}
					# But what if we find a non-N ambiguous base inside the middle gap region...
					# This could happen if the gap for this pair is a little smaller than expected
					# Should I set a flag/warning for this?  
				if ($gap && $res =~ m/[MRWSYKVHDB]/i && (($res_num >= $middle[0]) && ($res_num <= $middle[1])) ){		# non-N ambiguous character found in the gap!
					print STDERR "Non-N Ambiguous base $res found at position $res_num in the middle gap (gap from $middle[0] to $middle[1]) for id: $primerID.  Not treated as ambiguous.  The --gap parameter may be too wide.\n";
				}
			}
			
			# Save the residue to the consensus sequence array
			push @cons_seq, $res;  

								if ($debug){		
									if (scalar(@$max_key_value) > 1){		# If true, then there is an ambiguous base.  
#										$seqIO->write_aln($aln);
										print STDERR "id: $primerID\tPos: $res_num\n";
#										print STDERR "Max key(s), value:\n";
#										print STDERR Dumper($max_key_value);		
									}
								}
			$res_num++ unless ($res =~ m/-/);	# Residues that have a gap as a consensus (non-ambiguous) are usually indels.  They shouldn't be assigned a position in the final consensus sequence.  
		}	# Done reading each position in the alignment
		
		
#		exit if ($primerID eq "GTGTCCAGGC");		# This was for testing on a primerID group with an insertion in one of the sequences.
		
		# Concatenate the bases to make the consensus sequence.
		my $cons_seq = join "", @cons_seq; 
		
		# Check for any gaps, then collapse gaps
		if ($cons_seq =~ m/[-\.]/ && $verbose){
			# Then there was a gap in the consensus.  This could either be a real deletion, or else one of the reads had an insertion/homopolymer that the other reads did not have.
			my $gap_pos = $-[0]+1;
#			print STDERR "Gap in consensus at position $gap_pos for id: $primerID\nconsensus: $cons_seq\n";			
		}
		$cons_seq =~ s/[-\.]//g; 		# Delete any gaps in the consensus.  
		
		# Deal with ambiguous positions/reads (if any)
		my @ambig_pos = keys %$ambig_pos;

		if ($nongap_ambig_count <= $ambig){		
			# Then report/save read
			my @id_seq = ($primerID, $cons_seq);
			push @good, \@id_seq; 		# Increment the good kept read count
		}
		else {
								if ($debug){	
#									print STDERR "Too many ambiguities ($nongap_ambig_count) in sequence (id $primerID): $cons_seq\n";	
#									$seqIO->write_aln($aln);
#									print STDERR "gap range: $middle[0] to $middle[1]\n";
									print STDERR "ambiguous positions: ";
									print STDERR join " ", @ambig_pos, "\n";
									print STDERR Dumper($ambig_pos); 
								}
			my $id = $primerID . $ambig_summary_for_fasta_header_def;
			my @id_seq = ($id, $cons_seq);
			push @ambiguous, \@id_seq; 	# Increment the ambiguous tossed read count
			push @ambiguous_positions, @ambig_pos; 		# Add the whole array of ambiguous positions to the shared array of ambiguous positions
		}
		
	}
	); 	### End of Parallel Loop
	
	print STDERR "Done making consensus sequences. "; 
	&elapsed($start_time, 'Elapsed', $verbose);
	my $ambig_count = scalar(@ambiguous) || 0;
	my $good_count = scalar(@good) || 0;
	printf STDERR "Consensus reads not reported with > $ambig non-gap ambiguities (or N consensus base outside gap): %2d (%.3f)\n", $ambig_count, $ambig_count/($ambig_count + $good_count) if ($ambig_count);
	printf STDERR "Consensus reads reported: %2d (%.3f)\n", $good_count, $good_count/($ambig_count + $good_count) if($good_count);

	# Save the results
	if (@good){
		print STDERR "Writing consensus sequences to $out_file_good\n";
		my $writefh = open_to_write($out_file_good);
		foreach (@good){
			my ($id,$seq) = @{$_};
			$seq = uc($seq);		# fastx toolkit requires uppercase characters.  
			print $writefh ">$id\n$seq\n";
		}
		close($writefh);
	}
	
	# Write ambiguous file too
	if (@ambiguous){
		print STDERR "Printing ambiguous sequences to $out_file_ambig\n";
		my $ambigfh = open_to_write($out_file_ambig);
		foreach (@ambiguous){
			my ($id,$seq) = @{$_};
			$seq = uc($seq);		# fastx toolkit requires uppercase characters.  
			print $ambigfh ">$id\n$seq\n";
		}
		close($ambigfh);
	}
	
	# Print out a table of ambiguous positions
	if (@ambiguous_positions){
		my %ambig_pos;
		foreach (@ambiguous_positions){
			$ambig_pos{$_}++;
		}
		print STDERR "\nAmbiguous positions (not including N bases outside of gap):\n#Position\tcount\n";
#		foreach my $res_num (sort { $ambig_pos{$b} <=> $ambig_pos{$a} } keys %ambig_pos){		# Sort in descending order by the most prominent position. Or should I sort in ascending order by position?
		foreach my $res_num (sort { $a <=> $b } keys %ambig_pos){		# Sort in ascending order by position
			print STDERR "$res_num\t$ambig_pos{$res_num}";
			if (($res_num >= $middle[0]) && ($res_num <= $middle[1])){
				print STDERR "\tGAP";
			}
			print STDERR "\n";
		}
	}
	
}
#-----------------------------------------------------------------------------
sub read_bam_old {
	my $file = shift;
	my $read_tally;		# Tally of barcode groups.  An array of cleaned sequences for each barcode group.
	print STDERR "Processing file $file ...\n";
	my $sam = Bio::DB::Sam->new(-bam  =>"$file",
		#						-fasta=>"$fasta",		# At this point, I'm calling the reference using $a->dna, which uses the MD tag to calculate the reference sequence, so this fasta file is not required.  
								);
								
	# Process the reads one chromosome at a time
	my @targets = $sam->seq_ids;
	foreach my $chr (@targets){
		my @alignments = $sam->features(-filter  => sub {shift->qual > $mapq},		# This section of code is described on http://search.cpan.org/~lds/Bio-SamTools-1.37/lib/Bio/DB/Sam.pm
										-seq_id  => $chr,
										);
		READ: for my $a (@alignments) {
			# do something with the alignment
			my $seqid  = $a->seq_id;	# chromosome
			my $start  = $a->start;		# 1-based start position on chromosome
			my $end    = $a->end;		# end position on chromosome
			my $strand = $a->strand;	# either 1 (fwd) or -1 (rev)
			my $cig_str= $a->cigar_str;
			my $cigar  = $a->cigar_array;		# Array of arrays with each element of the cigar string, e.g., [ [H, 13], [M, 52], [S, 70], [H, 15] ]
						if ($debug){	print "\n$seqid\t$start\t$end\t$strand\t$cig_str\n";		} 	#print Dumper($cigar);	
			my $name   = $a->name;
			my $flag   = $a->flag;
			
			# where does the alignment start in the query sequence
			my $ref_dna   = uc($a->dna);        # reference sequence bases.  Uses MD tag to calculate the reference sequence.  Alternatively, I could get it right from the sequence using: $dna = $sam->seq($seqid,$start,$end)
			my $query_dna = $a->query->dna; # query sequence bases
			my @scores    = $a->qscore;     # per-base quality scores
			my $mapq      = $a->qual;       # quality of the match (i.e., mapq)
						if ($debug){	print "$seqid\t$start\t$end\n";	}
						if ($debug){	print "que: $query_dna\nref: $ref_dna\n";		}	#shift @$cigar; pop @$cigar; 		foreach (@$cigar){			print $_->[0]." ".$_->[1].", "; 		}		print "\n";
						if ($debug){	print join " ", @scores; print "\n"; print length($query_dna)." ".scalar(@scores)."\n";	}
							
			# Clean up the sequence:
				# Convert low-quality bases to N (based on --baseq option)
				# Remove soft-trimmed portion(s) of sequence from end of according to S elements in cigar string (ignore H elements in cigar)
				# Delete inserted bases (I elements in cigar string)
				# Convert deleted bases to X (from D elements in cigar string)
			
			my $clean_seq = $query_dna;
			# Ignore Hard trimming.  From SAM specification: "H can only be present as the first and/or last operation."
			if ($cigar->[0]->[0] eq "H"){
				shift @$cigar;
			}
			if ($cigar->[-1]->[0] eq "H"){
				pop @$cigar; 
			}
			
			# Low base quality conversion
				# Do this step first, before I change the sequence and make it out of sync with the @scores array
							if ($debug){	print "length: ". length($clean_seq) . ", # scores: ". scalar(@scores) . "\n"; }
			for (my $i=0; $i < @scores; $i++){
				if ($scores[$i] < $baseq){
					# Then replace it with an N.				
							if ($debug){	my $next = $i + 1;  print "replacing from position $i to $next because of score $scores[$i]\n";	}
							if ($debug){	print "bef: $clean_seq\n";	}
					substr($clean_seq, $i, 1) = "N";
							if ($debug){	print "aft: $clean_seq\n";	}
				}
			}
							if ($debug){	print "basq $clean_seq\n";	}

			# Soft-trim sequence removal.  From SAM specification: "S may only have H operations between them and the ends of the CIGAR string."
			if ($cigar->[0]->[0] eq "S"){
				$clean_seq = substr($clean_seq, $cigar->[0]->[1]);		# The behavior of substr is that the second argument is the number of characters to skip, so it should equal the S value. "If LENGTH is omitted, returns everything through the end of the string".  http://perldoc.perl.org/functions/substr.html
							if ($debug){	print "S 5' $cigar->[0]->[1]\n";	}
				shift @$cigar;
			}
			if ($cigar->[-1]->[0] eq "S"){
				$clean_seq = substr($clean_seq, 0, -1 * $cigar->[-1]->[1]); 	# "If LENGTH is negative, leaves that many characters off the end of the string."
							if ($debug){	print "S 3' $cigar->[-1]->[1]\n";	}
				pop @$cigar;
			}
			
			# Read other cigar elements and process.  For M, add to $good_seq_length variable.  For I, delete the bases.  For D, insert bases as X and add length to $good_seq_length.
				# After this, the sequence should be flush with the reference (no indels).
			my $good_seq_length = 0;
			my @clean_scores = @scores;
			foreach my $cigar_element (@$cigar){
							if ($debug){	print Dumper($cigar_element);		}
				if ($cigar_element->[0] eq "M"){	# Match to the reference
					$good_seq_length += $cigar_element->[1];
				}
				elsif($cigar_element->[0] eq "I"){	# Insertion to the reference.  Delete these bases from the read
					# From SAM specification: "Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ"
					# To delete the bases, concatenate the bases before the insertion and the bases after the insertion
					my $bases_before = substr($clean_seq, 0, $good_seq_length); 
					my $bases_after = substr($clean_seq, $good_seq_length + $cigar_element->[1]); 
					$clean_seq = $bases_before . $bases_after;
					
					# Now the scores
							if ($debug){	print join " ", @clean_scores, "\n", @clean_scores[0..$good_seq_length-1], "||", @clean_scores[$good_seq_length + $cigar_element->[1]..$#clean_scores], "\n", $good_seq_length, scalar(@clean_scores) - ($good_seq_length + $cigar_element->[1]), scalar(@clean_scores), "before","after", "total_original\n" ; }
					@clean_scores = @clean_scores[0..$good_seq_length-1,$good_seq_length + $cigar_element->[1]..$#clean_scores];
							if ($debug){	print "$good_seq_length : $bases_before"." "x $cigar_element->[1]."$bases_after before after\n";	}	# This step appears to be working based on this print statement.
				}
				elsif($cigar_element->[0] eq "D"){	# Deletion to the reference.  Put the bases back in as X.
					# To reinsert the bases in the correct position, concatenate the bases before the insertion point to X*N and that to the bases after the insertion point.
					my $bases_before = substr($clean_seq, 0, $good_seq_length); 
					my $bases_after = substr($clean_seq, $good_seq_length); 
							if ($debug){	print "$good_seq_length : $bases_before"." "x $cigar_element->[1]."$bases_after before after\n";	}	# This step appears to be working based on this print statement.
					$clean_seq = $bases_before . 'X' x $cigar_element->[1] . $bases_after;	# Does this work?  Does it add any of the bases around the insertion site twice?
					
					# Now the scores. Give it a score of 3
					my @old_clean_scores = @clean_scores;
					@clean_scores = @old_clean_scores[0..$good_seq_length-1];
					for (my $i=0; $i< $cigar_element->[1]; $i++){	push @clean_scores, '3'; 	}
					push @clean_scores, @old_clean_scores[$good_seq_length..$#old_clean_scores];
							if ($debug){	print join " ", @old_clean_scores, "\n", @old_clean_scores[0..$good_seq_length-1], "| " x $cigar_element->[1], @old_clean_scores[$good_seq_length..$#old_clean_scores], "\n", $good_seq_length, $cigar_element->[1], scalar(@old_clean_scores) - $good_seq_length, scalar(@old_clean_scores) + $cigar_element->[1], "before","inserted", "after", "total_new\n", @clean_scores, "\n" ; }

					$good_seq_length += $cigar_element->[1];
								
				}
				else {
					warn "Cigar element not recognized: " . $cigar_element->[0] . "\n";
				}
							if ($debug){	print "$good_seq_length : $clean_seq\n";	}
			}
							if ($debug){	print "fix: $clean_seq\n";	}	# This sequence should be flush with the reference now.  
	
			
			
			# Now populate the reads and quality scores in $read_tally.  Also get some other information necessary for creating the consensus BAM file, including the chromosome, start (earliest start), flag, read id... mapping qualities
			my @name = split(/:/, $name);
			my $barcode = $name[-1];	#saved as last part of name (read id) by filter_reads_by_barcode_length.pl
			my $id_prefix = $name[0];
			my @clean_seq = split(//, $clean_seq);
			push @{$read_tally->{$barcode}->{sequence}}, \@clean_seq;
			push @{$read_tally->{$barcode}->{quality}}, \@clean_scores;
			push @{$read_tally->{$barcode}->{id_prefix}}, $id_prefix;
			push @{$read_tally->{$barcode}->{chr}}, $seqid;
			push @{$read_tally->{$barcode}->{pos}}, $start;
			push @{$read_tally->{$barcode}->{flag}}, $flag;
			push @{$read_tally->{$barcode}->{mapq}}, $mapq;
			
			if (exists($read_tally->{$barcode}->{reference})){
				unless ($read_tally->{$barcode}->{reference} eq $ref_dna){
					my $warn = 0;
					# Check to see if one is a part of another sequence.  Even being from the same barcode, it might sequence a few more bases in one than the other.  
					if (length($read_tally->{$barcode}->{reference}) > length($ref_dna)){	# New sequence is longer than the current one. The current should be a part of the new one.
						if ($read_tally->{$barcode}->{reference} =~ m/$ref_dna/){
							$ref_dna = $read_tally->{$barcode}->{reference};
						} else {
							$warn++;
							$read_tally->{$barcode}->{bad}++;		# For now let's flag these primerID groups since these are usually truncated sequences anyway.  In the future it would be nice to be able to work with these some more to resolve them a little better since often there will be one bad sequence in a group of 3+ where you could at least use two of them to make a consensus.  
						}
						
					} else {
						if ($ref_dna =~ m/$read_tally->{$barcode}->{reference}/){		# If $read_tally->{$barcode}->{reference} is 
							# nothing.  current $ref_dna is longer so keep it.
						}	else {
							$warn++;
							$read_tally->{$barcode}->{bad}++;		# For now let's flag these primerID groups since these are usually truncated sequences anyway.  In the future it would be nice to be able to work with these some more to resolve them a little better since often there will be one bad sequence in a group of 3+ where you could at least use two of them to make a consensus.  
						}
					}
					if ($warn){		warn "Reference sequence determined to be different for two sequences with the same barcode!  (Keeping the first in the hash, although most likely skipping this primerID group in the merging step.) \n  primID:\t$barcode\n  ref_in_hash:\t$read_tally->{$barcode}->{reference}\n  ref_testing:\t$ref_dna\n";		}
				}
			}
			else {
				$read_tally->{$barcode}->{reference} = $ref_dna;
			}			
		}
	}
	return $read_tally;
}
#-----------------------------------------------------------------------------
sub make_consensus_old {
	# Takes the reads and quality scores (grouped by shared barcode sequence) and produces a single consensus sequence for each barcode group.  
	# I want to optionally produce a representative BAM file where there is one read per barcode group.   What to use for the quality score?  For now, I'll use fisher's statistic to combine the quality scores.  (Alternatively, I could be stringent and choose the weakest p-value.)  Only consider the p-values for the bases that match the consensus (others are considered outliers).
	my ($read_tally,$file) = @_; 	# HoHoA first key, $barcode; second keys 'sequence', 'quality', 'reference', 'id_prefix', 'chr', 'pos', 'flag', 'mapq'; AoA of sequence bases, AoA of quality score integers, or reference sequence as a string
	my $consensus_sequences;	# Arrayref of the consensus sequence

	#Initialize the new consensus bam file
	my ($consensus_bam_prefix,$dir,$ext) = fileparse($file,@suffixes);		
	my $consensus_sam_filename = $save_dir."/".$consensus_bam_prefix . '.consensus.sam';			# SAM file of the consensus sequence		# removed . ".min".$group."pergroup.sam"
	my $consensus_bam_filename = $save_dir."/".$consensus_bam_prefix . '.consensus' ;			# BAM file of the consensus sequence; will be created with samtools	 		# now I'm using samtools sort, so only need prefix, without extension.  	# Also removed ".min".$group."pergroup" . $ext
	my $original_bam = Bio::DB::Bam->open($file ,'r');		# Methods from here http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm#BAM_Files
	my $header = $original_bam->header();
	my $writefh = open_to_write($consensus_sam_filename);
	print $writefh $header->text; 

	foreach my $barcode (keys %$read_tally){
		next if (exists($read_tally->{$barcode}->{bad}));		# i.e., it was labeled as bad in a previous step, i.e., the reference sequence for the reads were different, which means the alignment was different for some of the reads.  
					if ($debug){	print "barcode: $barcode\n";	}
		my $number_seqs = scalar(@{$read_tally->{$barcode}->{sequence}});
		next if ($number_seqs < $min_reads);		# Only consider those barcode groups with at least 2 sequences.
		my @ref = split(//, $read_tally->{$barcode}->{reference});
#				print join " ", @ref; 
		my @consensus;		#Array containing the bases in the new consensus
		my @combined_qualities;	#Array containing the combined qualities scores in order
		
		for (my $i=0; $i < scalar(@ref); $i++){		# Walk through each base in all of the sequences
			
			# Get the majority rule base (A, C, T, G)
			my %base_tally;		# tally of bases at this position
			for (my $j=0; $j< $number_seqs; $j++){  # Look at the base at position $i, in each sequence, $j
				if ($read_tally->{$barcode}->{sequence}->[$j]->[$i]){
					my $base = $read_tally->{$barcode}->{sequence}->[$j]->[$i];
					$base_tally{$base}++;
				} else {
#							print Dumper($read_tally->{$barcode});
				}
			}
			if ($debug){	print Dumper(\%base_tally);		}
			my $max_keys = find_key_with_biggest_value(\%base_tally, 0, 1);	# Returns an arrayref of the bases with the maximum value.  
			my $good_base = "";
			if (scalar @$max_keys > 1){	
				# A tie for the most represented base Take REFERENCE sequence at this position
				# Could be a 2-way, 3-way or 4-way tie.  E.g., with three sequences and base tallies A 1, C 1, G 1.   
				#	push @consensus, $ref[$i];
				$good_base = $ref[$i];
				
				# Just out of curiosity, what would be the base if we used quality score to decide?
				my $best_quality_base = $read_tally->{$barcode}->{sequence}->[0]->[$i];	# Take base from first sequence to start.  Look for something with a better score.
				my $best_quality_score = $read_tally->{$barcode}->{quality}->[0]->[$i];
				for (my $a = 1; $a < scalar(@{$read_tally->{$barcode}->{quality}}); $a++){
					if ($read_tally->{$barcode}->{quality}->[$a]->[$i] && ($read_tally->{$barcode}->{quality}->[$a]->[$i] > $best_quality_score) ){
						$best_quality_base = $read_tally->{$barcode}->{sequence}->[$a]->[$i];
						$best_quality_score = $read_tally->{$barcode}->{quality}->[$a]->[$i];
					}
					elsif($read_tally->{$barcode}->{quality}->[$a]->[$i] && ($read_tally->{$barcode}->{quality}->[$a]->[$i] == $best_quality_score) ) {
						# Tie for highest so far...
							if ($debug){	print "Tie for highest quality score so far.  $best_quality_base $best_quality_score, $read_tally->{$barcode}->{sequence}->[$a]->[$i] $read_tally->{$barcode}->{quality}->[$a]->[$i]\n";		}
						# If this is the last one and it's a tie, and we are using quality scores to decide, we should probably refer to the REFERENCE base here.  If one of the bases is a reference base, then use it.  If not, then not sure what to do...
					}
				}
								if ($debug){	
									print "Multiple bases in reads: ", join " ", @$max_keys, "Taking reference base:", $ref[$i]; 
									print "\nBest base based on quality: $best_quality_base $best_quality_score\nquality scores: "; 
									for (my $j=0; $j< $number_seqs; $j++){ 
										print "$read_tally->{$barcode}->{sequence}->[$j]->[$i] $read_tally->{$barcode}->{quality}->[$j]->[$i], ";	
									}  
									print "\n";		
								}
			}
			else {
					#push @consensus, $max_keys->[0];
				$good_base = $max_keys->[0];
			}
			push @consensus, $good_base;	# Does this work okay?
			
			# Get the combined quality score of those that have this base.  Use 
			my @good_seq_indices;	# indexes for the reads that have the consensus base
			my @good_quality_scores;	# quality scores for the bases that have the consensus sequence
			for (my $j=0; $j< $number_seqs; $j++){  # Look at the base at position $i, in each sequence, $j
				if ($read_tally->{$barcode}->{sequence}->[$j]->[$i] && $read_tally->{$barcode}->{sequence}->[$j]->[$i] eq $good_base){
					push @good_seq_indices, $j;		
					push @good_quality_scores, $read_tally->{$barcode}->{quality}->[$j]->[$i];		
				}
			}
				# Now convert to p-values and combine ... 
			my @pvalues = map { convert_qual_to_pvalue($_) } @good_quality_scores; 
			my $quality_ascii = '$';	#default low quality, quality of 3, which is 36 in Sanger scale.
			if (@good_quality_scores){ 
				if ($debug){	print "qual scores: ", join " ", @good_quality_scores, "\npvalues:", @pvalues, "\n";	}
				$quality_ascii = convert_pvalue_to_ascii_qual(combine_pvalues(\@pvalues));	
			}
			else {
				if ($debug){						warn "No quality scores because none of the reads matched the consensus!  Assigning quality score of 3\n"; 				}	
			}			
			push @combined_qualities, $quality_ascii; 
					
		}
		
		# Get the other SAM items
		my @flags = @{$read_tally->{$barcode}->{flag}};
		my @chrs = @{$read_tally->{$barcode}->{chr}};
		my @id_prefix = @{$read_tally->{$barcode}->{id_prefix}};
		my @pos = @{$read_tally->{$barcode}->{pos}};
		my @mapq = @{$read_tally->{$barcode}->{mapq}};
		
		# Check those that should be the same for SAM
		foreach my $array (\@flags, \@chrs, \@id_prefix, \@pos){
			if (check_if_array_items_different($array)){
				warn "Some items of the same barcode group are different (barcode $barcode): ", join " ", @$array, "\n";
			}
		}
		my $flag = $flags[0];
		my $chr = $chrs[0];
		my $id_prefix = $id_prefix[0];
		
		# Earliest position for SAM
		my $pos = range(\@pos, 'min');	# this works when the read aligns to the forward strand; does it work for reverse strand?
		
		# Combined mapq for SAM id
		my $mapq = range(\@mapq, 'min');	# stringent
		my $mapq_concat = join ",", @mapq;
				
		# Now save the consensus bases and qualities
		my $consensus = join "", @consensus;
		push @$consensus_sequences, $consensus;		#** This is the part that we'll be returning from the subroutine
		my $qualities = join "", @combined_qualities;
		
		# Cigar, sequence, id for SAM
		my $id = $id_prefix . ":" . $barcode . ":" . scalar(@mapq) . ":" . $mapq_concat;
		my $cigar = length($consensus) . "M";
		(my $seq = $consensus) =~ s/X/N/g; 	# Does this work? Yes.
					if ($debug){	print "seq: \t$seq\ncons: \t$consensus\n";	}
					
		# Print out SAM record
		print $writefh join "\t", $id, $flag, $chr, $pos, $mapq, $cigar, '*', 0, 0, $consensus, $qualities;
		print $writefh "\n";
		
		
	}
	close($writefh);
	system($PWD."/samtools view -bS $consensus_sam_filename | samtools sort - $consensus_bam_filename ");
	system($PWD."/samtools index $consensus_bam_filename" . ".bam");
	return $consensus_sequences;
}
#-----------------------------------------------------------------------------
sub check_if_array_items_different {
	# Checks to see that the items are all the same
	# Reports 1 if they are different
	my $array = shift;
	my %hash;
	foreach my $item (@$array){		# alternatively, my %hash = map { $_ => 1} @$array; 
		$hash{$item}++;
	}
	my @values = sort {$a <=> $b} values %hash;
	if (scalar(@values) > 1){
		# There are more than one item
		return 1;
	}
	else {
		# All the same item
		return 0;
	}
}
#-----------------------------------------------------------------------------
sub combine_pvalues {
	# Combines multiple p-values into one combined p-value.  
	# Implementation of fisher's method
	# http://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/
	# R code: fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
	# http://en.wikipedia.org/wiki/Fisher's_method
		# log base e
	my $pvalues = shift;	# arrayref of pvalues.
	my $df = 2 * scalar (@$pvalues);
	my @log_pvalues = map {log($_)} @$pvalues;
			if ($debug){	print "pvalues: ", join " ", @$pvalues, "\nlog pvalues: ", @log_pvalues, "\n";	}
	my $sum_log_pvalues = total(\@log_pvalues);
	my $chisq = -2 * $sum_log_pvalues;
	my $combined_pvalue = Statistics::Distributions::chisqrprob ($df,$chisq);
#					print "pvalues: ", join " ", @$pvalues, "\nlog pvalues: ", @log_pvalues, "\n";
			if ($debug){	print "chisq: $chisq combined p-value: $combined_pvalue\n";	 }
	unless ($combined_pvalue){	
		# combined p-value is 0.  This happens when you have ~30 reads in a barcode group, so the combined p-value is too strong.  
		# This will produce an error when converting to ascii qual, which will take the log value of this.
		$combined_pvalue = 3.162e-26; 		# qual score of this p-value is 255.  
	}
	
	return $combined_pvalue;
}
#-----------------------------------------------------------------------------
sub convert_pvalue_to_ascii_qual {
	my $pvalue = shift;
	my $qual = round(-10 * log10($pvalue) );
	unless ( ($qual <= 41) && ($qual >= 0) ){
		if ($debug){	print "quality score out of expected range: $qual\n";	}
	}
	$qual += 33;	# Add 33 for Sanger
	# right now, from the combined pvalues, I'm getting quality scores of 60, 85, etc. (even before adding 33 for Sanger).  I wonder if that will mess up samtools or other downstream programs... 
	# I at least need to have a cap of 126 or else we'll get into weird symbols not on the keyboard
	$qual = 126 if ($qual > 126);		# Was 126 or 41
	my $quality_ascii = chr($qual);
			if ($debug){	print "ascii quality: $quality_ascii\n";	}
	return $quality_ascii;	
}
#-----------------------------------------------------------------------------
sub convert_ascii_qual_to_pvalue {
	# Takes a single ascii quality and returns the p-value.  Assumes Sanger scale.
	# To walk through an array of ascii quality scores, use map, e.g.,
	# my @pvalues = map { convert_ascii_qual_to_pvalue($_) } @good_quality_scores;
	my $quality_ascii = shift;		#takes a single ascii quality
	my $quality = ord($quality_ascii) - 33;		# using Sanger phred scale, which adds 33 to the score
	my $pvalue = convert_qual_to_pvalue($quality);

	return $pvalue;
}	
#-----------------------------------------------------------------------------
sub convert_qual_to_pvalue {
	# Takes quality scores in the 0-40 range and outputs the p-value
	my $quality = shift;
	my $pvalue = 10 ** ($quality / (-10) );
	return $pvalue; 
}
#-----------------------------------------------------------------------------
sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
#-----------------------------------------------------------------------------
sub iupac_ambiguities {
	# Sort the multiple bases present alphabetically, concatenate, and then lookup the ambiguity code in the hash.
	# Order: ACGT
	# http://en.wikipedia.org/wiki/Nucleic_acid_notation
	# http://droog.gs.washington.edu/parc/images/iupac.html
	my %hash = (
		'AC'	=>	'M',
		'AG'	=>	'R',
		'AT'	=>	'W',
		'CG'	=>	'S',
		'CT'	=>	'Y',
		'GT'	=>	'K',
		'ACG'	=>	'V',
		'ACT'	=>	'H',
		'AGT'	=>	'D',
		'CGT'	=>	'B',
		'ACGT'	=>	'N',
		'ac'	=>	'm',
		'ag'	=>	'r',
		'at'	=>	'w',
		'cg'	=>	's',
		'ct'	=>	'y',
		'gt'	=>	'k',
		'acg'	=>	'v',
		'act'	=>	'h',
		'agt'	=>	'd',
		'cgt'	=>	'b',
		'acgt'	=>	'n',
	);
	return \%hash;
}
#-----------------------------------------------------------------------------

