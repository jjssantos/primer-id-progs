#!/usr/local/bio_apps/perl-5.16.2/bin/perl
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
use Bio::DB::Sam;		#http://search.cpan.org/~lds/Bio-SamTools-1.37/lib/Bio/DB/Sam.pm
use Bio::SeqIO;
use Array::Utils qw(:all);
use File::Slurp qw( prepend_file );
use Parallel::Loops;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;
use File::Which;

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

my $PWD = pwd_for_hpc();
my $mafft_bin = $PWD.'/mafft'; # philip

my $save;
my $files;
my $verbose;
my $output;
my $gzip;
my $baseq;
my $mapq;
my $ref;
my $gff;
my $label;
my $debug;
my $fullpeptide;
my $prefix;
my $variant_threshold;
my $alignment_length;
my $cpu = 1; 
my $clustalw;
GetOptions('save=s' => \$save, 'output=s' => \$output, 'verbose' => \$verbose, 'files=s' => \$files, 'gzip' => \$gzip, 'baseq=s' => \$baseq, 'mapq=s' => \$mapq, 'ref=s' => \$ref, 'r=s' => \$ref, 'gff=s' => \$gff, 'label=s' => \$label, 'labels=s' => \$label, 'debug' => \$debug, 'fullpeptide=s' => \$fullpeptide, 'prefix=s' => \$prefix, 'variant_threshold=s' => \$variant_threshold, 'alignment_length=s' => \$alignment_length, 'cpu=s' => \$cpu, 'p=s' => \$cpu, 'clustalw' => \$clustalw);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "
convert_reads_to_amino_acid.pl takes one or more fasta files (e.g., sequencing reads), 
compares to a reference coding sequence and outputs frequency tables of nucleotide, 
codon, and amino acid sequence by position.  A merged table containing frequencies for 
all three is also produced, as well as phylip files for nucleotide and peptide sequence.  
All reads must be in the same orientation as the reference.  Reads are translated using 
the mammalian genetic code, including interpretation of wobble bases where the third 
position of the codon is an unknown base. In addition to amino acids, frequency for stop, 
unknown amino acid, and deletion in the read compared to reference are encoded as _, X, 
and *, respectively. 

Pre-requisites.
Programs required on PATH:
fasta_collapser.pl 
mafft
clustalw (if using --clustalw)

Non-standard Perl modules:
Bioperl, including Bio::DB::Sam and Bio::Tools::Run::Alignment::Clustalw
Array::Utils
Parallel::Loops


Required Arguments:
--files 	Input file(s), comma-delimited.  One fasta file per sample.  Required.
-r/--ref	Reference coding (CDS) sequence fasta file.  Required.  Used to determine 
		translation frame of reference.

Optional Arguments:
--label		Name(s) for sample(s), comma-delimited (no spaces) in the order in which 
		they were input into --files. Default: Sample_1,Sample_2, etc.
--prefix	Prefix for output files.  Default = Convert_reads
--mapq		Minimum read mapping quality to consider.  Default = 0.
--baseq		Minimum base quality to consider for translation.  Default = 0.
--save		Directory in which to save files. Default = pwd.  If folder doesn\'t exist, it 
		will be created.
-p/--cpu	Number of processors to use. Default = 1.
--clustalw	Use clustalw to run alignments. Default is to use mafft.  *Executables must be 
		on the PATH.*  Speed is about the same for clustalw and mafft (about 1.3x faster 
		for clustalw).
--alignment_length	Length for alignment from the starting position.  Default is the 
		majority read length.  Reads shorter than this will be removed and reads longer
		will be truncated.  
--variant_threshold	Minimum threshold frequency for output variant file.  Default = 0 
		(i.e., print all variants).  (Note that in calculate_linkage_disequilibrium.pl
		there is an option to filter for a variant_threshold as well.)

Example:
convert_reads_to_amino_acid.pl --files 30_S1.contigs.pid.btrim.0.majority.cons.fasta --ref CAL0409_HA_cds.fa --label 30_S1 --prefix 30_S1.contigs.pi.btrim.0.majority.cons -p 7

Notes:
1. Recommended h_vmem (RAM) value per thread is 6G, although it may be fine with less, depending on the size and complexity of your input file. 
";



# To Do:
# Add option to skip bases that have an indel (or to trim it at the 3' end right before the indel).
# Add ability to process genes/reads on - strand 
# I should also somehow store those reads or positions that affect the translation, such as early stop codons, changing initial methionine, or changing the stop codon to another amino acid (continued translation).  Perhaps add the option to ignore amino acids following an early stop and/or labeling the amino acids of a read as untranslated if the initial methioinine is mutated...  (unfortunately, the other reads from that same quasispecies will be counted as translated, since we have short sequences and they won't all contain the initial Met.), and/or to add the amino acids following an altered stop codon.
# Add columns Largest_Nonref_Count, Largest_Nonref_AA and Largest_Nonref_as_Fraction_of_Nonref
# Maybe simplify the table so there aren't so many columns.  Or, I could produce a summary table that has reference, consensus, major variant, minor variants by codon.  
# Query intergenic and non-coding regions too for nucleotides (right now, it is just looking at CDS regions)
# Maybe save any reads that were soft-trimmed to see what sequence is there... especially if the sequence quality is good
# Adapt script to get splice junction amino acids to use with mRNA-seq.  Right now it only works for genomic DNA or genomic RNA sequencing (i.e., introns are present).  
# Allow a read to overlap with multiple exons.  Right now it assumes that each read will only overlap with a single exon; it finds the first exon that overlaps and then trims the read to the ends of that exon. As is, it is fine for sequencing influenza H1N1 with <250bp reads, as the introns are ~475 and ~690 bp long.  As long as the template is genomic DNA, the reads won't overlap two exons; if the template is mRNA, then they could very easily overlap two exons.  Also, I would possibly have to align with TopHat instead of BWA.  Will says template is genomic RNA and probably always will be.  
	# Might be able to use gff output to splice together exon sequences.  Also cigar output (convert Ds to Is) for gaps, etc.  
# Make the clustalw Perl modules and some others optional based on input options (eval...)

# Output stats about each codon present as well.  Should probably be a separate file.  
	# Having the amino acid-level analysis will also be useful, but I can do this in excel for now. We would most certainly want to do  synonymous/non-synonymous comparisons for each codon at some point. 
	# Tally up amino acids.  
# Ways to potentially make it more amenable to a multi-threaded pipeline
	# Print out separate file of reads/peptides per gene and run linkage disequilibrium for each gene in parallel
	# Also run nuc, aa, and codon linkage disequilibrum and/or dN/dS in parallel
# * Print out graphs for frequencies so we can see variants in positions.  
# Use Exonerate to deal with exons.  Align CDS to reads with est2genome.
# Change --fasta to --ref
# Maybe let BAM be an option as input... (keep read_bam and some other subroutines that work with bam files)
# Maybe let FASTQ be an option as input too/instead?

# Change log
# 2013-01-15
# Allowing multiple genes on a chromosome.  The reads will be processed for each gene.
# Also allows genes with multiple exons.  
# 2013-01-21
# Removed phase portion and added code to parse codon_start attribute of GFF file to fix CDS start coordinate of first exon.
# 2013-01-22
# Fixed bug where intronic reads were being processed even though no overlap with exons.  Made an initial check to see that they intersect with an exon (at least 2bp) before cleaning up the sequence and processing.
# Fixed bug where stop codon amino acids (_) were not being tallied as "Ref" or "Nonref", resulting in inflated percentNonref values.  
# 2013-05-09
# Added --fullpeptide option.  
# 2013-05-28
# Lots of changes!  Nearly a full rewrite.  Including adding multiple reports (nuc, aa, codon frequency tables, merged table, variants above a threshold frequency, and cleanpeptides and cleanreads output files)
# Didn't change the method to clean up the reads (in read_bam).  That stayed the same.  
# Calculate consensus reference, make comparisons against that.    
# Count up the codons, report synonymous, nonsyn with respect to the consensus.  
# Make it fit this description: "tally up nucleotides, codons, and amino acids from BAM file reads with convert_reads_to_amino_acids.pl (add codon features) -> tables of nucleotide frequency, codon/amino acid frequency, and merged table"
	# Put "statistical comparison of samples from MpileupFrequency.pl" and "table with p-values of all two-way comparisons." in a separate script using the output of this script.  
# 2013-05-29
# Removed this option, made it default.
# 	--fullpeptide	Print out the full peptide sequence into this file.  e.g., 
#			--fullpeptide peptides.txt.  These will be cleaned with respect to the reference.
# 2013-06-05
# Added --variant_threshold option.  Changed default to 0.  These can be filtered in subsequent script, calculate_linkage_disequilibrium.pl.  
# Added back --save option
# 2013-06-05
# Added positionWithinCodon column for merged table.  
# Added --alignment_length option and $majority_length variable
# 2013-09-23
# Changed read_bam subroutine to use '|' as delimiter (instead of '_') for cleanpeptides.txt and cleanreads.txt output files to separate readID, gene, position.
# 2013-09-27
# Changed parse_gff to not skip 'X' reference amino acids and the corresponding codon (usually resulting from an ambiguous base in the reference).  
# Also added && $type ne 'aa' in the check for assigning the variant filter "VARIANT_N" and "CONSENSUS_N" so that we don't give this non-PASS status to amino acids.  
# 2013-12-13
# Major changes to accomodate changes such as using MiSeq reads and allowing indels now. Also, the previous script, merge_primerid_read_groups.pl was modified a lot, so the output it produces is different -- now fasta.  I'll still align it to the reference with bwa to get BAM files (unless it is not doing a very good job...) but there are no quality scores at the moment.   
# Changed default $mapq and $baseq to 0 since we've done some pretty good filtering up to this point.  
# 2013-12-17
# Continued some major changes. Simplified read cleanup step to not worry about trimming around exons.  Assuming the user has aligned the reads to the CDS sequences, so any introns will be soft-mapped on the ends of the reads and that there are NO full introns in the sequences.  
# 2014-01-02
# More major changes.  Re-aligning each of the reads to the CDS reference nucleotide sequence instead of using the CIGAR string to determine position of indels and frameshifts.  Now takes a fasta file for input reads (--files). New features in subroutine read_input_reads.  Meant to only process one gene/CDS, so no need for gff file.   
# Added --cpu for processing fasta files.
# Removed DNDS stuff.
# 2014-01-06 
# More changes to read_input_reads subroutine
# Removed read_input_bam subroutine (saved in a separate file)
# 2014-01-07
# More changes to read_input_reads subroutine
# Fixed the reports to read the new format for the data structure $nuc_aa_codon_tally
# Added sub get_cds_sequences to read and translate the reference sequence for storage in the global data structure $cds
# Added a column to the reports to note when the consensus base/codon/aa is different from the reference sequence
# 2014-01-09
# Changed temporary directory to be in /tmp/ instead of Cwd.  
# 2014-01-10
# Changed --fasta input to --ref (and changed one section where I had used variable $ref to $ref_sequence)
# Changed output extensions to use dots '.' instead of underscores.
# 2014-01-14
# Changed variants report to print out in order nuc codon aa to match the linkage report from calculate_linkage_disequilibrium.pl
# 2014-04-01
# Changed read_input_reads() sub so it won't start keeping track of nucleotide positions until we've reached the first full codon (to match the behavior with saving the nucleotide base) so the base position is synchronized with the corresponding base 
# Added fastx_collapser step in get_unique_seqs() sub so we can work with unique sequences (and keep track of their counts)
# Also added get_count_from_id() sub to work with fasta with increment by the read counts instead of 1 per sequence.  
# 2014-06-12
# Swapped out fastx_collapser for a custom script I wrote fasta_collapser.pl which allows for ambiguous characters, which sometimes creep into the consensus fasta file from merge_primerid_read_groups.pl (usually when they are in the middle gap).  
# 2014-07-18
# Added count_total_from_good_toss_arrays subroutine.  
# 2015-03-03
# Removed --gff option (not using it at the moment; using CDS reference sequence to get the frame of reference for translation).  
# Also removed --baseq and --mapq options because not using them now.
#--gff		GFF file indicating coding start and stop within the fasta file.
#		Must have \"CDS\" features.  One gene per chromosome.  Should be
#		sorted by coordinate.  (Used to get the gene names.)
#--mapq		Minimum read mapping quality to consider.  Default = 0.
#--baseq		Minimum base quality to consider for translation.  Default = 0.
# Edited usage statement.  Removed this:
# It is recommended that all reads begin with the same starting position.  For now, the script assumes there is one gene per chromosome and genes must be on the positive strand (e.g., viral genomes).  
# (Amino acid * usually results from deletions filled in with base X.).  

unless ($ref && ($files||$ARGV[0]) ){	print STDERR "$usage\n";	exit;	}	#fasta and gff are required
unless($files){	
	if (scalar(@ARGV)>1){	#This will allow you to use wildcard to pass files to the script from the command line. e.g, script.pl *.txt, and it will run through each.
		$files = join ",", @ARGV;
	}
	else{
		$files = $ARGV[0];	
	}
}

$mapq //= 0;	# zero is a valid value, hence //= ("defined or") operator.  Doesn't work in Perl 5.8.		 # Was 30
$baseq //= 0;	# zero is a valid value		# Was 13

print STDERR "baseq: $baseq\nmapq: $mapq\n"; 	# exit;
my $save_dir = $save || Cwd::cwd();
# print $save_dir;
# exit;
unless (-d $save_dir){	mkdir($save_dir) or warn "$!\n"	}
# warn "save directory: $save_dir\n";

my $start_time = time;
my @files = &aomisc::get_files($files);		#If allowing a directory, specify extension of the files in second argument, e.g., my @files = &get_files($files, 'bed');
my @suffixes = (qw(.bed .bed.gz .bed12 .bed12.gz .txt .txt.gz .BED .BED.gz .BED12 .BED12.gz .fasta .fa .FA .FASTA .FAS .fas), @SUFFIXES);	#for fileparse.  Feel free to add more accepted extensions.  @SUFFIXES comes from aomisc.pm.  

my @labels;
if ($label){
	@labels = split (/,/, $label);
}
else {
	for (my $i = 1; $i < (scalar(@files)+1); $i++){
		push @labels, "Sample_".$i;
	}
}

# Defaults for input options
$prefix ||= 'Convert_reads';
$variant_threshold //= 0;		# Was 0.0001 (1 variant out of 10000 sequences).  0 is a valid value, hence // instead of ||.


# Get reference sequence

my ($refseq,$gene);		# Reference coding nucleotide sequence
if ($ref){
	my $in = Bio::SeqIO->new(-file => $ref , -format => 'Fasta' );		# There should only be one sequence in this file, so it should be the first.   # This conflicts with the Bioperl module
	my $seq = $in->next_seq();
	$refseq = $seq->seq();
	$gene 	= $seq->display_id();
					if ($verbose){	print "refseq:\n$gene\n$refseq\n";	}
}

my %converter = (		# Modified from http://www.wellho.net/resources/ex.php4?item=p212/3to3
    'TCA' => 'S', # Serine
    'TCC' => 'S', # Serine
    'TCG' => 'S', # Serine
    'TCT' => 'S', # Serine
    'TCN' => 'S', # Serine wobble
    'TTC' => 'F', # Phenylalanine
    'TTT' => 'F', # Phenylalanine
    'TTA' => 'L', # Leucine
    'TTG' => 'L', # Leucine
    'TAC' => 'Y', # Tyrosine
    'TAT' => 'Y', # Tyrosine
    'TAA' => '*', # Stop
    'TAG' => '*', # Stop
    'TGC' => 'C', # Cysteine
    'TGT' => 'C', # Cysteine
    'TGA' => '*', # Stop (or Selenocysteine, U)
    'TGG' => 'W', # Tryptophan
    'CTA' => 'L', # Leucine
    'CTC' => 'L', # Leucine
    'CTG' => 'L', # Leucine
    'CTT' => 'L', # Leucine
    'CTN' => 'L', # Leucine	wobble
    'CCA' => 'P', # Proline
    'CCC' => 'P', # Proline
    'CCG' => 'P', # Proline
    'CCT' => 'P', # Proline
    'CCN' => 'P', # Proline wobble
    'CAC' => 'H', # Histidine
    'CAT' => 'H', # Histidine
    'CAA' => 'Q', # Glutamine
    'CAG' => 'Q', # Glutamine
    'CGA' => 'R', # Arginine
    'CGC' => 'R', # Arginine
    'CGG' => 'R', # Arginine
    'CGT' => 'R', # Arginine
    'CGN' => 'R', # Arginine wobble
    'ATA' => 'I', # Isoleucine
    'ATC' => 'I', # Isoleucine
    'ATT' => 'I', # Isoleucine
    'ATG' => 'M', # Methionine
    'ACA' => 'T', # Threonine
    'ACC' => 'T', # Threonine
    'ACG' => 'T', # Threonine
    'ACT' => 'T', # Threonine
    'ACN' => 'T', # Threonine wobble
    'AAC' => 'N', # Asparagine
    'AAT' => 'N', # Asparagine
    'AAA' => 'K', # Lysine
    'AAG' => 'K', # Lysine
    'AGC' => 'S', # Serine
    'AGT' => 'S', # Serine
    'AGA' => 'R', # Arginine
    'AGG' => 'R', # Arginine
    'GTA' => 'V', # Valine
    'GTC' => 'V', # Valine
    'GTG' => 'V', # Valine
    'GTT' => 'V', # Valine
    'GTN' => 'V', # Valine wobble
    'GCA' => 'A', # Alanine
    'GCC' => 'A', # Alanine
    'GCG' => 'A', # Alanine
    'GCT' => 'A', # Alanine
    'GCN' => 'A', # Alanine wobble
    'GAC' => 'D', # Aspartic Acid
    'GAT' => 'D', # Aspartic Acid
    'GAA' => 'E', # Glutamic Acid
    'GAG' => 'E', # Glutamic Acid
    'GGA' => 'G', # Glycine
    'GGC' => 'G', # Glycine
    'GGG' => 'G', # Glycine
    'GGT' => 'G', # Glycine
    'GGN' => 'G', # Glycine wobble
    );

# Parse GFF file to get start and stop, and translated ORF
my $genes;
#$genes = parse_gff($gff);

# Instead of parsing GFF, use the reference CDS sequence to get the coding nucleotides, codons, and amino acids by position.
my $cds; 
$cds = get_cds_sequences($refseq); 

#		print Dumper($genes); exit;

my @NUC = qw(A C G T N);
my @AA = qw(A C D E F G H I K L M N P Q R S T V W Y * X -); 	# All amino acids, and * for stop, X for unknown, and - for deletion in read compared to reference
my @CODON;
my $codon_hash;
foreach my $codon (sort keys %converter){
	push @{$codon_hash->{$converter{$codon}}}, $codon . ":" . $converter{$codon}; 
}
foreach my $aa (@AA){		# Make the codons in the report in the same order as the amino acids in @AA, and group codons by the amino acid they encode
	if (exists($codon_hash->{$aa})){
		foreach my $codon_aa (@{$codon_hash->{$aa}}){
			push @CODON, $codon_aa;
		}
	}
}

#			print join " ", @CODON; print "\n"; exit;

# Parse the BAM files and print reports

my $nucleotide_report 	= $save_dir . "/" . $prefix . ".nuc.tally.xls";
my $codon_report 		= $save_dir . "/" . $prefix . ".codon.tally.xls";
my $amino_acid_report 	= $save_dir . "/" . $prefix . ".aa.tally.xls";
my $merged_report 		= $save_dir . "/" . $prefix . ".merged.tally.xls";
my $variants			= $save_dir . "/" . $prefix . ".variants.minfreq".$variant_threshold.".xls";		# Maybe add threshold to the name.  Or in the header of the file.

my $nucleotide_report_fh 	= open_to_write("$nucleotide_report");
my $codon_report_fh 		= open_to_write("$codon_report");
my $amino_acid_report_fh 	= open_to_write("$amino_acid_report");
my $merged_report_fh		= open_to_write("$merged_report");
my $variants_fh 			= open_to_write("$variants");

print $nucleotide_report_fh "#name\t";
print $nucleotide_report_fh join "\t", qw(gene nucleotidePosition refNucleotide consensusNucleotide refDiffConsensus Sample coverageDepth);
print $nucleotide_report_fh "\tnum"; 
print $nucleotide_report_fh join "\tnum", @NUC, "Other";
# print $nucleotide_report_fh "\t"; 
# print $nucleotide_report_fh join "\t", qw(del SumACGT numNonzeroACGT Ref Nonref percentNonref);
print $nucleotide_report_fh "\n";

print $codon_report_fh "#name\t";
print $codon_report_fh join "\t", qw(gene codonPosition refCodon consensusCodon refDiffConsensus Sample coverageDepth);
print $codon_report_fh "\tnum"; 
print $codon_report_fh join "\tnum", @CODON, "Other";
# print $codon_report_fh "\t"; 
# print $codon_report_fh join "\t", qw(del SumACGT numNonzeroACGT Ref Nonref percentNonref);		# Maybe syn, nonsyn; or maybe save that for the summary report.  
print $codon_report_fh "\n";

print $amino_acid_report_fh "#name\t";
print $amino_acid_report_fh join "\t", qw(gene aminoAcidPosition refAminoAcid consensusAminoAcid refDiffConsensus Sample coverageDepth); 
print $amino_acid_report_fh "\tnum"; 
print $amino_acid_report_fh join "\tnum", @AA, "Other";  		# Needed to add "num" before amino acids, since R converts _ to "X_" and * to "X."
# print $amino_acid_report_fh "\t"; 
# print $amino_acid_report_fh join "\t", qw(SumKnownAminoAcids numNonzeroKnownAminoAcids Ref Nonref percentNonref);
print $amino_acid_report_fh "\n"; 

print $merged_report_fh "#name\t";
print $merged_report_fh join "\t", qw(gene nucleotidePosition positionWithinCodon refNucleotide consensusNucleotide refDiffConsensus aminoAcidPosition refCodon consensusCodon refAminoAcid consensusAminoAcid Sample coverageDepth); 
foreach my $type (qw(Nucleotide Codon AminoAcid)){
	foreach my $rank (qw(top second third other)){
		print $merged_report_fh "\t$rank"."$type";		#topNucleotide secondNucleotide otherNucleotide countTopNucleotide countSecondNucleotide
	}
	foreach my $rank (qw(top second third other)){
		print $merged_report_fh "\tcount".ucfirst($rank)."$type";
	}
}
print $merged_report_fh "\n";		#	#name	chr	gene	nucleotidePosition	aminoAcidPosition	refAminoAcid	consensusAminoAcid	Sample	coverageDepth	topNucleotide	secondNucleotide	thirdNucleotide	otherNucleotide	countTopNucleotide	countSecondNucleotide	countThirdNucleotide	countOtherNucleotide	topCodon	secondCodon	thirdCodon	otherCodon	countTopCodon	countSecondCodon	countThirdCodon	countOtherCodon	topAminoAcid	secondAminoAcid	thirdAminoAcid	otherAminoAcid	countTopAminoAcid	countSecondAminoAcid	countThirdAminoAcid	countOtherAminoAcid


print $variants_fh "## Variants above threshold frequency $variant_threshold\n";
print $variants_fh "#Type\tgene\tposition\tconsensus\tvariant\tcount\tcoverage\tfrequency\tfilter\n";


for (my $i = 0; $i < @files; $i++){
	my $file = $files[$i];
	
	my $unique_seqs_file = get_unique_seqs($file);
	my $nuc_aa_codon_tally = read_input_reads($unique_seqs_file);
#			print Dumper($nuc_aa_codon_tally); exit;
	print_reports($nuc_aa_codon_tally, $labels[$i]);		
	print_merged_report($nuc_aa_codon_tally, $labels[$i]);	
	print_variants($nuc_aa_codon_tally, $labels[$i]);
}


close($nucleotide_report_fh);
close($codon_report_fh);
close($amino_acid_report_fh);
close($merged_report_fh);
close($variants_fh);

&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub get_unique_seqs {
	my $file = shift;
	
	# Get the unique sequences in a fasta file with counts and a hash with number of times the sequence was found in the file = "count".  
	# This is designed to speed up the computation significantly, since many of the reads will be identical.  
	# Store the count for the sequence in the definition so it can be retrieved and stored in the seq_info hash (inside the parallel loop)
	# In fact, use custom script to do this
	# fasta_collapser.pl -i some.fasta -o some.unique.fasta
	# Output example:
	# >1-16595
	# ATTCGAAAGATTCAAAATATTTCCCAAAGAAAGCTCATGGCCCGACCACAACACAAACGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCNNNNNNNTAACAGTAAGGAACAACAGAATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
	# >2-870
	# ATTCGAAAGATTCAAAATATTTCCCAAAGAAAGCTCATGGCCCGACCACAACACAACCGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCNNNNNNNTAACAGTAAGGAACAACAGAATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
	# >3-657
	# ATTCGAAAGATTCAAAATATTTCCCAAAGAAAGCTCATGGCCCGACCACAACACAAGCGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCNNNNNNNTAACAGTAAGGAACAACAGAATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
	# ...
	# The id has seq number - count ; it's sorted by most common to least common sequence

	my ($filename,$dir,$ext) = fileparse($file,@SUFFIXES);
#	my $unique_seqs_fasta = $dir . $filename . ".uniq" . $ext; 		# Make it a real file, not a temp file
	my $unique_seqs_fasta = $save . $filename . ".uniq" . $ext; 		# Make it a real file, not a temp file
	
#	my $check_for_fasta_collapser = which("fasta_collapser.pl");		# returns undef if not found on system.
#	if ($check_for_fasta_collapser){
		# Found fasta_collapser.pl
		print STDERR "Saving unique sequences to $unique_seqs_fasta\n";
		my $cmd = $PWD."/fasta_collapser.pl -i $file -o $unique_seqs_fasta";
		system($cmd);
		return $unique_seqs_fasta;
	# }
	# else {
	# 	print STDERR "Can't find fasta_collapser.pl to make unique sequences.  Will run each sequence in the input file. (I suggest you cancel and install fastx toolkit to make this step faster.)\n";
	# 	return $file;
	# }
}	
#-----------------------------------------------------------------------------
sub read_input_reads {
	my $file = shift;
	my $nuc_aa_codon_tally; 
	print STDERR "Processing file $file ...\n";
	
	
	my $in = Bio::SeqIO->new(-file => $file , -format => 'Fasta' );		
	
	
	my $pl = Parallel::Loops->new($cpu);

	my $count = 0;		
	# Beginning of parallel loop.  See http://search.cpan.org/~pmorch/Parallel-Loops-0.03/lib/Parallel/Loops.pm
	# Syntax:   $pl->while($conditionSub, $childBodySub)
	# Conceptually similar to 
	#	while($conditionSub->()) {
	#		$childBodySub->();
	#	}
	my $seq = "";
	
	# Set up shared variables to store data.
	my (@good,@toss);	 # arrays to save information for the reads that were kept or tossed. Each element of either array is a hashref with information about the read, including sequence, location of early stop codons, etc.  
	$pl->share(\@good, \@toss);
	
	$pl->while( sub { $count++; $seq = $in->next_seq();  },	sub {			# was $pl->foreach( \@primerIDs,	sub {	

		# Report progress
		my $denominator = 1000;
		$denominator = 10 if ($verbose);
		if ($count % $denominator == 0){
			print STDERR "Processed $count reads. "; 
			&elapsed($start_time, 'Elapsed', $verbose);
		}
		
		# Create temporary fasta file
		my $cwd = Cwd::cwd();
		my $tempdir = File::Temp->newdir( "/tmp/temp_dir_seq_".$count."_XXXXXXXX" );		# CLEANUP => 0 
		my $temp_fasta = $tempdir."/temp.fa";
		my $temp_aln 	= $tempdir."/temp.aln";
#		my $temp_aln 	= "temp.aln";
		my $fh = open_to_write($temp_fasta, 0, 0, 1);
		print $fh ">ref\n$refseq\n>read ". $seq->display_id(). "\n".$seq->seq()."\n";
		close($fh);
	
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
								if ($debug){	print STDERR "Running command: $cmd\n";		} 
			system($cmd);
			my $aln_in = Bio::AlignIO->new(
				-file 	=> $temp_aln,
				-format	=> 'fasta',
			);
			$aln = $aln_in->next_aln();
		}
				
								if ($debug){	my $seqIO = Bio::AlignIO->new(-fh => \*STDERR, -format => 'clustalw');  print STDERR "ALIGNMENT: \n"; $seqIO->write_aln($aln);	}
		
		# Get the aligned reference and read strings (with gaps)
		my $ref_aligned  = $aln->get_seq_by_id("ref");
		my $read_aligned = $aln->get_seq_by_id("read");

		# Get a matrix to represent gaps in the alignment. 
		my $gap_char = "-";		# dash as default, for mafft alignment output
		$gap_char = "." if ($clustalw);	# gaps are represented by '.' in clustalw alignment output.  
#		my $mat = $aln->gap_col_matrix($gap_char);		# This makes an AoH where each position of the array represents a position in the alignment and the keys of the hash are sequence ids, with value of 1 if there is gap at this position or '' if no gap.  (Not using...)
#								if ($verbose){	print STDERR "gap mat:\n";		print STDERR Dumper($mat);	}

		# Walk through the alignment column by column, starting at the first base of the read in the alignment
#		my $len = $aln->length();
		my $start = $aln->column_from_residue_number("read", 1);	# Get position of the beginning of the read in the alignment.
		my $read_length = length($seq->seq() );
		my $end   = $aln->column_from_residue_number("read", $read_length);	# Get position of the end of the read in the alignment.
		my $ref_nuc_pos = $start - 1;		# Count of reference bases considered so far (including current position). (For getting $ref_frame and $read_nuc_pos).  Was $ref_base_count.  $read_nuc_pos is the position to use to store the bases in the hash.
		
		my $first_codon = 0; 	# Will become 1 when we've entered the initial first position within a codon (i.e., first full codon; to avoid translating 5' partial codons in the read)
		my $frameshift = 0;		# Will become read gap size if there is a gap (unless resolved).  If gap in read, add to $frameshift.  If gap in ref at this position, deduct from $frameshift.  If ($frameshift % 3) == 0, then in frame, otherwise out of frame.  
		my $out_of_frame = 0; 	# Will become 1 when frameshift initially occurs.  Will return to 0 when frameshift resolves 
		my @frameshift_positions;	# Array of arrays to store the positions where frameshift occurs (for when we hit a stop, to know when the last frameshift began).  Each element is an array [$position, type] where type is 'insertion_in_read' or 'deletion_in_read'.  Each entry represents the initial position that the read sequence goes out of frame; e.g., if there is a 2bp insertion, there will only be an entry for the first base position, but not for the second base.  Position is $ref_nuc_pos, not read_nuc_pos.  (should it be read_nuc_pos?)
		my $ref_codon;		# Used to get reference codons for translation
		my $read_codon;		# will be used to get codon as three nucleotides become available.  flushed after translating and storing.  Only add to it once partial_codon_bases_left is at zero.
#		my ($ref_gaps,$read_gaps) = (0,0);	# Count the number of gaps in the aligned region (in either the read or the reference).  Used to get the nuc_pos with respect to the reference sequence 
		
		my $seq_info; # HoH hashref with first level keys 'nuc', 'codon', 'aa', 'frameshift_positions', 'type', 'full_seqs'; second level keys are positions with value as the base, codon, or aa for the read.  Value for key frameshift_positions is a reference to @frameshift_positions.  Value for key 'type' is a reference to \%type.  Value for 'full_seqs' is hash with keys 'nuc' and 'aa'
		my %type; 	# Codon types.  keys 'normal', 'in-frame_indel', 'out-of-frame_indel', 'unknown'
		my $toss = 0;		# Will become 1 if the read has a premature stop codon due to a frameshift.  Then the local hash $seq_info will be pushed to @toss instead of @good.
		
		my ($full_read_nuc,$full_read_aa);
		
		NUC: for (my $aln_pos = $start; $aln_pos <= $end; $aln_pos++){
			my $ref_res  = $ref_aligned->subseq($aln_pos,$aln_pos); 
			my $read_res = $read_aligned->subseq($aln_pos,$aln_pos);
			
			# Check for gap in read or reference  (can't have gap in both read and ref, so if gap in one then the other is a base and can be saved)
			if ($read_res =~ m/$gap_char/){
				$frameshift++;
#				$read_gaps++;
				
			}
			elsif($ref_res =~ m/$gap_char/){
				$frameshift--;
#				$ref_gaps++;
			}

			# Get nucleotide positions for this base
			$ref_nuc_pos++ unless ($ref_res =~ m/$gap_char/);		# Is this always the same value as $ref_aligned->location_from_column($aln_pos);
			my $read_nuc_pos = $ref_nuc_pos - $frameshift;			# This is the position for storing the nucleotide in the local hash.  **Is this working?  In all situations, including frameshifts that are resolved a few codons later?
						
			# Check if we have reached the first full codon yet or if we are still in an initial partial codon
			my $ref_frame = $ref_nuc_pos % 3;	 # Translation frame of reference sequence for the "current" reference base.  1 if first base in a codon, 2 if second base in a codon, 0 if in third base in a codon.  
			unless($first_codon){
				$first_codon++ if ($ref_frame == 1);		# Must be the first position in a codon to start adding bases to codons.  
			}

			
			next NUC unless ($first_codon);




			# We have reached/passed the first full codon.  Add base(s) to read_codon string (not gap characters).  Also start making the full nucleotide sequence for the cleanreads.txt file (**calculate_linkage_disequilibrium script assumes the first triplet is in frame**, e.g., @seq = unpack("(a3)*", $read_seq);	# codon, split into triplets).  Therefore, don't start saving nucleotide positions in seq_info->{nuc} until here either. 

			# Save nucleotide information in local hash
			#*** check the format we need for individual and merged tables.  In report, I need a hashref of this format: $nuc_aa_codon_tally->{'nuc'}->{$nuc_pos}->{$nuc} = count;  For now, I'll make an array of these: {'nuc'/'aa'/'codon'}->{$pos}->{$nuc/$aa/$codon} = count
			$seq_info->{'nuc'}->{$read_nuc_pos} = uc($read_res); 

			$full_read_nuc .= $read_res;
			$read_codon = $read_codon . $read_res 	 unless ($read_res =~ m/$gap_char/);
			
			# Check that we have a full codon in read to translate
			next NUC unless (length($read_codon) == 3);
			
			
			# There are three bases in the read_codon if we reach this point
			
			# Save position if a frameshift is found.  (Do this only after we have reached a full codon to first get IN frame.)  [Is this working?]
			unless ($out_of_frame){		# If $out_of_frame is already set, then don't save the frameshift position again
				if ($ref_frame != 0){	 
					$out_of_frame++;
								if ($verbose){	print STDERR "read_nuc_pos: $read_nuc_pos frameshift: $frameshift out_of_frame: $out_of_frame ref_nuc_pos: $ref_nuc_pos ref_frame: $ref_frame\n"; 	}
					my $indel = "deletion_in_read";		# default, deletion in read  
					if ($ref_frame < 0){
						$indel = "insertion_in_read";
					}
					my @frameshift = ($ref_nuc_pos, $indel);
					push @frameshift_positions, \@frameshift; 	
				}				
			}
			if ($ref_frame == 0){		# If not out of frame, make sure that $out_of_frame is not set.  
				$out_of_frame = 0;
			} 

			# Now translate codon 
			my $read_aa_pos = $read_nuc_pos / 3;		# Does this work?  Is it always a multiple of three?   
			$read_codon = uc($read_codon);
			my $read_aa = ""; 
			if (exists($converter{$read_codon})){
				$read_aa = $converter{$read_codon};
			}
			else {
				$read_aa = "X";
			}
			
			$full_read_aa .= $read_aa;
			
			# Save the codon type at this position.  (Not sure how useful this will be...) 
				# if ($frameshift == 0), then they are an in-frame codons, so translate read_codon and store it in the local hash.  
				# if (($frameshift != 0) && ( ($frameshift % 3) = 0) ), then it is an in frame insertion/deletion, which is fine 
				# if ( ($frameshift % 3) == 0), then it is a frameshifted codon.  Translate it and flag the read; if we reach a stop codon further down, then mark the read to be tossed.  
					# If it's a STOP (*) then mark the read so the local hash will be added to the the global array @toss instead of @good.  Save the position of the last 
			if ($frameshift == 0){		# In-frame codons with no indel present.  Most often, we'll be here.  
				$type{'normal'}++;
			}
			elsif(($frameshift != 0) && ( ($frameshift % 3) == 0) ){		# In-frame codons with an insertion/deletion present (i.e., length of insertion or deletion is a multiple of 3)
				$type{'in-frame_indel'}++;
			}
			elsif( ($frameshift % 3) != 0){		# Frameshift codon.  Translate it 
				$type{'out-of-frame_indel'}++;
			}
			else {
				$type{'unknown'}++;
							if ($verbose){	print STDERR "out_of_frame: $out_of_frame read_nuc_pos: $read_nuc_pos read_aa_pos: $read_aa_pos read_codon: $read_codon frameshift: $frameshift\tmodulus: ". $frameshift % 3 . "\n";	}
				# Shouldn't be able to arrive here...
			}
			
			# Check for premature stop codons due to frameshift
			if ( ($read_aa =~ m/^\*$/) && $out_of_frame ){		# Stop codon found when frameshift is present.  Mark the read to be tossed
				$toss = 1;
							if ($verbose){	print STDERR "Premature stop in a frameshift region!\nFull nuc seq (so far): $full_read_nuc\nFull aa seq (so far): $full_read_aa\ncodon:$read_codon\n";		}
			}
			
			# Save the codon and amino acid sequences to local hash
			$seq_info->{'codon'}->{$read_aa_pos} = $read_codon;
			$seq_info->{'aa'}->{$read_aa_pos} = $read_aa;
			

			# Flush codon sequence
			$read_codon = "";

			last NUC if ($toss);		# No need to continue if we are going to toss this.  

		
		}	
		
		# Save @frameshift_positions to local hash
		$seq_info->{'frameshift_positions'} = \@frameshift_positions;
		
		# Save %type to local hash.  TESTING
		$seq_info->{'type'} = \%type;
		
		# Save the full sequences to local hash
		$seq_info->{'full_seqs'}->{'nuc'} = uc($full_read_nuc);
		$seq_info->{'full_seqs'}->{'aa'} = uc($full_read_aa);
		$seq_info->{'full_seqs'}->{'id'} = $read_aligned->desc();
		
		# Save the seqID to local hash
		$seq_info->{'id'} = $seq->display_id();		# Important for when processing unique sequences.  e.g., 2-520 means sequence 2, count = 520.  		
		
		# Save $seq_info local hash to shared array @good or @toss
		if ($toss){
			push @toss, $seq_info;
		}
		else {
			push @good, $seq_info;
		}
			
	});		### End of Parallel Loop

	# Count good and tossed sequence
	my ($good_unique_count,$good_total_count) = count_total_from_good_toss_arrays(\@good);
	my ($toss_unique_count,$toss_total_count) = count_total_from_good_toss_arrays(\@toss);
	print STDERR "Good reads:\t$good_unique_count unique ($good_total_count total)\nTossed reads:\t$toss_unique_count unique ($toss_total_count total)\n";

	# Now walk through the sequences in @good
		# tally up nucleotides, codons, and amino acids. [Should this be in a separate script?  It could just read the output aligned cleaned reads.]
			# Save in $nuc_aa_codon_tally with first keys 'nuc', 'codon', or 'aa'; second, position; third, base/aa/codon -> count
		# print out cleanreads and cleanpeptides
		# Save frameshift positions and types?  Or at least tally them up? (if it's useful)
	
	print STDERR "Processing good reads to make cleanread and cleanpeptide PHYLIP files. " if ($good_unique_count);
	&elapsed($start_time, 'Elapsed', $verbose);
	# Prepare the output files of cleaned reads and full peptides
	my ($filename,$dir,$ext) = fileparse($file,@SUFFIXES);
	
	my $newfilename = $save_dir . "/" . $filename . ".cleanreads.txt";
	my $clean_reads_fh = open_to_write("$newfilename");
	print $clean_reads_fh "#readID|gene|positionForFirstNucleotideBase\tcleanRead\n";
	
	my $fullpeptide = $save_dir . "/" . $filename . ".cleanpeptides.txt";
#					print STDERR "Full peptide alignment will be stored in this file: $fullpeptide\n";
	my $pepfh = open_to_write("$fullpeptide");
	print $pepfh "#readID|gene|positionForFirstAminoAcid\tcleanPeptide\n";

#	my $printed; 	#Hashref to keep track of the sequences that have been printed already.  
	my ($maxlength,$majority_length,$peplength) = get_max_sequence_length_from_fasta($file);
	
	foreach my $seq_info (@good){
							if ($debug){	print Dumper($seq_info->{full_seqs});	}
							if ($debug){	print Dumper($seq_info);				}
		
		my $seq_count = get_count_from_id( $seq_info->{'id'} );
		
		# Save tally of codons, aa, nuc
		foreach my $type (qw(codon aa nuc) ){
			my @pos = sort {$a <=> $b} keys (%{$seq_info->{$type}});
			my $min_pos = range(\@pos,'min');
			foreach my $pos (@pos){
				my $base_codon_aa = $seq_info->{$type}->{$pos};
				$nuc_aa_codon_tally->{$type}->{$pos}->{$base_codon_aa} += $seq_count;
			}
			# Print out the full sequences		[Are these printing out properly?  What about sequences with insertions in the read?  Should they be removed?]
			# Do I need to print out the individual reads that are identical?  If so, do I need the original primerid in the name?  
			# Is the fact that I have ids such as 1-520 and 2-32 -- i.e., different lengths -- messing up the PHYLIP format?  Before I had the primerid, which was all the same length.  http://www.phylo.org/index.php/help/relaxed_phylip 
			if ($type eq 'aa'){
				print $pepfh 			$seq_info->{'full_seqs'}->{id}	. "|$gene|$min_pos " . $seq_info->{'full_seqs'}->{aa} . "\n";
			}
			elsif ($type eq 'nuc'){
				print $clean_reads_fh	$seq_info->{'full_seqs'}->{id}	. "|$gene|$min_pos " . $seq_info->{'full_seqs'}->{nuc} . "\n";
			}
		}
	}
	close($clean_reads_fh);
	close($pepfh);


	# Prepend the PHYLIP header to the cleanreads.txt and cleanpeptide.txt files
	foreach my $fix_file ($newfilename, $fullpeptide){
		my ($count, $aln_length) = get_count_and_aln_length($fix_file);
								if ($debug){	print STDERR "file: $fix_file count: $count aln_length: $aln_length\n";	}
		prepend_file( "$fix_file", "  $count $aln_length\n" );
	}

		
	# Walk through the sequences in @toss, tally up last frameshift positions (for deletions and insertions separately)
	my $indels;	 # hashref to tally the positions of insertions and deletions.  first keys deletion_in_read, insertion_in_read; second keys positions, value = count
	foreach my $seq_info (@toss){
		#print Dumper($seq_info->{frameshift_positions}); exit;
		my $seq_count = get_count_from_id( $seq_info->{'id'} );
		foreach my $frameshift (@{$seq_info->{frameshift_positions}}){
			my ($pos,$indel) = @$frameshift;
			$indels->{$indel}->{$pos} += $seq_count;
		}
	}
	print STDERR "Location of indels in tossed reads:\n" if ($toss_unique_count);
	foreach my $indel (keys %$indels){
		print STDERR "Ref_Nuc_Pos\t$indel\n";
		foreach my $pos (sort {$a <=> $b} keys %{$indels->{$indel}}){
			print "$pos\t$indels->{$indel}->{$pos}\n";
		}
	}
	

	
	
	return $nuc_aa_codon_tally;
	
}
#-----------------------------------------------------------------------------
sub count_total_from_good_toss_arrays {
	my ($array,$good_or_toss) = @_;
	# Array is an AoH where a pertinent key in the has is the sequence id, which has format such as  2-520, meaning sequence 2, count = 520.
	my ($unique_count,$total_count);
	foreach my $seq (@$array){
		$unique_count++;
		my $id = $seq->{'id'};
		my ($seq_num,$count) = split(/-/, $id);
		$total_count += $count;
	}
	return ($unique_count,$total_count);
}
#-----------------------------------------------------------------------------
sub get_count_from_id {
	# Get the sequence id for counting/amplification purposes.  
	# For example, if the id says 2-520, that means sequence 2, count = 520.  (This is what fasta_collapser.pl output ids look like).  
	# If the id doesn't match this format, then the count will remain as 1.
	my $id = shift;		# From $seq_info->{'id'}
	my ($num,$seq_count) = (1,1);
	if ($id =~ m/^\d+-\d+$/){
		($num,$seq_count) = split(/-/, $id);
	}
	return $seq_count;
}
#-----------------------------------------------------------------------------
sub get_max_sequence_length_from_bam {
	my $file = shift;
	my $pep_length = 0;
	my $max_length = 0;
	my $majority_length = 0;
	
	# Read through the first BAM file (up to 50K alignments) to get the maximum length of a read for the peptide and clean_reads alignment
	# Also get the majority read length (to use in the alignment)
	my $bam = Bio::DB::Bam->open("$file");
	my $header = $bam->header;
	my $n = 0;
	my %lengths;
	BAM: while(my $align = $bam->read1){
		my $length = $align->calend - $align->pos;
		$max_length = $length if ($length > $max_length);
		$lengths{$length}++;
		$n++;
		last BAM if ($n > 50000);
	}
#	print Dumper(\%lengths);
	$majority_length = find_key_with_biggest_value(\%lengths);
#	$pep_length = int(($max_length + 2)/3);
	$pep_length = int(($majority_length + 2)/3);
					#print STDERR "$n max: $max_length\tpeplength: $peplength\n"; exit;

	return ($max_length,$majority_length,$pep_length);
}
#-----------------------------------------------------------------------------
sub get_max_sequence_length_from_fasta {
	my $file = shift;
	my $pep_length = 0;
	my $max_length = 0;
	my $majority_length = 0;
	
	# Read through the first BAM file (up to 50K alignments) to get the maximum length of a read for the peptide and clean_reads alignment
	# Also get the majority read length (to use in the alignment)
	my $in = Bio::SeqIO->new(-file => $file , -format => 'Fasta' );
	my $n = 0;
	my %lengths;
	FASTA: while(my $seq = $in->next_seq() ){
		my $length = length($seq->seq() );
		$max_length = $length if ($length > $max_length);
		$lengths{$length}++;
		$n++;
		last FASTA if ($n > 50000);
	}
#	print Dumper(\%lengths);
	$majority_length = find_key_with_biggest_value(\%lengths);
#	$pep_length = int(($max_length + 2)/3);
	$pep_length = int(($majority_length + 2)/3);
					#print STDERR "$n max: $max_length\tpeplength: $peplength\n"; exit;

	return ($max_length,$majority_length,$pep_length);
}
#-----------------------------------------------------------------------------
sub get_count_and_aln_length {
	my $file = shift;
	my $count = 0;
	my $aln_length = 0;
	my $readfh = open_to_read($file);
	while(<$readfh>){
		next if ( (m/^#/) || (m/^\s+$/) ); 	# skip header lines and empty lines
		unless ($count){
			# First line
			my @line = split(/\s+/);
			$aln_length = length($line[1]);
		} 
		$count++;
	}
	
	close($readfh);
	return ($count, $aln_length);
}
#-----------------------------------------------------------------------------
sub parse_gff {
	# Takes a gff file and a genome multi-fasta file of chromosomes and outputs a hashref containing the CDS starts, stops, strand, nucleotide sequence and protein sequence for each gene.  Also, nucleotide, amino acid and codons tables for lookup by position for each gene.  (Is it possible to have multiple CDS per gene in a gff file?  How would that affect this?)
	my $file = shift;
	my $genes;	# hashref to store the coordinates for the CDS of genes on the chromosomes.  first keys, chr; second keys, gene id; third keys start, stop, starts, stops, strand, nucleotide, protein, nuc->pos, aa->pos, codon->pos.
	my $readfh = open_to_read($file);
	
	while(<$readfh>){
		next if (m/^#/);
		my @line = split(/\t/);
		next unless (scalar(@line) > 8);
		next unless ($line[2] eq "CDS");
		my ($chr, $start, $stop, $strand) = ($line[0], $line[3], $line[4], $line[6]);
		my $gene_id = get_gff_attribute($line[8], "ID");
		my $codon_start = get_gff_attribute($line[8], "codon_start");
					if ($debug){	print "id: $gene_id\tcodon_start: $codon_start\n";	}
		if ($codon_start && ($codon_start > 1)){
			unless ($genes->{$chr}->{$gene_id}->{starts}){		# This only applies to the first exon.  For some reason, the codon_start attribute occurs on all CDS exon features of a CDS, all referring to the start codon of the first CDS exon.
				# Then this is the first CDS exon.
				$start += $codon_start - 1;
			}
		}
		
#		next unless ($strand eq "+");		# Assuming + strand for now because my annotation is all on + strand.		

		$genes->{$chr}->{$gene_id}->{strand} = $strand;	
		if (exists($genes->{$chr}->{$gene_id}->{start})){
			# One portion of the CDS is already present.  For the purpose of the TSS and TTS, get the earliest start and latest stop
			$genes->{$chr}->{$gene_id}->{start} 	= ($start < $genes->{$chr}->{$gene_id}->{start})	? $start 	: $genes->{$chr}->{$gene_id}->{start}; 		# If the value present is larger, then reassign to the smaller one.
			$genes->{$chr}->{$gene_id}->{stop} 	= ($stop  > $genes->{$chr}->{$gene_id}->{stop}) 	? $stop 	: $genes->{$chr}->{$gene_id}->{stop}; 		# If the value present is smaller, then reassign to the larger one.
		}
		else {
			$genes->{$chr}->{$gene_id}->{start} = $start;
			$genes->{$chr}->{$gene_id}->{stop} = $stop;				
		}
		
		push @{ $genes->{$chr}->{$gene_id}->{starts} }, $start;				# Some influenza proteins are spliced.  In that case, multiple CDS lines from GFF file need to be concatenated...
		push @{ $genes->{$chr}->{$gene_id}->{stops} }, $stop;	
	
	}
	
	close ($readfh);

#				return $genes;

	# Now take CDS entries and get sequence for them from --fasta input file, then translate that sequence
	my $in  = Bio::SeqIO->new(-file => "$ref" ,
                           -format => 'Fasta');
	while ( my $seq = $in->next_seq() ) {
		my $chr = $seq->display_id();	# chr name		#print "id: $chr\n";
		next unless exists ($genes->{$chr});		# Skip this fasta record if there are no annotated genes on this chromosome.
		foreach my $gene (keys %{$genes->{$chr}}){		# Look at each gene on this chromosome.  Assuming + strand since my GFF file only has + strand genes.  Actually, I can just reverse complement the final concatenated sequence if it is the - strand.  
			# Get the first CDS segment
			my $seqstr_nuc = $seq->subseq( $genes->{$chr}->{$gene}->{starts}->[0], $genes->{$chr}->{$gene}->{stops}->[0] );		#	print "$chr $gene: ".length($seqstr_nuc)."\n";
			if (scalar(@{$genes->{$chr}->{$gene}->{starts}})>1){		# If there are additional CDS segments
				# Then get the rest of the sequences.	I can concatenate the rest to the first segment $seqstr_nuc because I've assumed above that they are all on the + strand.
				for (my $i = 1; $i < @{$genes->{$chr}->{$gene}->{starts}}; $i++){
					$seqstr_nuc .= $seq->subseq( $genes->{$chr}->{$gene}->{starts}->[$i], $genes->{$chr}->{$gene}->{stops}->[$i] );	#	print "$chr $gene: ".length($seqstr_nuc)."\n";
				}
			}
			if ($genes->{$chr}->{$gene}->{strand} eq '-'){
				my $seqstr_nuc = Fastq_utils::rev_comp($seqstr_nuc);	# Not tested.  
			}
			my $prot = translate_seq($seqstr_nuc);
			$genes->{$chr}->{$gene}->{nucleotide} = $seqstr_nuc;
			$genes->{$chr}->{$gene}->{protein} = $prot;
			
			# Now save the nucleotides, amino acids and codons by position.  
			my @nuc = split(/|/, $seqstr_nuc);
			my @aa = split(/|/, $prot);
			my @codon = unpack("(a3)*", $seqstr_nuc);		# http://stackoverflow.com/questions/372370/how-can-i-split-a-string-into-chunks-of-two-characters-each-in-perl
				# Nucleotides first
			for (my $i = 0; $i < @nuc; $i++){
				my $nt = uc($nuc[$i]);
				my $pos = $i + 1;
				if (exists($genes->{$chr}->{$gene}->{'nuc'}->{$pos})){
					warn "This gene + nucleotide position already present: $gene, $pos, $nt\n";
				}
				else {
					$genes->{$chr}->{$gene}->{'nuc'}->{$pos} = $nt;
				}
			}
			
				# Now save amino acids and codons
			for (my $i = 0; $i < @aa; $i++){
				my $aa = uc($aa[$i]);
				my $codon = uc($codon[$i]);
				my $pos = $i + 1;
#				unless ($aa =~ m/X/i){			# Commented this out because if we skip these amino acids, it will cause a frameshift.
					# Amino acid
					if (exists($genes->{$chr}->{$gene}->{'aa'}->{$pos})){
						warn "This gene + amino acid position already present: $gene, $pos, $aa\n";
					}
					else {
						$genes->{$chr}->{$gene}->{'aa'}->{$pos} = $aa;
					}
					# Codon
					if (exists($genes->{$chr}->{$gene}->{'codon'}->{$pos})){
						warn "This gene + codon position already present: $gene, $pos, $codon\n";
					}
					else {
						$genes->{$chr}->{$gene}->{'codon'}->{$pos} = $codon;
					}
#				}
			}	
		}
	}
	
	return $genes;
}
#-----------------------------------------------------------------------------
sub get_gff_attribute {
	# gets the gene ID.  "IDs must be unique within the scope of the GFF file"	http://gmod.org/wiki/GFF
	my $attributes = shift;
	my $tag = shift;
	my @attributes = split(/;/, $attributes);
	my $value;
	foreach my $pair (@attributes){
		my @pair = split(/=/, $pair);
		if ($pair[0] eq "$tag"){
			$value = $pair[1];
		}
		else {
#			print "Not ID: $pair[0]\n";
		}
	}
	unless ($value){
		die "No $tag found in the GFF attributes: \n$attributes\n";
	}
	return $value;
}
#-----------------------------------------------------------------------------
sub get_cds_sequences {
	# Takes a nucleotide sequence representing the coding region and stores the nucleotides, codons, and amino acids by position.  
	my $seq = shift;
	my $cds; 	# hashref to store sequences by position.  First keys 'nuc' 'codon' 'aa' (type); second keys position; value is the base/aa/codon
	
	my $num_codons = length ($seq) / 3; 
	
	my @seq = split(/|/, $seq);
	NUC: for (my $i = 0; $i < @seq; $i++){
		my $nuc_pos = $i + 1;
		$cds->{'nuc'}->{$nuc_pos} = $seq[$i];
	}
	
	CODON: for (my $i = 0; $i < $num_codons; $i++){		 # Could use unpack as in get_reference_protein_hash if desired.  
		my $codon = substr($seq, 0, 3); 
		my $aa_pos = $i + 1;
		next CODON if (length($codon) < 3);		# partial codons at the 3' end of the sequence  
		my $aa = "";
		if (exists($converter{uc $codon})){
			$aa = $converter{uc $codon};
		}
		else {
			 warn "Reference codon position $i: no codon found $codon ! Assigning X as amino acid.\n";	
			 $aa = 'X';		# Symbol X is used for amino acids that are unidentified.
		}
		$cds->{'codon'}->{$aa_pos} 	= $codon;
		$cds->{'aa'}->{$aa_pos} 	= $aa;

		$seq = substr($seq, 3);		# delete the codon we just looked at.
	}
		
	return $cds;
}
#-----------------------------------------------------------------------------
sub translate_seq {
	# Assumes the sequence is in frame.  Translates to amino acid sequence
	my $seq = shift;	#nucleotide sequence
	my @aa;	
	my $num_codons = length ($seq) / 3; 
	CODON: for (my $i = 0; $i < $num_codons; $i++){		 # Could use unpack as in get_reference_protein_hash if desired.  
		my $codon = substr($seq, 0, 3); 
		next CODON if (length($codon) < 3);		# partial codons at the 3' end of the sequence  
		my $aa = "";
		if ($codon =~ m/X/i){		# Codons with a deletion in the read compared to reference; give them aa '*', like the symbol in mpileup.  Some codons with N will work in the converter hash (e.g., wobble third position), so not skipping them here.
			$aa = '*';
		}
		elsif (exists($converter{uc $codon})){
			$aa = $converter{uc $codon};
		}
		else {
			 if ($verbose){	warn "no codon found $codon ($i) ! Assigning X as amino acid.\n";	}
			 $aa = 'X';		# Symbol X is used for amino acids that are unidentified.
		}
		push(@aa, $aa);
		$seq = substr($seq, 3);		# delete the codon we just looked at.
	}
	my $prot = join "", @aa;	# full protein / amino acid sequence
	return $prot;
}
#-----------------------------------------------------------------------------
sub print_reports {
	my ($nuc_aa_codon_tally, $label) = @_;
	# Print Reports to filehandles $nucleotide_report_fh, $codon_report_fh, and $amino_acid_report_fh
	#name(chr:gene:aminoAcidPosition:consensusAminoAcid)	chr	gene	aminoAcidPosition	refAminoAcid	consensusAminoAcid	Sample	coverageDepth	numA	numC	numG	numT	numN	[del	SumACGT	numNonzeroACGT	Ref	Nonref]

	my @all_codons = @CODON;
	for (@all_codons){
		s/:\w+//;
	}
		
	
	foreach my $type (qw(nuc codon aa)){		#e.g., $nuc_aa_codon_tally->{'aa'}->{$pos}->{$aa}++;
		foreach my $pos (sort {$a <=> $b} keys %{$nuc_aa_codon_tally->{$type}}){
			my $ref_sequence = $cds->{$type}->{$pos};
			unless ($ref_sequence){
				print STDERR "Error in print_reports: no reference $type: $gene $pos\n";	
				if ($debug){	print Dumper($cds->{$type});		}
				if ($debug){	print Dumper($nuc_aa_codon_tally->{$type}->{$pos});	}
			}
			my $consensus = find_key_with_biggest_value($nuc_aa_codon_tally->{$type}->{$pos});		# Assuming there won't be a tie for most abundant
			my $diff = "";
			$diff = "DIFF" if ($consensus ne $ref_sequence);
			my $name = $gene.":".$pos.":".$consensus;
			my $coverage = total(values(%{$nuc_aa_codon_tally->{$type}->{$pos}})); 

			# Now print the counts, number of amino acids represented, ref aa counts, and nonref aa counts.
			if ($type eq 'nuc'){
				print $nucleotide_report_fh join "\t", $name, $gene, $pos, $ref_sequence, $consensus, $diff, $label, $coverage;
				print_tally_line($nucleotide_report_fh, 	$nuc_aa_codon_tally->{$type}->{$pos}, \@NUC);
			}
			elsif($type eq 'codon'){
				print $codon_report_fh join "\t", $name, $gene, $pos, $ref_sequence, $consensus, $diff, $label, $coverage;
				print_tally_line($codon_report_fh, 		$nuc_aa_codon_tally->{$type}->{$pos}, \@all_codons);
			}
			else {
				print $amino_acid_report_fh join "\t", $name, $gene, $pos, $ref_sequence, $consensus, $diff, $label, $coverage;
				print_tally_line($amino_acid_report_fh, 	$nuc_aa_codon_tally->{$type}->{$pos}, \@AA);
			}
		}
	}
}
#-----------------------------------------------------------------------------
sub print_tally_line {
	my ($writefh, $tally, $all_array) = @_;
	my @seen;
	foreach my $nuc_codon_aa (@$all_array){
		if (exists($tally->{$nuc_codon_aa})){
			my $count = $tally->{$nuc_codon_aa};
			print $writefh "\t".$count;
			push @seen, $nuc_codon_aa;
		}
		else {
			print $writefh "\t0";
		}
	}
	my @keys = keys %$tally;
	my @other = array_minus( @keys, @seen );		#http://stackoverflow.com/questions/2933347/comparing-two-arrays-using-perl  # get items from array @a that are not in array @b: my @minus = array_minus( @a, @b );
	if (@other > 0){
		my $other_sum = 0;
		foreach (@other){
			$other_sum += $tally->{$_};
		}
		print $writefh "\t$other_sum";		# Maybe I should find a way to record what these other ones are... a separate report somewhere?  probably a lot of Xs, is my guess.  
		my $others;
		foreach (@other){
			$others .= "$_:$tally->{$_},";

		}
		$others =~ s/,$//;
		print $writefh "($others)";

	}
	else {
		print $writefh "\t0";
	}
	print $writefh "\n";
}
#-----------------------------------------------------------------------------
sub print_merged_report {
	my ($nuc_aa_codon_tally, $label) = @_;
	# Print Reports to filehandles $nucleotide_report_fh, $codon_report_fh, and $amino_acid_report_fh
	#name(gene:aminoAcidPosition:consensusAminoAcid)	gene nucleotidePosition refNucleotide consensusNucleotide aminoAcidPosition refCodon	consensusCodon	refAminoAcid	consensusAminoAcid	Sample	coverageDepth	topNucleotide	secondNucleotide	thirdNucleotide	otherNucleotide	countTopNucleotide	countSecondNucleotide	countThirdNucleotide	countOtherNucleotide	topCodon	secondCodon	thirdCodon	otherCodon	countTopCodon	countSecondCodon	countThirdCodon	countOtherCodon	topAminoAcid	secondAminoAcid	thirdAminoAcid	otherAminoAcid	countTopAminoAcid	countSecondAminoAcid	countThirdAminoAcid	countOtherAminoAcid
	# $nuc_aa_codon_tally format example: $nuc_aa_codon_tally->{'aa'}->{$pos}->{$aa}++;
	my @all_codons = @CODON;
	for (@all_codons){
		s/:\w+//;
	}
		
	
	my $first_nuc_pos = (sort {$a <=> $b} keys %{ $nuc_aa_codon_tally->{'nuc'} })[0]; 
#						print "first nuc position: $first_nuc_pos\n"; exit; 
	my $first_aa_pos =  (sort {$a <=> $b} keys %{ $nuc_aa_codon_tally->{'aa'} })[0]; 
#						print "first aa position: $first_aa_pos\n"; exit;
	my $last_aa_codon_line = "\t"x16 ."\n";			# For repeating the amino acid and codon information for each nucleotide base of the codon.  Default will be an empty line (keep tabs for proper spreadsheet spacing) in case there are a few nt before the first codon.  
	foreach my $nuc_pos (sort {$a <=> $b} keys %{$nuc_aa_codon_tally->{'nuc'}}){		# 
		my $aa_pos = int(($nuc_pos+2)/3);		# Need to add 2 first before dividing by 3.  E.g., ACGTAG, for nucl position 1, to get amino acid position, int((1+2)/3) = 1.  For nuc position 3, int((3+2)/3)=1. For nuc position 4, int((4+2)/3) = 2.
		my $position_in_codon = get_position_in_codon($nuc_pos);
		my $ref_nuc 		= $cds->{'nuc'}->{$nuc_pos};
		my $consensus_nuc 	= find_key_with_biggest_value($nuc_aa_codon_tally->{'nuc'}->{$nuc_pos});		# Assuming there won't be a tie for most abundant
		my $diff = "";
		$diff = "DIFF" if ($consensus_nuc ne $ref_nuc);
		my $ref_codon 		= $cds->{'codon'}->{$aa_pos};
		my $consensus_codon	= find_key_with_biggest_value($nuc_aa_codon_tally->{'codon'}->{$aa_pos}) || ""; 
		my $ref_aa			= $cds->{'aa'}->{$aa_pos};
		my $consensus_aa	= find_key_with_biggest_value($nuc_aa_codon_tally->{'aa'}->{$aa_pos}) || "";	
		my $name = $gene.":".$aa_pos.":".$consensus_aa;
		my $coverage = total(values(%{$nuc_aa_codon_tally->{'nuc'}->{$nuc_pos}})); 

		# Print reference and consensus information
		print $merged_report_fh join "\t", $name, $gene, $nuc_pos, $position_in_codon, $ref_nuc, $consensus_nuc, $diff, $aa_pos, $ref_codon, $consensus_codon, $ref_aa, $consensus_aa, $label, $coverage;

		# Print Nucleotide top hits, then Codon, then Amino Acid, e.g., topNucleotide	secondNucleotide	thirdNucleotide	otherNucleotide	countTopNucleotide	countSecondNucleotide	countThirdNucleotide	countOtherNucleotide	
		my $arrayRef = [ ['nuc', $nuc_pos], ['codon', $aa_pos], ['aa', $aa_pos] ];
		foreach my $pairs (@$arrayRef){
			my ($type, $pos) = @$pairs;
#								print "$type $pos\n";
			if (exists($nuc_aa_codon_tally->{$type}->{$pos})){
				print_top_hits($nuc_aa_codon_tally->{$type}->{$pos});
			}
			else {
				print $merged_report_fh "\t"x8; 	# print blank fields.  Shouldn't have any because the nucleotide sequences are trimmed to the first codon before tallying up nucleotides and amino acids.  Although maybe not trimmed at the 3' end of reads...
			}
		}
		
		# Now print the counts, number of amino acids represented, ref aa counts, and nonref aa counts. 
		
		print $merged_report_fh "\n"; 
	}
}
#-----------------------------------------------------------------------------
sub get_position_in_codon {
	my $nuc_pos = shift;
	my $position_in_codon = 0;
	my $remainder = $nuc_pos % 3; 
	# If remainder is 0, then position is 3; if remainder is 2, position is 2; if remainder is 1, position is 1.  
	$position_in_codon = $remainder || 3;
	return $position_in_codon;
}
#-----------------------------------------------------------------------------
sub print_top_hits {
	# Prints to $merged_report_fh
	# E.g., topNucleotide	secondNucleotide	thirdNucleotide	otherNucleotide	countTopNucleotide	countSecondNucleotide	countThirdNucleotide	countOtherNucleotide	
	# Keep in mind that some will only have one or two nucleotides (or aa or codons) in the tally.  In that case, just keep the other columns blank. 
	my ($tally) = @_;
#				print Dumper($tally);
	my $the_rest = "";
	my $the_rest_count = 0;
	my @top_three;
	my @top_three_counts;
	my $i = 0;
	foreach my $key (sort {$tally->{$b} <=> $tally->{$a}} keys %$tally){
		if ($i<3){
			push @top_three, $key;
			push @top_three_counts, $tally->{$key};
		}
		else {
			$the_rest .= $key.",";
			$the_rest_count += $tally->{$key};
		}		
		$i++;
	}
	$the_rest =~ s/,$//;
	push @top_three, $the_rest;
	push @top_three_counts, $the_rest_count;
#	print Dumper(\@top_three, \@top_three_counts);
	foreach my $array (\@top_three, \@top_three_counts){
		for (my $i = 0; $i < 4; $i++){		# Added "the_rest" to these arrays, so now there are 4 elements
			if ($array->[$i]){
				print $merged_report_fh "\t$array->[$i]";
			}
			else {
				print $merged_report_fh "\t";
			}
		}
	}
}
#-----------------------------------------------------------------------------
sub print_top_hits_alternative {
	# Prints to $merged_report_fh
	# E.g., topNucleotide	secondNucleotide	thirdNucleotide	otherNucleotide	countTopNucleotide	countSecondNucleotide	countThirdNucleotide	countOtherNucleotide	
	# Keep in mind that some will only have one or two nucleotides (or aa or codons) in the tally.  In that case, just keep the other columns blank. 
	my ($tally) = @_;
#				print Dumper($tally);
	my $the_rest = "";
	my $the_rest_count = 0;
	my @top_three;
	my @top_three_counts;
	my $i = 0;
	foreach my $key (sort {$tally->{$b} <=> $tally->{$a}} keys %$tally){
		if ($i<3){
			push @top_three, $key;
			push @top_three_counts, $tally->{$key};
		}
		else {
			$the_rest .= $key.",";
			$the_rest_count += $tally->{$key};
		}		
		$i++;
	}
	$the_rest =~ s/,$//;
	push @top_three, $the_rest;
	push @top_three_counts, $the_rest_count;
#					print Dumper(\@top_three, \@top_three_counts);
	for (my $i = 0; $i < 4; $i++){		# Added "the_rest" to these arrays, so now there are 4 elements
		foreach my $array (\@top_three, \@top_three_counts){
			if ($array->[$i]){
				print $merged_report_fh "\t$array->[$i]";
			}
			else {
				print $merged_report_fh "\t";
			}
		}
	}
}
#-----------------------------------------------------------------------------
sub print_variants {
	my ($nuc_aa_codon_tally, $label) = @_;
	# Filters the variants (nuc, aa, codon) for those that are above a certain threshold in frequency and then performs linkage disequilibrium in a pairwise manner for all variants at level of nuc, codon, or aa.
	
	# First get the variants that are above a frequency threshold
		#Type\tgene\tposition\tconsensus\tvariant\tcount\tcoverage\tfrequency\n";
	 	# type is nucleotide, codon or aa.  
	 	# position is either aa or nuc position
	 	# consensus is the one that is highest
	 	# variant is the one with frequency above threshold
	 	# Count is the count for this particular variant
	 	# Coverage is the number of reads at this position
	 	# Frequency is count/coverage
	print STDERR "Printing out variants report. ";
	&elapsed($start_time, 'Elapsed', $verbose);

	my $variants_pass;	#hashref of same structure as $nuc_aa_codon_tally, only has variants to consider for linkage disequilibrium.  
	my $max_coverage = 0; 	
	foreach my $type (qw(nuc codon aa)){	# was foreach my $type (keys %$nuc_aa_codon_tally){			#e.g., $nuc_aa_codon_tally->{'aa'}->{$pos}->{$aa}++;
		my $i = 0;
		foreach my $pos (sort {$a <=> $b} keys %{$nuc_aa_codon_tally->{$type}}){
			my $consensus = find_key_with_biggest_value($nuc_aa_codon_tally->{$type}->{$pos}) || "";
			my $coverage = total(values(%{ $nuc_aa_codon_tally->{$type}->{$pos} }));
#								print "type: $type\tposition: $pos\tcoverage: $coverage\n";
			$max_coverage = $coverage if ($coverage > $max_coverage);
#					next unless ($coverage >= $max_coverage * 0.05);		# Ignore low-coverage regions.  5% of max is just arbitrary.  This could be a problem if the reads start out with low coverage.  With amplicon sequencing, that won't be a problem.  
			my $filter = "";
			$filter .= "LOW_COVERAGE," 	if ($coverage < $max_coverage * 0.05);
			$filter .= "CONSENSUS_N," 	if ($consensus =~ m/N/ && $type ne 'aa');		# N is asparagine amino acid, which is valid.  
			$filter .= "CONSENSUS_X," 	if ($consensus =~ m/X/);
			$filter .= "CONSENSUS_DEL,"	if ($consensus =~ m/-/);
			my $threshold_count = $variant_threshold * $coverage; 
			foreach my $variant (keys (%{ $nuc_aa_codon_tally->{$type}->{$pos} })){
				next if ($variant eq $consensus);
								if ($verbose){	print STDERR "pos: $pos consensus: $consensus variant: $variant\n";	}
				my $variant_filter = $filter;
				$variant_filter .= "VARIANT_N," 	if ($variant =~ m/N/ && $type ne 'aa');		# N is asparagine amino acid, which is valid. Only assign this message if you are looking at codons or individual nucleotides.
				$variant_filter .= "VARIANT_X," 	if ($variant =~ m/X/);
				$variant_filter .= "VARIANT_DEL,"	if ($variant =~ m/-/);
				$variant_filter = "PASS" unless ($variant_filter);		# After all filters checked.
				my $count = $nuc_aa_codon_tally->{$type}->{$pos}->{$variant};
				if ($count >= $threshold_count){
					my $frequency = sprintf("%.5f",  $nuc_aa_codon_tally->{$type}->{$pos}->{$variant} / $coverage);
					print $variants_fh join "\t", $type, $gene, $pos, $consensus, $variant, $count, $coverage, $frequency, $variant_filter;
					print $variants_fh "\n";
#							if ($variant_filter eq "PASS"){
#								$variants_pass->{$type}->[$i]->{'pos'} = $pos; 
#								$variants_pass->{$type}->[$i]->{'variant'} = $variant;
#								$variants_pass->{$type}->[$i]->{'filter'} = $variant_filter;  # Save the filter so I can check it later.  Or I could just save only those that have "PASS" filter...
#								$i++;	# increment the variant count.
#							}
				}
			}
		}
	}
}
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
