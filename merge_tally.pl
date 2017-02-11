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
use primerid;

use Cwd;
use diagnostics; 
use Getopt::Long;
use Data::Dumper;
#use File::Basename;
#use Array::Utils qw(:all);
#use File::Slurp qw( prepend_file );
#use File::Which;


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
#my $mafft_bin = $PWD.'/mafft'; # philip
my $mafft_bin = 'mafft';

my $save;
my $dir;
my $verbose;
my $output;
my $gzip;
my $ref;
my $debug;
my $fullpeptide;
my $prefix;
my $variant_threshold;
my $use_sample_col;
my $sample;
GetOptions('save=s' => \$save, 
	'output=s' => \$output, 
	'verbose' => \$verbose, 
	'dir|d|i|input|files=s' => \$dir,
	'gzip' => \$gzip, 
	'ref|r=s' => \$ref, 
	'debug' => \$debug, 
	'fullpeptide=s' => \$fullpeptide, 
	'prefix=s' => \$prefix, 
	'variant_threshold=s' => \$variant_threshold, 
	'sample=s' => \$sample,
	'use_sample_col' => \$use_sample_col,
);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "
merge_tally.pl takes a directory containing output from convert_reads_to_amino_acid.pl  
and merges the results.  The merged results are  output as frequency tables of nucleotide, 
codon, and amino acid sequence by position -- the same format as the output to 
convert_reads_to_amino_acid.pl.  A merged table containing frequencies for all three is also 
produced.  For now, Phylip files for nucleotide and peptide sequence are not merged.

Required Arguments:
--d/--dir/-i/--input 	Directory containing files to merge.
--sample	Sample Id for all of the reads.  This will force merge all of the tally files 
		in the input directory, regardless of the 'Sample' column.  No spaces or commas 
		allowed in the id.
OR
--use_sample_col	Use the 'Sample' column in the output of convert_reads_to_amino_acid.pl 
		to group reads into separate sets of output files, one set of files for each unique 
		sample id.  

Optional Arguments:
--prefix	Prefix for output files.  Default = Merge_convert_reads.  Sample id will be 
		added to the prefix.  If the prefix value contains a directory, the directory must exist.
--variant_threshold	Minimum threshold frequency for output variant file.  Default = 0 
		(i.e., print all variants).  (Note that in calculate_linkage_disequilibrium.pl
		there is an option to filter for a variant_threshold as well.)
--verbose	Print out verbose progress/error messages.

Example:
merge_tally.pl --dir 30_S1_tallies --sample Parental --prefix 30_S1.contigs.pi.btrim.0.majority.cons 

";



# To Do:


# Change log
# 2017-02-04
# Created script, based on convert_reads_to_amino_acid.pl
# 2017-02-11
# Removed --save option.  Just include any directory in the --prefix.  ($prefix can contain absolute or relative path).
# --save		Directory in which to save files. Default = pwd.  If folder doesn\'t exist, it 
#			will be created.


unless ($dir||$ARGV[0]){	print STDERR "$usage\n";	exit;	}	#fasta and gff are required
unless ($sample || $use_sample_col){	print STDERR "Please specify --sample or --use_sample_col.\n";	exit;	}
if ($use_sample_col && $sample){
	print STDERR "Provide a sample id (--sample) or use the Sample column in the tally files (--use_sample_col), not both.\n"; 
	exit 1;
}
if ($sample && $sample =~ m/[\s+,\|]/){
	print STDERR "Sample id invalid: $sample\n";
	exit 1;
}


my $start_time = time;

my @suffixes = (qw(.bed .bed.gz .bed12 .bed12.gz .txt .txt.gz .BED .BED.gz .BED12 .BED12.gz .fasta .fa .FA .FASTA .FAS .fas), @SUFFIXES);	#for fileparse.  Feel free to add more accepted extensions.  @SUFFIXES comes from aomisc.pm.  



# Defaults for input options
$prefix ||= 'Merge_convert_reads';
$variant_threshold //= 0;		# Was 0.0001 (1 variant out of 10000 sequences).  0 is a valid value, hence // instead of ||.


# Make hash to convert codons to amino acids (maybe this should go into primerid.pm too?)

my %converter = make_codon_to_aa_hash();

# Get the tally files
my @all_files = &aomisc::get_files($dir, '(nuc|aa|codon).tally.xls'); 	# This will add $dir to each of the file names.  Was: #opendir(DIR, "$dir"); my @all_files  = grep(/\.(nuc|aa|codon).tally.xls$/,readdir(DIR));  closedir(DIR);

my @nuc_files	= grep(/\.nuc.tally.xls$/,@all_files);
my @codon_files	= grep(/\.codon.tally.xls$/,@all_files);
my @aa_files 	= grep(/\.aa.tally.xls$/,@all_files);



	# 	print STDERR join "\n", @all_files, "\n";

# Get the gene name.  Note that all of the file should be for the same gene.  Check to make sure.  Die if not.
my $gene = get_gene(\@all_files);

my @NUC = qw(A C G T N);
my @AA = qw(A C D E F G H I K L M N P Q R S T V W Y * X -); 	# All amino acids, and * for stop, X for unknown, and - for deletion in read compared to reference
my @CODON = get_codons("$codon_files[0]");	# Take the first file named *.codon.tally.xls in the directory and read the codon order from there.

	#	print join " ", @CODON; print "\n"; exit;


# Read all files and store in data structure.

# Reference CDS sequence -- coding nucleotides, codons, and amino acids by position.
	# Global variable to store the reference sequence
	# $cds->{'aa'}->{$aa_pos} 	= $aa;
	# Other first keys are codon, nuc
	# Store this as we're reading the input files

my ($sample_nuc_aa_codon_tally, $cds) = read_input_tally(\@nuc_files, \@codon_files, \@aa_files);
	#			print Dumper($nuc_aa_codon_tally); exit;
print STDERR "Done processing " . scalar(@all_files) . " files.";
	&elapsed($start_time, ' Elapsed', $verbose);

my @samples = keys %$sample_nuc_aa_codon_tally;

	# 	print Dumper(\@samples); exit;

for (my $i = 0; $i < @samples; $i++){
	my $prefix_with_sample = $prefix . "." . $samples[$i];
	print_reports($sample_nuc_aa_codon_tally->{$samples[$i]}, $samples[$i], $gene, \@NUC, \@CODON, \@AA, $cds, \%converter, $prefix_with_sample, $verbose, $debug);		
	print_merged_report($sample_nuc_aa_codon_tally->{$samples[$i]}, $samples[$i], $gene, \@CODON, $cds, $prefix_with_sample);	
	print_variants($sample_nuc_aa_codon_tally->{$samples[$i]}, $samples[$i], $gene, $variant_threshold, $prefix_with_sample, $start_time, $verbose);
}



&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub get_gene {
	my $all_file_list = shift;
	# Reads all of the files, getting gene from column 2
	# If there is only one gene across files, then it returns the gene name
	# If there is more than one gene, it will exit.

	# #name	gene	nucleotidePosition	refNucleotide	consensusNucleotide	refDiffConsensus	Sample	coverageDepth	unambigCoverageDepth	unambigConsensus	numUnambigConsensus	numUnambigNonConsensus	majorAltAllele	numMajorAltAllele	numA	numC	numG	numT	numN	numOther
	# HA_full:37:A	HA_full	37	A	A		151_78_S2.0.000	9922	9922	A	9889	33	G	23	9889	3	23	7	0	0
	# HA_full:38:G	HA_full	38	G	G		151_78_S2.0.000	9922	9922	G	9918	4	A	3	3	0	9918	1	0	0

	my $genes;	# hashref to store all gene id(s).  Should only have one gene across all files, theoretically.
	print STDERR "Checking to see if same gene listed for all files...\n";

	foreach my $file (@$all_file_list){
		my $readfh = open_to_read("$file");
		while(<$readfh>){
			chomp;
			next if (m/^#/);

			my @F = split(/\t/);
			next unless (@F > 1);

			$genes->{$F[1]}++;	# gene is second column.
		}
		close($readfh);
	}
	
	#		print Dumper($genes);
	my @genes = keys %$genes;
	if (@genes > 1){
		print STDERR "More than one gene present in files!  Exiting...\n";
		exit;
	}
	
	print STDERR "Gene: $genes[0]\n";

	return $genes[0];
}
#-----------------------------------------------------------------------------
sub get_codons {
	my $codon_file = shift;

	#	print STDERR "$codon_file\n";

	my @codons;

	my @header = get_header($codon_file);
	# #name	gene	codonPosition	refCodon	consensusCodon	refDiffConsensus	Sample	coverageDepth	unambigCoverageDepth	unambigConsensus	numUnambigConsensus	numUnambigNonConsensus	majorAltAllele	numMajorAltAllele	numGCA:A	numGCB:A	numGCC:A	numGCD:A	numGCG:A	numGCH:A	numGCK:A	numGCM:A	numGCN:A	numGCR:A	numGCS:A	numGCT:A	numGCV:A	numGCW:A	numGCY:A	numTGC:C	numTGT:C	numTGY:C	numGAC:D	numGAT:D	numGAY:D	numGAA:E	numGAG:E	numGAR:E	numTTC:F	numTTT:F	numTTY:F	numGGA:G	numGGB:G	numGGC:G	numGGD:G	numGGG:G	numGGH:G	numGGK:G	numGGM:G	numGGN:G	numGGR:G	numGGS:G	numGGT:G	numGGV:G	numGGW:G	numGGY:G	numCAC:H	numCAT:H	numCAY:H	numATA:I	numATC:I	numATH:I	numATM:I	numATT:I	numATW:I	numATY:I	numAAA:K	numAAG:K	numAAR:K	numCTA:L	numCTB:L	numCTC:L	numCTD:L	numCTG:L	numCTH:L	numCTK:L	numCTM:L	numCTN:L	numCTR:L	numCTS:L	numCTT:L	numCTV:L	numCTW:L	numCTY:L	numTTA:L	numTTG:L	numTTR:L	numATG:M	numAAC:N	numAAT:N	numAAY:N	numCCA:P	numCCB:P	numCCC:P	numCCD:P	numCCG:P	numCCH:P	numCCK:P	numCCM:P	numCCN:P	numCCR:P	numCCS:P	numCCT:P	numCCV:P	numCCW:P	numCCY:P	numCAA:Q	numCAG:Q	numCAR:Q	numAGA:R	numAGG:R	numAGR:R	numCGA:R	numCGB:R	numCGC:R	numCGD:R	numCGG:R	numCGH:R	numCGK:R	numCGM:R	numCGN:R	numCGR:R	numCGS:R	numCGT:R	numCGV:R	numCGW:R	numCGY:R	numAGC:S	numAGT:S	numAGY:S	numTCA:S	numTCB:S	numTCC:S	numTCD:S	numTCG:S	numTCH:S	numTCK:S	numTCM:S	numTCN:S	numTCR:S	numTCS:S	numTCT:S	numTCV:S	numTCW:S	numTCY:S	numACA:T	numACB:T	numACC:T	numACD:T	numACG:T	numACH:T	numACK:T	numACM:T	numACN:T	numACR:T	numACS:T	numACT:T	numACV:T	numACW:T	numACY:T	numGTA:V	numGTB:V	numGTC:V	numGTD:V	numGTG:V	numGTH:V	numGTK:V	numGTM:V	numGTN:V	numGTR:V	numGTS:V	numGTT:V	numGTV:V	numGTW:V	numGTY:V	numTGG:W	numTAC:Y	numTAT:Y	numTAY:Y	numTAA:*	numTAG:*	numTGA:*	numOther
	@codons = @header[14..$#header -1];	# Not "numOther".  This is added by print_reports (via print_tally_line())

	s/num// for @codons;
		print STDERR "codons: ", join ", ", @codons, "\n" if $verbose;

	return @codons;		# Return the array (not arrayref)
}
#-----------------------------------------------------------------------------
sub read_input_tally {
	my ($nuc_files, $codon_files, $aa_files) = @_;  # lists of files inside $dir (need to add "$dir/" to the file name to access the file)
	my $sample_nuc_aa_codon_tally;  # data structure for saving all of the count information. e.g., $sample_nuc_aa_codon_tally->{$sample}->{$type}->{$pos}->{$base_codon_aa} += $seq_count;
	my $cds;	# data structure for saving reference sequence as nuc, codon, aa.  e.g., $cds->{'nuc'}->{$pos}		= $nuc;  first keys: nuc, aa, codon
		# Global variable to store the reference sequence
		# $cds->{'aa'}->{$aa_pos} 	= $aa;
		# Other first keys are codon, nuc
		# Store this as we're reading the input files


	## Example input.
	# nuc:
	# #name	gene	nucleotidePosition	refNucleotide	consensusNucleotide	refDiffConsensus	Sample	coverageDepth	unambigCoverageDepth	unambigConsensus	numUnambigConsensus	numUnambigNonConsensus	majorAltAllele	numMajorAltAllele	numA	numC	numG	numT	numN	numOther
	# HA_full:37:A	HA_full	37	A	A		151_78_S2.0.000	9922	9922	A	9889	33	G	23	9889	3	23	7	0	0
	# HA_full:38:G	HA_full	38	G	G		151_78_S2.0.000	9922	9922	G	9918	4	A	3	3	0	9918	1	0	0
	
	# codon:
	# #name	gene	codonPosition	refCodon	consensusCodon	refDiffConsensus	Sample	coverageDepth	unambigCoverageDepth	unambigConsensus	numUnambigConsensus	numUnambigNonConsensus	majorAltAllele	numMajorAltAllele	numGCA:A	numGCB:A	numGCC:A	numGCD:A	numGCG:A	numGCH:A	numGCK:A	numGCM:A	numGCN:A	numGCR:A	numGCS:A	numGCT:A	numGCV:A	numGCW:A	numGCY:A	numTGC:C	numTGT:C	numTGY:C	numGAC:D	numGAT:D	numGAY:D	numGAA:E	numGAG:E	numGAR:E	numTTC:F	numTTT:F	numTTY:F	numGGA:G	numGGB:G	numGGC:G	numGGD:G	numGGG:G	numGGH:G	numGGK:G	numGGM:G	numGGN:G	numGGR:G	numGGS:G	numGGT:G	numGGV:G	numGGW:G	numGGY:G	numCAC:H	numCAT:H	numCAY:H	numATA:I	numATC:I	numATH:I	numATM:I	numATT:I	numATW:I	numATY:I	numAAA:K	numAAG:K	numAAR:K	numCTA:L	numCTB:L	numCTC:L	numCTD:L	numCTG:L	numCTH:L	numCTK:L	numCTM:L	numCTN:L	numCTR:L	numCTS:L	numCTT:L	numCTV:L	numCTW:L	numCTY:L	numTTA:L	numTTG:L	numTTR:L	numATG:M	numAAC:N	numAAT:N	numAAY:N	numCCA:P	numCCB:P	numCCC:P	numCCD:P	numCCG:P	numCCH:P	numCCK:P	numCCM:P	numCCN:P	numCCR:P	numCCS:P	numCCT:P	numCCV:P	numCCW:P	numCCY:P	numCAA:Q	numCAG:Q	numCAR:Q	numAGA:R	numAGG:R	numAGR:R	numCGA:R	numCGB:R	numCGC:R	numCGD:R	numCGG:R	numCGH:R	numCGK:R	numCGM:R	numCGN:R	numCGR:R	numCGS:R	numCGT:R	numCGV:R	numCGW:R	numCGY:R	numAGC:S	numAGT:S	numAGY:S	numTCA:S	numTCB:S	numTCC:S	numTCD:S	numTCG:S	numTCH:S	numTCK:S	numTCM:S	numTCN:S	numTCR:S	numTCS:S	numTCT:S	numTCV:S	numTCW:S	numTCY:S	numACA:T	numACB:T	numACC:T	numACD:T	numACG:T	numACH:T	numACK:T	numACM:T	numACN:T	numACR:T	numACS:T	numACT:T	numACV:T	numACW:T	numACY:T	numGTA:V	numGTB:V	numGTC:V	numGTD:V	numGTG:V	numGTH:V	numGTK:V	numGTM:V	numGTN:V	numGTR:V	numGTS:V	numGTT:V	numGTV:V	numGTW:V	numGTY:V	numTGG:W	numTAC:Y	numTAT:Y	numTAY:Y	numTAA:*	numTAG:*	numTGA:*	numOther																							
	# HA_full:13:AGG	HA_full	13	AGG	AGG		151_78_S2.0.000	9922	9922	AGG	9877	45	GGG	23	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	23	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	7	9877	0	0	0	0	0	3	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	7	0	0	0	0	0	0	0																							
	# HA_full:14:CAA	HA_full	14	CAA	CAA		151_78_S2.0.000	9922	9922	CAA	9888	34	CAG	15	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	9888	15	0	0	0	0	9	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	7	0	0	0																							

	# aa:
	# #name	gene	aminoAcidPosition	refAminoAcid	consensusAminoAcid	refDiffConsensus	Sample	coverageDepth	unambigCoverageDepth	unambigConsensus	numUnambigConsensus	numUnambigNonConsensus	majorAltAllele	numMajorAltAllele	numA	numC	numD	numE	numF	numG	numH	numI	numK	numL	numM	numN	numP	numQ	numR	numS	numT	numV	numW	numY	num*	numX	num-	numOther
	# HA_full:13:R	HA_full	13	R	R		151_78_S2.0.000	9922	9922	R	9887	35	G	23	0	0	0	0	0	23	0	0	3	0	1	0	0	0	9887	1	0	0	7	0	0	0	0	0
	# HA_full:14:Q	HA_full	14	Q	Q		151_78_S2.0.000	9922	9922	Q	9903	19	R	9	0	0	0	1	0	0	0	0	1	1	0	0	0	9903	9	0	0	0	0	0	7	0	0	0

	my @type = qw(nuc codon aa);
	my @file_lists_by_type = ($nuc_files, $codon_files, $aa_files);

	# Tally up all "num..." columns except "numOther" -- "numOther" is added by print_reports (actually print_tally_line)"
	# Need to remove the "num" part

	# I can either assume that the order of @NUC, @AA, and @CODON is the same order as the columns (which it should be), or I can lookup the header for each file and use that.
	# I would hate to be off by one for some reason, so I'll just lookup the header for each file.

	my $files_done = 0;

	for (my $i = 0; $i < @type; $i++){
		my $type = $type[$i];
		my $file_list = $file_lists_by_type[$i];

		my @header = get_header("$file_list->[0]"); 	# get header array from first file

		s/^num// for @header;
		s/:.+// for @header;

		foreach my $file (@$file_list){
				print STDERR "reading file: $file\n" if $verbose;
			my $wc = `cat $file | wc -l`;
			# 	print STDERR "wc: $wc\n";
			if ($files_done && ($files_done % 100) == 0){
				my $total_files = scalar(@all_files);
				print "Done processing $files_done of $total_files files.";
				&elapsed($start_time, ' Elapsed', $verbose);
			}
			$files_done++;
			next if $wc < 2;

			my $readfh = open_to_read("$file");
			while(<$readfh>){
				chomp;
				next if (m/^#/);

				my @F = split(/\t/);
				#	my ($pos,$ref_residue,$sample) = ($F[2],$F[3],$F[6]);
				my ($pos,$ref_residue) = ($F[2],$F[3],"SV12_parental");
				my $sample = $sample ? $sample : $F[6];		# If user specifies $sample, use that as the sample id for all reads.  Otherwise, use the 'Sample' column value.  

				# Populate $cds for reference
				if (exists($cds->{$type}->{$pos}) && $cds->{$type}->{$pos} ne $ref_residue){
					print STDERR "Reference residue $ref_residue different from stored value for this type ($type) and position ($pos): $cds->{$type}->{$pos}\n";
				} else {
					$cds->{$type}->{$pos} 	= $ref_residue;
				}

				for (my $c = 14; $c < (scalar(@header) - 1); $c++){		 # i.e., start in column index 14, go all the way to the second to last column.  The last column is numOther, which we don't want to capture
					$sample_nuc_aa_codon_tally->{$sample}->{$type}->{$pos}->{$header[$c]} += $F[$c] if ($F[$c] > 0);
				}
			}
	#				print Dumper($sample_nuc_aa_codon_tally); exit;

			close($readfh);

		}

	}


	return ($sample_nuc_aa_codon_tally, $cds);

}
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
