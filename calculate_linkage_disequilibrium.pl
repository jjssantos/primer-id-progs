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
use File::Temp;
use Statistics::R;
use Parallel::Loops;


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



#Print out the options
if (@ARGV){		print STDERR "Arguments: ", join " ", @ARGV, "\n";	}

my $save;
my $files;
my $verbose;
my $output;
my $gzip;
my $prefix;
my $nofilter;
my $variant_threshold = 0.001;
my $cpu = 2;
my $force;
my $label = "Sample";
my $group_id = "Group";
my $table_only;
GetOptions('save=s' => \$save, 
	'output=s' => \$output, 
	'verbose' => \$verbose, 
	'files=s' => \$files, 
	'gzip' => \$gzip, 
	'prefix=s' => \$prefix, 
	'nofilter' => \$nofilter, 
	'variant_threshold=s' => \$variant_threshold, 
	'cpu|p=s' => \$cpu, 
	'force' => \$force, 
	'label|labels=s' => \$label, 
	'group_id=s' => \$group_id, 
	'table_only' => \$table_only,
	);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "
calculate_linkage_disequilibrium.pl Takes a file of variants and sequences and checks for 
linkage between each pair of variants.  Variants are queried at the level of nucleotide, 
codon, and amino acid.  Linkage is computed with fisher's exact test and is reported as 
p-value, FDR-adjusted p-value, and Odds ratio (OR).  This script requires R on the PATH.  

Usage: calculate_linkage_disequilibrium.pl <variants.xls> <cleanreads.txt> [<cleanpeptides.txt>]

variants.xls, cleanreads.txt, and cleanpeptides.txt are output files from 
convert_reads_to_amino_acid.pl.  
variants.xls and cleanreads.txt are required; cleanpeptides.txt is recommended as well.
e.g., 
Convert_reads.variants.xls 
## Variants above threshold frequency 0.005
#Type	gene	position	consensus	variant	count	coverage	frequency	filter
codon	AAM75158.1	159	AGN	AGT	424	19003	0.0223	CONSENSUS_N,
codon	AAM75158.1	164	TTG	NTG	1838	19003	0.0967	VARIANT_N,
codon	AAM75158.1	166	TGG	TNG	135	19003	0.0071	VARIANT_N,
codon	AAM75158.1	167	CTG	CTN	108	19003	0.0057	VARIANT_N,
codon	AAM75158.1	176	GAG	AAG	127	19003	0.0067	PASS

head *cleanreads.txt
  24002 431
#readID:gene:positionForFirstNucleotideBase	cleanRead
Y3YRC_ATAAGCTT_3|AAM75158.1|158	AGCAGNTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGNCCTTGTACTGTGGGGTATTCATCACCCGCCTAACAGTAA
Y3YRC_TAGAGCTT_6|AAM75158.1|158	AGCAGNTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGNCCTTGTACTGTGGGGTATTCATCACCCGCCTAACAGTAA
Y3YRC_GGTAATCG_2|AAM75158.1|158	AGCAGNTTTTACAGAAATNTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTNAACAAAAAAAGGAAAGAAGNCCTTGTACTGTGGGGTATTCATCACCCGCCTAACAGTAA

head *cleanpeptides.txt
  24002 143
#readID:gene:positionForFirstAminoAcid	cleanPeptide
Y3YRC_ATAAGCTT_3|AAM75158.1|472	SXFYRNLLWLTKKESSYPELKNSYVNKKRKEXLVLWGIHHPPNS
Y3YRC_TAGAGCTT_6|AAM75158.1|472	SXFYRNLLWLTKKESSYPELKNSYVNKKRKEXLVLWGIHHPPNS
Y3YRC_GGTAATCG_2|AAM75158.1|472	SXFYRNXLWLTKKESSYPELKNSYVNKKRKEXLVLWGIHHPPNS

Example output:
cat Calculate_linkage_disequilibrium_report.xls 
#group	sample	type	gene	comparison	pos1	c1	v1	v1freq	pos2	c2	v2	v2freq	c1c2	c1v2	v1c2	v1v2	p-value	OR	FDR
Group	Sample	nuc	AAM75158.1	526:G:A:556:A:G	526	G	A	0.006468859	556	A	G	0.015343491	18329	255	89	32	8.62E-31	25.82045	8.62E-31
Group	Sample	codon	AAM75158.1	176:GAG:AAG:186:AGG:GGG	176	GAG	AAG	0.006420546	186	AGG	GGG	0.015355805	18315	255	88	32	6.59E-31	26.09373	6.59E-31
Group	Sample	aa	AAM75158.1	176:E:K:186:R:G	176	E	K	0.006419859	186	R	G	0.015354162	18317	255	88	32	6.57E-31	26.09658	6.57E-31

group is the Group name, specified by --group_id
sample is the label for the sample, specified by --label
type is the type of variant being considered, either nuc, codon, or aa.
gene is the gene name from the variants.xls input file that the variant corresponds to
comparison is a colon-delimited combination of pos1:c1:v1:pos2:c2:v2 that can be used as a unique identifier
c1 and v1 refer to consensus and variant, respectively, in the first position
c2 and v2 refer to consensus and variant, respectively, in the second position
v1freq and v2freq are the variant frequencies from the variants.xls input file
pos1 and pos2 refer to positions (either amino acid/codon or nucleotide) of the first and second variants, respectively
c1c2 is the count of reads that have consensus at both positions
c1v2 is the count of reads that have consensus at the first position and variant at the second position
v1c1 is the count of reads that have variant at the first position and consensus at the second position
v1v2 is the count of reads that have variant at both positions
OR is the odds ratio.  In this case, it is the likelihood that the sequence at the first position 
	affects the sequence at the second position.  
	OR > 1 indicates that variants for both positions or consensus sequences for both positions are often found together
	OR < 1 indicates that having a variant at one position and a consensus at the other is the common trend.
	OR = 1 indicates that the two positions are not linked.
	* Note that the OR should not be considered meaningful unless the FDR value shows that the relationship is significant (e.g., < 0.05).  
p-value is the result from the fisher's exact test
FDR is the adjusted p-value from the Benjamini & Hochberg method

Replicates
If multiple replicates for a sample are available, it is recommended that these be 
processed together.  The data for different replicates is not merged, but information for
all replicates is used to include all variants that pass the threshold in at least one
replicate.  (The data can be merged by feeding the output of this script to 
combine_linkage_values.pl .)  To process replicates together, provide a comma-delimited 
list of each argument type.  If you are inputting files for replicates, the --group_id and 
--label parameters are recommended.  For example:
calculate_linkage_disequilibrium.pl --group_id Immunized --label Sample1,Sample2,Sample3 variants1.xls,variants2.xls,variants3.xls Sample1.cleanreads.txt,Sample2.cleanreads.txt,Sample3.cleanreads.txt Sample1.cleanpeptides.txt,Sample2.cleanpeptides.txt,Sample3.cleanpeptides.txt

OPTIONS:
	--save		Directory in which to save files. Default = pwd.  If folder doesn't exist, 
			it will be created.
	--prefix	Prefix for output files. '<--group_id value>.<label[i]>' will be added to
			the prefix.
	--variant_threshold	Minimum threshold to consider for a variant.  Default is 0.001.     
	-p/--cpu	Number of CPUs to use for parallel processing.  Default = 2.  Speed-up is 
			~0.5x per cpu.  E.g., for 8 cpu, I get ~4x speedup compared to using 1 cpu. 
			8-12 cpus is optimal (or number of physical cores) and time estimate is most
			accurate in that range as well.  
	--label		Single label for sample or a comma-delimited list of labels if 
			replicates are used. Default is 'Sample'.  Default when replicates are detected
			is a list of numbers such as '1,2,3', depending on the number of replicates
			detected.
	--group_id	Group id.  A name to give to this sample (or samples) to use for the Group 
			column.  Useful for naming output files that were processed together, especially
			if you are going to combine multiple files afterwards.  Default is 'Group'.  
	--table_only	Don't calculate p-values, just print out the table of counts.
	--force		Force the analysis to run, even if expected runtime is > 24h.  Default is
			to only run if low range of expected is <= 24h. 
	--nofilter	Compare all variants in variants.xls file.  Default is to only consider 
			variants that have PASS status.  
";

# To do
# Maybe add another option to specify which filters to allow (comma-delimited)
# Add threads to 
# Compute quasi-cliques with DENSE.  (or Cocain?).  Add description of requirements (DENSE/mgqce) to usage.
# Add fdr threshold option for quasi-clique variants. Default = 0.05
# Add part to count the number of pairs with fdr <= 0.05 while printing.
# Redo time estimate after changes on 2017-02-26.

# Change log
# 2013-06-05 
# Added $variant_threshold option.  This way I can run this script multiple times on the same variant file with different frequency thresholds, rather than running convert_reads_to_amino_acid every time I want a different threshold.  
# 2013-09-23
# Changed get_reads subroutine to use '|' as delimiter for cleanpeptides.txt and cleanreads.txt input files to separate readID, gene, position.  
# Changed get_reads sub to skip PHYLIP header lines:  next if (m/^\s+\d+\s+\d+/);	
# 2013-09-27
# Added option --cpu and made a requirement for Parallel::Loops module.  This now allows parallel processing for calculating p-values in R.  Speedup is # cpu/2 (e.g., for 8 cpus, I get about 4-fold speedup).    
# Added --force
# 2014-01-07
# Added -p to input to --cpu.
# Modified read_variants to adjust the count to remove comparisons of alternative alleles of the same position and also modified linkage_disequilibrium sub to make sure that if two variants have the same position that they won't be compared to one another. 
# Adjust time estimation to account for the number of reads in the file.  
	# Ran read file size from 1000 - 9000 (every 1000) and -p from 2 to 20 (every 2) on niaid-1-11 for datapoints.  Used variants file with only codons, using threshold 0.005, which gave 9 variants (36 comparisons). Computed time_per_comparison = real_time (from "time" command) * p / 36.  
	# Also ran comparison with the same variant file and read file sizes but -p 1 to 4 (every 1) on laptop.
	# Equation from niaid-1-11 test: 	Y =  0.0700371 + number_reads * 0.0001235 + p * 0.1297   (Very good fit)
	# Equation from laptop test: 		Y = -0.0520766 + number_reads * 0.0001522 + p * 0.2966441 (Very good fit)
	# Equation for merged data: 		Y =  0.3099219 + number_reads * 0.0001317 + p * 0.1115561  		
# 2014-01-10
# Changed file output name to use dots '.' instead of underscores.
# 2014-04-02
# Modified get_reads() sub to also get the read_count from the id (to accomodate files that have been uniqued with fastx_collapser) and then modified the linkage_disequilibrium sub to use that count in incrementing the genotypes, e.g., $AB += $read_count;
# 2014-06-12
# Modified script to store the frequency of the variants (called by convert_reads_to_amino_acid.pl) and include that in the output report as well.  
# 2014-09-06
# Modified to make it work with replicates.  Added options --label --group_id, --table_only.  
# Now arguments can either be regular or comma-delimited lists if replicates are provided.
# Added "group", "sample", and "comparison" columns to the output file
# Added get_common_variant_list subroutine to get a list of variants that pass the threshold in at least one of the replicates (if replicates are provided).
# ("Process samples in groups of replicates so that I consider all variants that pass the threshold in at least one replicate."  )
# Modified linkage_disequilibrium subroutine to take the $variants as well (for frequency?)
# 2014-11-15
# Fixed time calculation in (and put it in get_common_variant_list) to work for multiple-file input with set of variants that are above threshold in at least one of the files.  Adds up separate times based on the size of the cleanreads.txt files and then gives a total estimate.  
# Fixed labels to default to 1, 2, 3, etc. if no --labels option is given and there are multiple files
# Add labels for group and Sample ids
# Changed default prefix to '<--group_id value>.<label[i]>' for multiple file input.  
# Fix loop to write a separate output file for each input file if multiple file input.
# Updated example output in $usage
# Modify the parallelization to be more efficient by sending a pair of variants to compare to a thread, instead of setting up threads for only the second variant.  (And change the time estimation calculation to match)
# Changed the structure of the $variants hashref in get_common_variant_list to be all hashes HoHoHoH (?) instead of HoHoAoH.  It had an array in there before because I wanted to walk through each of the variants that were present and do the comparisons that way.  But now, I need to do comparisons for some variants that are NOT present in a sample's variants.xls file.  It will look up the c1v1, c1v2, c2v1, c2v2 counts from the cleanreads.txt and cleanpeptides.txt files anyway.  (That is in the case where it was above threshold in only one of the replicates and was not found in one of the other replicates.) Now I just need to look up the frequency from that data structure.  Maybe it's not necessary to keep it around....
# Changed get_common_variant_list to make an AoA of all comparisons to be performed.  This way, each comparison can be separate
# Updated equation.  I couldn't quite remember what I used to fit the graph previously, but I used Excel this time, threshold of 0.004, which had 20 variants (187 codons).  Fit to equation Time_per_comparison = A + B * read_count + C * cpu.  Allow A, B, and C change.  Make column of values for "predicted_time_per_comparison" with that formula.  Then another column for diff^2 (Time_per_comparison-predicted_time_per_comparison)^2.  Make sum of the diff^2 column and set that to minimize.  For the test, I did p = 2 to 20 (by 2) and read_counts from 1000 to 9000 (by 1000) like last time.  
	# Equation from niaid-1-11 test: Y = 0.54984528 + number_reads * 8.6017*10^-6 + p * 8.96667*10^-9
	# Plotted results and found that optimal speed-up is found with -p 8-12.  
# 2017-02-26
# Updated linkage_disequilibrium sub to compute p-values outside of the Parallel Loop in order to avoid errors with checking version of R in Statistics::R
# From To do: "Fix bug - not carrying over frequency from input variants.xls file."  I checked and it seems to be fine.
# 2017-02-28
# Change prefix to be group.id even for single-sample, unless prefix is provided by user.

unless ($ARGV[1]){	print STDERR "$usage\n";	exit;	}
unless ($ARGV[2]){	print STDERR "No peptide file found.  Will only look at nucleotide and codon variants.\n";	}

my $start_time = time;

my @variants_files = split(/,/, $ARGV[0]);
my @clean_reads_files = split(/,/, $ARGV[1]);
my @clean_peptides_files = split(/,/, $ARGV[2]) if $ARGV[2];
my @labels = split(/,/, $label);

### Check for replicates and fix the --label default if so (see $usage)
# Throw an error if number of variants files is different than number of cleanreads files or cleanpeptide files. Required arguments
my $num_var_files = scalar(@variants_files);
my $num_clean_reads_files = scalar(@clean_reads_files);
my $num_clean_peptide_files = scalar(@clean_peptides_files) if $ARGV[2];
my $warned = 0;
if ($num_clean_reads_files != $num_var_files){
	print STDERR "Number of cleanreads.txt files doesn't match the number of variants files!\n";
	$warned++;
}
elsif($ARGV[2] && ($num_clean_peptide_files != $num_var_files) ){
	print STDERR "Number of cleanpeptides.txt files doesn't match the number of variants files!\n";
	$warned++;
}
exit if ($warned);
# Be more lenient with labels because not required argument
if ($label =~ m/^Sample$/ && $num_var_files > 1){
	# User didn't supply labels
	# Set new defaults for @labels
	@labels = (1..$num_var_files); 
}

# Get list of variants passing threshold in at least one file.  
	# If some are found, an estimate of the total processing time will be computed.  
		# If that value is > 24 hours, the script will die unless --force is also specified
	# If no variants passing threshold are found, this sub will also kill the script.
my ($variants_to_consider, $comparisons) = get_common_variant_list(\@variants_files);		# Get a list/hash of variants that pass the threshold in at least one replicate.  Hashref, $variants->{$type}->{$gene}->{$position}->{$variant_combined_name}++; 
		#  print "Variants to consider: \n", Dumper($variants_to_consider), Dumper($comparisons); exit;


my $save_dir = $save || Cwd::cwd();
unless (-d $save_dir){	mkdir($save_dir) or warn "$!\n"	}
warn "save directory: $save_dir\n";


# Walk through each input variant file and 
for (my $i = 0; $i < @variants_files; $i++){

	print STDERR "Processing file $variants_files[$i]...\n\t";
	
	my $filebase;
	if ($prefix){
		$filebase = $prefix . "." . $group_id . "." . $labels[$i];				
	}
	else {
		$filebase = $group_id . "." . $labels[$i];
	}
	my $linkage_disequil	= $save_dir . "/" . $filebase . ".linkage.minfreq".$variant_threshold.".xls";			# was $prefix . "_linkage_disequilibrium_report_minfreq".$variant_threshold.".xls";
	$linkage_disequil =~ s/_linkage.linkage/_linkage/;
	my $linkage_disequil_fh		= open_to_write("$linkage_disequil");
	print $linkage_disequil_fh "#group\tsample\ttype\tgene\tcomparison\tpos1\tc1\tv1\tv1freq\tpos2\tc2\tv2\tv2freq\tc1c2\tc1v2\tv1c2\tv1v2\tp-value\tOR\tFDR\n";		

	my $variants = read_variants($variants_files[$i]);		# HoHoAoH; first keys, type, i.e., 'nuc', 'codon', 'aa'; second keys gene; array of variants; hash with keys 'frequency', 'position', 'variant', 'consensus'
					#		 print Dumper($variants);  # exit;

	# Now read the nucleotide file *cleanreads.txt into memory and get linkage for variants of type 'nuc' and 'codon'
	my $clean_reads = get_reads($clean_reads_files[$i]);		#	print Dumper($clean_reads); exit; 
	my $stats = linkage_disequilibrium($clean_reads,$variants,'nuc',$labels[$i],$linkage_disequil_fh);		# print Dumper($stats); exit; 
	$stats = linkage_disequilibrium($clean_reads,$variants,'codon',$labels[$i],$linkage_disequil_fh);	# print Dumper($stats); exit; 	# I'll probably need a special step for codon

	# Flush the variable to free up memory
	undef($clean_reads);

	# Now, if a cleanpeptide.txt file is available, read it into memory 
	exit unless $ARGV[2];
	$clean_reads = get_reads($clean_peptides_files[$i]);
	$stats = linkage_disequilibrium($clean_reads,$variants,'aa',$labels[$i],$linkage_disequil_fh);		# Calculates pairwise linkage between variants.  

	close($linkage_disequil_fh);

	&elapsed($start_time, "\tElapsed", $verbose);
}


&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub get_common_variant_list {
	# Looks in the input file(s) to see which variants are above the threshold.  
	# If a variant is above the threshold frequency in at least one of the input files, 
	# it is retained in the list.  This will be all of the variants considered for the 
	# comparison.  This list will be used to analyze each of the input files.  Therefore
	# the time to analyze them can be calculated for one file and then multiplied by the 
	# number of input files.
	# Also enumerates every pair of variants to be compared in an AoA. 
	
	
	my $files = shift;
	my @files = @$files;	# list of variant files
	my $variants;		# Hashref of variants above the threshold in at least one of the replicates.   Format: $variants->{$type}->{$gene}->{$variant_combined_name}++;
	
	# Analyze each file separately
	for (my $i = 0; $i < @files; $i++){
		my $file = $files[$i];
					if ($verbose){	print STDERR "Reading file to get list of common variants: $file\n"; }
		my $readfh = open_to_read($file);
		while(<$readfh>){
			#	#Type	gene	position	consensus	variant	count	coverage	frequency	filter
			#	codon	AAM75158.1	159	AGN	AGT	424	19003	0.0223	CONSENSUS_N,
			chomp;
			next if (m/^#/);
			my ($type,$gene,$position,$consensus,$variant,$count,$coverage,$freq,$filter) = split(/\t/);
			unless($nofilter){
				next unless ($filter eq "PASS");
			}
			next unless ($freq > $variant_threshold);
			my $variant_combined_name = join ":", $position, $consensus, $variant; 		# position:consensus:variant
			$variants->{$type}->{$gene}->{$position}->{$variant_combined_name}++; 	
		}
		close ($readfh); 
	}

	# Did any variants pass the threshold in at least one file?  If not, then no need to continue.  
	unless ($variants){
		print STDERR "No variants passing threshold.\n";
		 exit;
	}

					if ($verbose){		print "variants\n", Dumper($variants); 	}
		
	# Estimate the time to run for each file (partly based on the size of the cleanreads.txt files)
	# Initialize variables for time estimation
	my $total_time = 0;
	my $total_comparisons = 0;
	my $counts;				# hashref of counts by type
	my $variant_positions; 	# hashref of counts by type and position.  I don't want to compare linkage for variants that have the same position, so I need to count the number of comparisons that would be done with variants of the same position and subtract that number from the total comparisons for each type.  


	foreach my $type (keys %$variants){
		foreach my $gene (keys %{$variants->{$type}}){
			foreach my $position (keys %{$variants->{$type}->{$gene}}){
				my $num_variants_at_position = scalar(keys %{$variants->{$type}->{$gene}->{$position}});
				
				# Increment the type and count the number of variants at a single position
				$counts->{$type} += $num_variants_at_position;
				$variant_positions->{$type}->{$position} += $num_variants_at_position;
			}	  			
		}
	}

						if ($verbose){		print "var positions\n", Dumper($variant_positions);	print "counts\n", Dumper($counts); 	}
	
	# Estimate time for each file by passing cleanreads for each file (and $counts, $variant_positions, which is the same for all files)
	for (my $i = 0; $i < @files; $i++){
		my $file = $files[$i];
		my $file_num = $i + 1;
		print STDERR "\nFile $file_num:\t$file\n";
		my ($time_this_file,$num_comparisons_this_file) = calculate_time_for_file($counts, $variant_positions, $clean_reads_files[$i]);		# $counts, $variant_positions will be same for all files, but need to pass the cleanreads file so it knows how many reads are in this file. 
		$total_time += $time_this_file; 
		$total_comparisons += $num_comparisons_this_file;
	}

	
	# Now check the total time it will take for all files to see if it is less than 24 hours. 
	my $calculated_time_str = hours_and_minutes($total_time);
	print STDERR "\n$total_comparisons total comparisons in all files.  Estimated to take about $calculated_time_str with $cpu cpu(s).\n\n";
	unless ($force){
		if ($total_time > 86400){	# 24 hours = 86400 seconds
			print STDERR "Too much time.  Please choose a different variant_threshold.\n";
			exit;	
		}
	}	

	
	# Now enumerate the comparisons that will be done for each variant type 
	my $comparisons; 	# HoHoAoA, first key type; second key gene; array of 2-way comparisons
	my $enumerated_comparisons; 
	foreach my $type (keys %$variants){ 		# Format: $variants->{$type}->{$gene}->{$position}->{$variant_combined_name}++; 	
		foreach my $gene (keys %{$variants->{$type}}){
			my @positions = sort {$a <=> $b} keys %{$variants->{$type}->{$gene}};
#								print STDERR "type: $type positions: ", "@positions", "\n"; 
			for (my $i = 0; $i < (@positions - 1); $i++){	# First variant position
#								print STDERR "i:$i\n";
				foreach my $first_var (keys %{$variants->{$type}->{$gene}->{$positions[$i]}}){		# Walk through all of the variants of the first position
#								print STDERR "\t\tfirst: $positions[$i] $first_var\n";
					for (my $j = $i + 1; $j < @positions; $j++){	# Second variant position.  start with one position higher than the first variant position 
#								print STDERR "\tj:$j\n";
						foreach my $second_var (keys %{$variants->{$type}->{$gene}->{$positions[$j]}}){	 	# walk through all of the variants of the second position. 
#								print STDERR "\t\t\tsecond: $positions[$j] $second_var\n";
							push @{$comparisons->{$type}->{$gene}}, [ $first_var, $second_var ]; 
							$enumerated_comparisons->{$type}++;
						}
					}
				}
			}			
		}
	}
	
	my $total_enumerated_comparisons = total(values %$enumerated_comparisons) * scalar(@files);
				if ($verbose){	print STDERR "Enumerated comparison count: $total_enumerated_comparisons\n", Dumper($enumerated_comparisons); 	}
	unless ($total_enumerated_comparisons == $total_comparisons){
		print STDERR "Number of comparisons don't match: $total_enumerated_comparisons (enumerated) vs. $total_comparisons (calculated)\n";
	}
	
	return ($variants,$comparisons);		# Maybe not necessary to return $variants...  it's used in read_variants right now to only save a variant to the hash if it's in this list of variants that pass the threshold.  Maybe not necessary; could save all, since it will be flushed.  
}
#-----------------------------------------------------------------------------
sub calculate_time_for_file {
	# Takes a hashref of counts (with type of variant as key and number of variants as value),
	# calculates the total number of comparisons to be performed, then adjusts the counts 
	# to remove self-comparisons (because they are not performed), then applies a formula
	# to get an estimate for time based on previous tests.
	# In the future, it might be a good idea to give an estimate based on the current performance,
	# e.g., as it is computing, get a rate of seconds per comparison and use the number of 
	# comparisons done and number remaining to get an estimate of time remaining and an 
	# adjusted total time as it is going (like in GATK).  

	# It might be simpler to just get rid of this sub altogether and use the enumerated_comparisons above for the calculation....  
	
	
	my ($counts, $variant_positions, $clean_reads_file) = @_;
	my %counts = %$counts;

	# Calculate the total number of comparisons by type (including self-comparisons, i.e., comparisons of variants at the same position)
	my %comparisons;
	foreach my $type (keys %counts){
		my $k = 2;
		if ($counts{$type} > 1){
			$comparisons{$type} = nCk($counts{$type}, $k);
		}
		else {
			$comparisons{$type} = 0;
		}
	}

	# Adjust the number of comparisons by type (in %comparisons) by deducting the comparisons between variants of the same position (because these will not be performed)
	foreach my $type (keys %$variant_positions){
		foreach my $pos (keys %{$variant_positions->{$type}}){
			my $this_count = $variant_positions->{$type}->{$pos};
			if ($this_count > 1){
				my $nCk = nCk($this_count,2);	# Count the number of self-comparisons that would be done for the variants at this position
				$comparisons{$type} -= $nCk;
						if ($verbose){		print STDERR "Removing $nCk comparison(s) for $type\n"; }
			}
		}
	}
						if ($verbose){		print Dumper(\%comparisons);				}

	my @grthn2;		# Array of variant types that have at least two positions.  
	my $total_comparisons = 0; 
	my $calculated_time = 0;
	if ($variant_positions){
		print STDERR "Variants passing threshold:\n";
		foreach my $type (sort keys %$variant_positions){
			my $num_variants = $counts{$type};
			my $num_positions = scalar (keys %{$variant_positions->{$type}});
			my $num_comparisons = $comparisons{$type};
			print STDERR "$type	$num_variants\tat $num_positions positions\t($num_comparisons comparisons)\n";
			push (@grthn2, $type) if ($num_positions >= 2);
		}
		if (@grthn2){
			print STDERR "Will compare linkage at the level of " . join " ", @grthn2; 
			print STDERR "\n";
			unless($table_only){		# No need to give a time estimate if we aren't going to compute p-values.				
				# Calculate time based on number of comparisons (range from ~15-25 seconds per comparison -- that probably depends a lot on the number of reads in the cleanreads.txt file...)
				my $read_count = `cat $clean_reads_file | wc -l` - 2;
								if($verbose){	print STDERR "read count: $read_count\n";	}
				foreach my $type (keys %comparisons){
					my $this_comparisons = $comparisons{$type};
					$total_comparisons += $this_comparisons;
					my $n = $counts{$type};
					my $effective_cpus = range([$n, $cpu], 'min'); 	# Will be different depending on the number of variants of this type (i.e., if you ask for 8 cpus, but there are only variants at 5 positions, you can only run 5 comparisons at a time -- 5 "effective" cpus).  ** Note that when I've changed the script to run each comparison in parallel (instead of each position in parallel), I'll need to change this estimation method.  It will be faster.  
#					my $this_time = $this_comparisons * (0.3099219 + $read_count * 0.0001317 + $effective_cpus * 0.1115561) / $effective_cpus;		# The equation calculates the estimated user time that will be used; need to divide by the number of cpus to get the wall time.  From 2014-01-07
#					my $this_time = $this_comparisons * (0.54984528 + $read_count * 8.6017e-6 + $effective_cpus * 8.96667e-9) / $effective_cpus; 	# From 2014-11-15 New equation using Excel Solver for all datapoints.
					my $this_time = $this_comparisons * (0.31 + $read_count * 4e-5) / $effective_cpus; 	# From 2014-11-15 just using all data from p 8, p 10, p 12 -> linear regression for x = read_count, y = time_per_comparison.  R^2 value 0.97566
					$calculated_time += $this_time;
									if($verbose){	print "effective_cpus: $effective_cpus\n";	}
									if($verbose){	print "Time for $type: $this_time\tRead counts: $read_count\n";	}
				}

			}
		}
		else {
			print STDERR "Not enough variants to make any comparisons.\n";
			 exit;
		}
	}
	else {
		print STDERR "No variants passing threshold.\n";
		 exit;
	}

									if($verbose){	print "Time for this file: $calculated_time\n";	}
	
	return ($calculated_time, $total_comparisons);

#	my ($high_range, $low_range) = (0,0);
#				$high_range += $nCk * 25 / $effective_cpus;		
#				$low_range 	+= $nCk * 8 / $effective_cpus; 
#			my ($high_range_str, $low_range_str) = ( hours_and_minutes($high_range), hours_and_minutes($low_range) );
#			print STDERR "$comparison_count total comparisons.  May take between $low_range_str and $high_range_str to complete with $cpu cpus.\n";


}
#-----------------------------------------------------------------------------
sub nCk {
	my ($n,$k) = @_;
	my $nCk = factorial($n)/(factorial($k)*factorial($n-$k));  	# n!/(k!*(n-k)!)
	return $nCk;
}
#-----------------------------------------------------------------------------
sub factorial {
	# From http://www.catonmat.net/blog/perl-one-liners-explained-part-three/
	use bignum;	# Changed from bigint to bignum so I can do float-point math with the results.
	my $num = shift;
	my $f = 1; 
	$f *= $_ for 1..$num; 
	return $f;
}
#-----------------------------------------------------------------------------
sub enumerate_comparisons {
	# Takes an arrayref of variants of format "position:consensus:variant"
	# Returns an AoA of all comparisons to perform.  Each sub-array will have two values
	my $variants = shift;
	
	my $comparisons;	# AoA to return
	for (my $i = 0; $i < @$variants; $i++){	# First variant
		for (my $j = $i + 1; $j < (@$variants); $j++){
			push @{$comparisons}, [ $variants->[$i], $variants->[$j] ]; 
		}
	}
	return $comparisons;
}
#-----------------------------------------------------------------------------
sub read_variants {
	my $file = shift;
	my $variants;		# HoHoA to store the variants being read from the file.  
	my $readfh = open_to_read($file);
	while(<$readfh>){
		#	#Type	position	consensus	variant	count	coverage	frequency	filter
		#	codon	AAM75158.1	159	AGN	AGT	424	19003	0.0223	CONSENSUS_N,
		chomp;
		next if (m/^#/);
		my ($type,$gene,$position,$consensus,$variant,$count,$coverage,$freq,$filter) = split(/\t/);
		unless($nofilter){
			next unless ($filter eq "PASS");
		}
		my $variant_combined_name = join ":", $position, $consensus, $variant; 		# position:consensus:variant
			#	print STDERR "$variant_combined_name\n";
		next unless (exists($variants_to_consider->{$type}->{$gene}->{$position}->{$variant_combined_name})); 		# Check to see that it was above threshold in at least one of the replicates
			#	print STDERR "\tpass\n";
		# Store the variant info in a hash
		my %hash = (
			'position' => $position,
			'variant' => $variant,
			'consensus' => $consensus,
			'frequency' => $freq,
		);
		
		# Add variant to the array.
#		push @{$variants->{$type}->{$gene}}, \%hash;
		
		# Add variant to the hash.  
		$variants->{$type}->{$gene}->{$variant_combined_name} = \%hash;
		
		
	}
	close ($readfh);
	
	return $variants;
}
#-----------------------------------------------------------------------------
sub get_reads {
	my $file = shift;
	my $reads; 	# Save the reads in an AoH.  Hashes will have keys 'gene', 'pos', 'seq'
	my $readfh = open_to_read($file);
	while(<$readfh>){
		#	#readID:gene:positionForFirstNucleotideBase	cleanRead
		#	Y3YRC:ATAAGCTT:3|AAM75158.1|158	AGCAGNTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGNCCTTGTACTGTGGGGTATTCATCACCCGCCTAACAGTAA
		#	#readID:gene:positionForFirstAminoAcid	cleanPeptide
		#	Y3YRC:ATAAGCTT:3|AAM75158.1|472	SXFYRNLLWLTKKESSYPELKNSYVNKKRKEXLVLWGIHHPPNS
		# OR
		#	#readID|gene|positionForFirstNucleotideBase	cleanRead
		#	1-42|Seq12_093009_HA_cds|382 TTCGAAAGATTCAAAATATTTCCCAAAGAAAGCTCATGGCCCGACCACAACACAAACGGAGTAACGGCAGCATGCTCCCATGAGGGGAAAAACAGTTTTTACAGAAATTTGCTATGGCTGACGAAGAAGGAGAGCTCATACCCAGAGCTGAAAAATTCTTATGTGAACAAAAAAAGGAAAGAAGTCCTTGTACTGTGGGGTATTCATCNNNNNNNTAACAGTAAGGAACAACAGAATCTCTATCAGAATGAAAATGCTTATGTCTCTGTAGTGACTTCAAATTATAACAGGAGATTTACCCCGGAAATAGCAGAAAGACCCAAAGTAAAAGGTCAAGCTGGGAGGATGAACTATTACTGGACCTTGCTAAAACCCGGAGACACAATAATATTTGAGGCAAATGGAAATCTAATAGCACCAATGTATGCTTTC
		#	#readID|gene|positionForFirstAminoAcid	cleanPeptide
		#	1-42|Seq12_093009_HA_cds|128 FERFKIFPKESSWPDHNTNGVTAACSHEGKNSFYRNLLWLTKKESSYPELKNSYVNKKRKEVLVLWGIHXXXNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVKGQAGRMNYYWTLLKPGDTIIFEANGNLIAPMYAF
	
		chomp;
		next if (m/^#/);
		next if (m/^\s+\d+\s+\d+/);		# First line of file, prepended as a PHYLIP header
		my ($name,$read) = split(/\s+/);
		my @name = split(/\|/,$name);	# Was split(/_/,$name);
		my $position = pop(@name);
		my $gene = pop(@name);
		
		# Get sequence count from id.  If the id has format number - number, then it will be assumed the second number is the read count.  
		my $id = shift(@name);
		my ($num,$seq_count) = (1,1);
		if ($id =~ m/^\d+-\d+$/){
			($num,$seq_count) = split(/-/, $id);
		}

		my %hash = (
			'gene' => $gene,
			'position' => $position,
			'seq' => $read,
			'count' => $seq_count,
		);
		push @$reads, \%hash;
	}

	close ($readfh);
	return $reads;
}
#-----------------------------------------------------------------------------
sub linkage_disequilibrium {
	my ($reads, $variants, $type, $label, $linkage_disequil_fh) = @_;  
		# $reads is a AoH where hashes have keys 'gene', 'pos', 'seq'
		# $type will decide how to split up the read/sequence (either into single characters or multiple
	
	print STDERR "Processing variants of type '$type'\n"; 
	
	my $this_variants = $variants->{$type};		#HoAoH, where first key is gene, and keys of deepest hash are 'position', 'variant', 'consensus'
	
	my $this_variants_to_consider = $variants_to_consider->{$type}; 		# Format of $variants_to_consider : $variants->{$type}->{$gene}->{$variant_combined_name}++;
	
	my $maxProcs = $cpu;
	
	my $pl = Parallel::Loops->new($maxProcs);
	
	#					 print Dumper($this_variants); 		# exit;
	my @stats; 	# AoA to save data.  Need to save it all first so I can then calculate the FDR from the p-values. 
	
	$pl->share(\@stats);  # Make the array shared between fork children
	 
	foreach my $gene (keys %$this_variants_to_consider){
		# Walk through each variant doing all pairwise comparisons 
		my @variants = sort keys %{$this_variants_to_consider->{$gene}}; 
		my $count = 0;
		my $comparison = "";	# initialize.  Single array of two elements.  	
		$pl->while( sub {$count++; $comparison = $comparisons->{$type}->{$gene}->[$count - 1]; },	sub {			# was $pl->foreach( \@primerIDs,	sub {	
			my @variants_to_compare = @$comparison;
	#								if($verbose){		print Dumper(\@variants_to_compare);	}
			my ($first_pos,$first_cons,$first_var) 		= split(/:/, $comparison->[0]);
			my ($second_pos,$second_cons,$second_var) 	= split(/:/, $comparison->[1]);
			my ($first_var_freq,$second_var_freq) = (0,0);	# Default freq if not found in the $this_variants hash
			
			$first_var_freq  = $this_variants->{$gene}->{$comparison->[0]}->{frequency} if (exists($this_variants->{$gene}->{$comparison->[0]}));
			$second_var_freq = $this_variants->{$gene}->{$comparison->[1]}->{frequency} if (exists($this_variants->{$gene}->{$comparison->[1]}));
			
			
			next if ($first_pos == $second_pos);		# This shouldn't ever happen, but just in case.
			
			# Now walk through the $reads AoH and count up the variants at each of these positions
			my $AB = 0;		# consensus at both positions. A is first position, B is second position.
			my $Ab = 0;		# consensus at first position, variant at second position
			my $aB = 0; 	# variant at first position, consensus at second position
			my $ab = 0; 	# variant at both positions
			foreach my $read (@$reads){
#							 print calculate_linkage_disequilibrium.pl Convert_reads_variants_minfreq0_codon.xls($read);	print "gene: $gene\n";
				next unless ($gene eq $read->{'gene'});		# Only consider a read if it was mapped to the gene that the variant came from, otherwise the positions won't be accurate.  
				my ($read_pos,$read_seq,$read_count) = ($read->{'position'},uc($read->{'seq'}),$read->{'count'});
				if ($type eq 'codon'){
					$read_pos = ($read_pos + 2) / 3;		# Convert nucleotide position to amino acid/codon position.  (Add 2 before dividing by 3 because the position refers to the first base of the codon, so dividing that by 3 will not give a whole number.)
				}
						# print "read: $read->{'gene'}, $read_pos, $read_seq\n";
				my @seq;
				if ($type =~ m/nuc|aa/){
					@seq = split(/|/,$read_seq);		# nuc or aa, split into single characters
				}
				else {
					@seq = unpack("(a3)*", $read_seq);	# codon, split into triplets
				}
						# print join " ", @seq, "\n";
				my $first_pos_in_read = $first_pos - $read_pos;
				my $second_pos_in_read = $second_pos - $read_pos;
						# print "positions in read: $first_pos_in_read, $second_pos_in_read\n";	
				next unless($seq[$first_pos_in_read] && $seq[$second_pos_in_read]);
			
				unless ($seq[$first_pos_in_read] && $first_cons && $seq[$second_pos_in_read] && $second_cons && $first_var && $second_var){
					print STDERR "One of the positions or variants is not right:\n";
					print STDERR join ", ", $first_pos_in_read, $seq[$first_pos_in_read], $second_pos_in_read, $seq[$second_pos_in_read], $first_cons, $second_cons, $first_var, $second_var; print STDERR "\n$read_seq\n"; exit;
				}
				if ( ($seq[$first_pos_in_read] eq $first_cons) && ($seq[$second_pos_in_read] eq $second_cons) ){
					$AB += $read_count;
				}
				elsif( ($seq[$first_pos_in_read] eq $first_cons) && ($seq[$second_pos_in_read] eq $second_var) ){
					$Ab += $read_count;
				}
				elsif( ($seq[$first_pos_in_read] eq $first_var) && ($seq[$second_pos_in_read] eq $second_cons) ){
					$aB += $read_count;
				}
				elsif( ($seq[$first_pos_in_read] eq $first_var) && ($seq[$second_pos_in_read] eq $second_var) ){
					$ab += $read_count;
				}
				else {
					# This is possible to happen, I suppose.  Doesn't match these variants or consensus.  Don't count it. 
					my $first_pos_in_read_1_based = $first_pos_in_read + 1;
					my $second_pos_in_read_1_based = $second_pos_in_read + 1;
					print STDERR "Either variants ($first_pos: $first_var, $second_pos: $second_var) or consensus ($first_pos: $first_cons, $second_pos: $second_cons) doesn't match, so not counting. Actual is $first_pos ($first_pos_in_read_1_based in read): $seq[$first_pos_in_read], $second_pos ($second_pos_in_read_1_based in read): $seq[$second_pos_in_read].  Probably a read with a variant that was below your threshold or else a position with more than 2 variants.  Or, you have multiple samples and this variant is present in one sample but not the other.\n", Dumper($read) if ($verbose);
				}
			}
						# print "AB: $AB, Ab: $Ab, aB: $aB, ab: $ab\n"; 
			
			my $comparison_combined_name = join ":", $first_pos, $first_cons, $first_var, $second_pos, $second_cons, $second_var; 		# position:consensus:variant
			#group	sample	type	gene	comparison	pos1	c1	v1	v1freq	pos2	c2	v2	v2freq	c1c2	c1v2	v1c2	v1v2	p-value	OR	FDR
				# (p-value	OR	FDR) columns will be added outside of the parallel loop
			my @array = ( $group_id, $label, $type, $gene, $comparison_combined_name, $first_pos, $first_cons, $first_var, $first_var_freq, $second_pos, $second_cons, $second_var, $second_var_freq, $AB, $Ab, $aB, $ab );
			push @stats, \@array;

			
			if ($count % 100 == 0){
				print STDERR "Processed $count comparisons. "; 
				&elapsed($start_time, 'Elapsed', 1);   # $verbose
			}				
		});		# End of parallel loop
	}	
	

	my $stats = \@stats;

	# Compute p-values from contingency table counts
	my $R = Statistics::R->new();		# http://search.cpan.org/~fangly/Statistics-R-0.34/lib/Statistics/R.pm
	my ($p,$OR);
	for (my $i = 0; $i < @$stats; $i++){
		if ($table_only){
			# Don't compute p-values
			($p,$OR) = ("NA","NA");
		}
		else {
			# Now get the p-value using Fisher's exact test in R	
			my ($AB,$Ab,$aB,$ab) = @{$stats->[$i]}[-4..-1];
			$R->set('x', [$AB,$Ab,$aB,$ab]);
			$R->run(q`y = matrix(x,nrow=2)`);
			$R->run(q`temp = fisher.test(y)`);
			$R->run(q`p = temp$p.value`);
			$R->run(q`or = temp$estimate`);
			$OR = $R->get('or');
			$p = $R->get('p');
		}
		push @{$stats->[$i]}, ($p,$OR);
	}



	# Now get the FDR-adjusted p-value using Benjamini & Hochberg method
							# print Dumper($stats);
	my @pvalues;
	foreach (@$stats){
		push @pvalues, $_->[-2];
	}
							# print join " ", @pvalues; print "\n";
	if (scalar(@pvalues)>1){
		$R->set('x',\@pvalues);		
		$R->run(q`y=p.adjust(x,method="BH",n=length(x))`);	# 	z$adj.p.value = p.adjust(z$p.value, method="BH",n = sum(!is.na(z$p.value)))
		my $fdr_array = $R->get('y');
								# print join " ", @$fdr_array; print "\n";
		die "Error in calculating FDR adjusted p-values: different number of p-values and FDR values\n", join " ", @pvalues, "\n", @$fdr_array, "\n" unless (scalar(@$fdr_array) == scalar(@$stats));
		for (my $i=0; $i < @$fdr_array; $i++){
			push @{$stats->[$i]}, $fdr_array->[$i]; 
		}
	}
	else {
		# There is only one item in @$stats array
		# Put the same value for FDR as for p-value.  
		push @{$stats->[0]}, $stats->[0]->[-2];
	}
	
	# Print out the results
	foreach (@$stats){
		if (scalar(@{$_}) > 1){
			print $linkage_disequil_fh join "\t", @{$_};
			print $linkage_disequil_fh "\n";
		}
	}
	
	
	
	return $stats;		# Not really necessary since this variants have already been printed, but might be nice to return this value for debugging purposes.
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

