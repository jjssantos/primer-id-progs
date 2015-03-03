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

my $save;
my $files;
my $verbose;
my $output;
my $gzip;
my $prefix;
my $nofilter;
my $variant_threshold;
my $cpu;
my $force;
GetOptions('save=s' => \$save, 'output=s' => \$output, 'verbose' => \$verbose, 'files=s' => \$files, 'gzip' => \$gzip, 'prefix=s' => \$prefix, 'nofilter' => \$nofilter, 'variant_threshold=s' => \$variant_threshold, 'cpu=s' => \$cpu, 'p=s' => \$cpu, 'force' => \$force);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "
calculate_linkage_disequilibrium.pl Takes a file of variants and sequences and checks for 
linkage between each pair of variants.  Variants are queried at the level of nucleotide, 
codon, and amino acid.  Linkage is computed with fisher's exact test and is reported as 
p-value, FDR-adjusted p-value, and Odds ratio (OR).  This script requires R on the PATH.  

Usage: calculate_linkage_disequilibrium.pl <variants.xls> <cleanreads.txt> <cleanpeptides.txt>

variants.xls, cleanreads.txt, and cleanpeptides.txt are output files from 
convert_reads_to_amino_acid.pl
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
#type	gene	pos1	c1	v1	pos2	c2	v2	c1c2	c1v2	v1c2	v1v2	p-value	OR	FDR
nuc	AAM75158.1	526	G	A	556	A	G	18329	255	89	32	8.622218e-31	25.82045	8.622218e-31
codon	AAM75158.1	176	GAG	AAG	186	AGG	GGG	18315	255	88	32	6.588369e-31	26.09373	6.588369e-31
aa	AAM75158.1	176	E	K	186	R	G	18317	255	88	32	6.566664e-31	26.09658	6.566664e-31

type is the type of variant being considered, either nuc, codon, or aa.
c1 and v1 refer to consensus and variant, respectively, in the first position
c2 and v2 refer to consensus and variant, respectively, in the second position
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

OPTIONS:
	--save		Directory in which to save files. Default = pwd.  If folder doesn't exist, 
			it will be created.
	--prefix	Prefix for output files.  Default = Calculate_linkage
	--variant_threshold	Minimum threshold to consider for a variant.  Default is 0.001.     
	--nofilter	Compare all variants in variants.xls file.  Default is to only consider 
			variants that have PASS status.  
	-p/--cpu		Number of CPUs to use for parallel processing.  Default = 2.
	--force		Force the analysis to run, even if expected runtime is > 24h.  Default is
			to only run if low range of expected is <= 24h. 
";

# To do
# Maybe add another option to specify which filters to allow (comma-delimited)
# Add threads to 
# Compute quasi-cliques with DENSE.  (or Cocain?).  Add description of requirements (DENSE/mgqce) to usage.
# Add fdr threshold option for quasi-clique variants. Default = 0.05
# Add part to count the number of pairs with fdr <= 0.05 while printing.
	
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

unless ($ARGV[1]){	print STDERR "$usage\n";	exit;	}
unless ($ARGV[2]){	print STDERR "No peptide file found.  Will only look at nucleotide and codon variants.\n";	}

my $save_dir = $save || Cwd::cwd();
unless (-d $save_dir){	mkdir($save_dir) or warn "$!\n"	}
warn "save directory: $save_dir\n";

$cpu ||= 2;

my $start_time = time;

$variant_threshold //= 0.001;

my $variants = read_variants($ARGV[0]);		# HoAoH; first keys, type, i.e., 'nuc', 'codon', 'aa'; array of variants; hash with keys 'gene', 'pos', 'var', 'consensus', 'frequency'	 (This isn't quite accurate, I think I changed the structure a little since I first wrote this comment.)
#						print Dumper($variants);


$prefix ||= "Calculate_linkage";
my $linkage_disequil	= $save_dir . "/" . $prefix . ".linkage.minfreq".$variant_threshold.".xls";			# was $prefix . "_linkage_disequilibrium_report_minfreq".$variant_threshold.".xls";
$linkage_disequil =~ s/_linkage.linkage/_linkage/;
my $linkage_disequil_fh		= open_to_write("$linkage_disequil");
print $linkage_disequil_fh "#type\tgene\tpos1\tc1\tv1\tv1freq\tpos2\tc2\tv2\tv2freq\tc1c2\tc1v2\tv1c2\tv1v2\tp-value\tOR\tFDR\n";		


# Now read the nucleotide file *cleanreads.txt into memory and get linkage for variants of type 'nuc' and 'codon'
my $clean_reads = get_reads($ARGV[1]);
my $stats = linkage_disequilibrium($clean_reads,'nuc');
$stats = linkage_disequilibrium($clean_reads,'codon');	# I'll probably need a special step for codon

# Flush the variable to free up memory
undef($clean_reads);

# Now, if a cleanpeptide.txt file is available, read it into memory 
exit unless $ARGV[2];
$clean_reads = get_reads($ARGV[2]);
$stats = linkage_disequilibrium($clean_reads,'aa');		# Calculates pairwise linkage between variants.  


&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub read_variants {
	my $file = shift;
	my $variants;
	my $readfh = open_to_read($file);
	my %count;				# hash of counts by type
	my $variant_positions; 	# hashref of counts by type and position.  I don't want to compare variants that have the same position, so I need to count the number of comparisons that would be done with variants of the same position and subtract that number from the total comparisons for each type.  
	while(<$readfh>){
		#	#Type	position	consensus	variant	count	coverage	frequency	filter
		#	codon	AAM75158.1	159	AGN	AGT	424	19003	0.0223	CONSENSUS_N,
		chomp;
		next if (m/^#/);
		my ($type,$gene,$position,$consensus,$variant,$count,$coverage,$freq,$filter) = split(/\t/);
		unless($nofilter){
			next unless ($filter eq "PASS");
		}
		next unless ($freq > $variant_threshold);
		my %hash = (
			'position' => $position,
			'variant' => $variant,
			'consensus' => $consensus,
			'frequency' => $freq,
		);
		push @{$variants->{$type}->{$gene}}, \%hash;
		$count{$type}++;
		$variant_positions->{$type}->{$position}++;
	}
	close ($readfh);
						if ($verbose){		print Dumper($variant_positions);	}
	# Adjust the number of comparisons by type (in %count) by deducting the comparisons between variants of the same position
	foreach my $type (keys %$variant_positions){
		foreach my $pos (keys %{$variant_positions->{$type}}){
			my $count = $variant_positions->{$type}->{$pos};
			if ($count > 1){
				my $nCk = nCk($count,2);	# Count the number of comparisons that would be done for these variants.
				$count{$type} -= $nCk;
			}
		}
	}
						if ($verbose){		print Dumper(\%count);				}
	
	my @grthn2;
	if (%count){
		print STDERR "Variants passing threshold:\n";
		foreach my $type (sort keys %count){
			print STDERR "$type	$count{$type}\n";
			push (@grthn2, $type) if ($count{$type} >= 2);
		}
		if (@grthn2){
			print STDERR "Will compare linkage at the level of " . join " ", @grthn2; 
			print STDERR "\n";
			# Calculate the number of comparisons to estimate the time it should take (range from ~15-25 seconds per comparison -- that probably depends a lot on the number of reads in the cleanreads.txt file...)
			my $comparison_count = 0;
			my $read_count = `cat $ARGV[1] | wc -l` - 2;
								if($verbose){	print STDERR "read count: $read_count\n";	}
			my ($high_range, $low_range) = (0,0);
			my $calculated_time = 0;
			foreach my $type (@grthn2){
				my $k = 2;
				my $n = $count{$type};
				my $nCk = nCk($n, $k);				
				$comparison_count += $nCk;
				my $effective_cpus = range([$n, $cpu], 'min'); 
								if($verbose){	print "effective_cpus: $effective_cpus\n";	}
#				$high_range += $nCk * 25 / $effective_cpus;		
#				$low_range 	+= $nCk * 8 / $effective_cpus; 
				$calculated_time += $nCk * (0.3099219 + $read_count * 0.0001317 + $effective_cpus * 0.1115561) / $effective_cpus;		# The equation calculates the estimated user time that will be used; need to divide by the number of cpus to get the wall time.  
			}
#			my ($high_range_str, $low_range_str) = ( hours_and_minutes($high_range), hours_and_minutes($low_range) );
#			print STDERR "$comparison_count total comparisons.  May take between $low_range_str and $high_range_str to complete with $cpu cpus.\n";
			my $calculated_time_str = hours_and_minutes($calculated_time);
			print STDERR "$comparison_count total comparisons.  Estimated to take about $calculated_time_str with $cpu cpu(s).\n";
			unless ($force){
				if ($calculated_time > 86400){	# 24 hours = 86400 seconds
					print STDERR "Too much time.  Please choose a different variant_threshold.\n";
					exit;	
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
	my ($reads, $type) = @_;  
		# $reads is a AoH where hashes have keys 'gene', 'pos', 'seq'
		# $type will decide how to split up the read/sequence (either into single characters or multiple
	my $this_variants = $variants->{$type};		#HoAoH, where first key is gene, and keys of deepest hash are 'position', 'variant', 'consensus'
	
	my $maxProcs = $cpu;
	
	my $pl = Parallel::Loops->new($maxProcs);
	
	#					 print Dumper($this_variants); 		# exit;
	my @stats; 	# AoA to save data.  Need to save it all first so I can then calculate the FDR from the p-values. 
	
	$pl->share(\@stats);  # Make the array shared between fork children
#	$pl->share($reads);
	 
	foreach my $gene (keys %$this_variants){
		# Walk through each variant doing all pairwise comparisons [Maybe this could be threaded to speed things up...]

#		my $i = -1;		
#		$pl->while( sub { $i++ < ( @{$this_variants->{$gene}}) }, sub {	# First variant
#			for (my $j = $i+1; $j < @{$this_variants->{$gene}}; $j++){	# Second variant
		for (my $i = 0; $i < @{$this_variants->{$gene}}; $i++){	# First variant
			my $j = $i;
			$pl->while( sub { $j++ < ( @{$this_variants->{$gene}} - 1) }, sub {	# Second variant.  Compare first variant to all of these second variants in parallel.  I tried making the first loop parallel but that doesn't seem to speed it up at all in my quick tests, and the output is messier than this.
				my ($first_pos,$first_cons,$first_var,$first_var_freq) 		= ($this_variants->{$gene}->[$i]->{'position'},$this_variants->{$gene}->[$i]->{'consensus'},$this_variants->{$gene}->[$i]->{'variant'},$this_variants->{$gene}->[$i]->{'frequency'});
				my ($second_pos,$second_cons,$second_var,$second_var_freq) 	= ($this_variants->{$gene}->[$j]->{'position'},$this_variants->{$gene}->[$j]->{'consensus'},$this_variants->{$gene}->[$j]->{'variant'},$this_variants->{$gene}->[$j]->{'frequency'});
							if ($verbose){	print "i: $i\t j: $j\n";	}
							if ($verbose){	print "first: $first_pos, $first_cons/$first_var ($first_var_freq)\tsecond: $second_pos, $second_cons/$second_var ($second_var_freq)\n";		}
				unless ($first_pos == $second_pos){
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
							# Doesn't match these variants or consensus.  Don't count it.  
						}
					}
								# print "AB: $AB, Ab: $Ab, aB: $aB, ab: $ab\n"; 
					# Now get the p-value using Fisher's exact test in R
					my $R = Statistics::R->new();		# http://search.cpan.org/~fangly/Statistics-R-0.31/lib/Statistics/R.pm
					$R->set('x', [$AB,$Ab,$aB,$ab]);
					$R->run(q`y = matrix(x,nrow=2)`);
					$R->run(q`temp = fisher.test(y)`);
					$R->run(q`p = temp$p.value`);
					$R->run(q`or = temp$estimate`);
					my $OR = $R->get('or');
					my $p = $R->get('p');
								# print "p-value: $p\tOddsRatio: $OR\n";
					#type	gene	pos1	c1	v1	v1freq	pos2	c2	v2	v2freq	c1c2	c1v2	v1c2	v1v2	p-value	OR	FDR
					my @array = ( $type, $gene, $first_pos, $first_cons, $first_var, $first_var_freq, $second_pos, $second_cons, $second_var, $second_var_freq, $AB, $Ab, $aB, $ab, $p, $OR );
					push @stats, \@array;
				}
			});		# End of parallel loop
		}
	}	
	# Now get the FDR-adjusted p-value using Benjamini & Hochberg method
							# print Dumper($stats);
	my @pvalues;
	my $stats = \@stats;
	foreach (@$stats){
		push @pvalues, $_->[-2];
	}
							# print join " ", @pvalues; print "\n";
	if (scalar(@pvalues)>1){
		my $R = Statistics::R->new();
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
	
	
	
#	return $stats;
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
sub nCk {
	my ($n,$k) = @_;
	my $nCk = factorial($n)/(factorial($k)*factorial($n-$k));  	# n!/(k!*(n-k)!)
	return $nCk;
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub read_and_write {
	my $file = shift;
	my $readfh = open_to_read($file);
	my ($filename,$dir,$ext) = fileparse($file,@SUFFIXES);		# fileparse($file,qr/\.[^.]*/);
#	my $newfilename = $save_dir.'/'.$filename.'__'.$variable.'.bed';
#	my $writefh = open_to_write($newfilename, $gzip);
	while (<$readfh>){
		
	}
	
#	close ($writefh);
	close($readfh);
}
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
