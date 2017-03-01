#-------------------------------------------------------------------------------
#----                                MODULE NAME                            ----
#-------------------------------------------------------------------------------
package primerid;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use IO::File;
use IO::Zlib; 		#qw(:gzip_external 1);
use Carp;
use Data::Dumper;
use Cwd;
use File::Basename;
use Array::Utils qw(:all);
use Statistics::R;

# Custom
use aomisc;

our @ISA = qw(Exporter);
our @EXPORT = qw(
	make_codon_to_aa_hash
	print_reports
	print_merged_report
	print_variants
);


# Andrew J. Oler, PhD
# Created at
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

#-------------------------------------------------------------------------------
#----------------------------------- CHANGELOG ---------------------------------
#-------------------------------------------------------------------------------
# 2017-02-05
# Created module.  Moved some subroutines from convert_reads_to_amino_acid.pl to 
# this module so they can be used by merge_tally.pl as well.
# 2017-02-11
# Removed $save_dir as input, using $prefix alone ($prefix can contain absolute 
# or relative path).
# 2017-02-15
# Added cMAF cMAF_95%_CI_low cMAF_95%_CI_high medianUnambigFreq 
# median95%ConfIntHigh passThreshold columns to the output tally files for nuc, 
# codon, aa.  Added 95%_CI_low 95%_CI_high columns to variant report.  Moved 
# some code from print_reports to a new sub get_unambig_values, and added new
# subs get_confidence_interval and get_median_unambig_freq.  



#-------------------------------------------------------------------------------
#----------------------------------- FUNCTIONS ---------------------------------
#-------------------------------------------------------------------------------
sub make_codon_to_aa_hash {
	# This sub takes no input.  
	# Returns a hash where keys are codon sequences and values are amino acid sequences.

	my %converter = (		# Modified from http://www.wellho.net/resources/ex.php4?item=p212/3to3
	    'TCA' => 'S', # Serine
	    'TCC' => 'S', # Serine
	    'TCG' => 'S', # Serine
	    'TCT' => 'S', # Serine
	    'TCN' => 'S', # Serine wobble
	    'TTC' => 'F', # Phenylalanine
	    'TTT' => 'F', # Phenylalanine
	    'TTY' => 'F', # Phenylalanine ambiguous C or T
	    'TTA' => 'L', # Leucine
	    'TTG' => 'L', # Leucine
	    'TTR' => 'L', # Leucine ambiguous A or G
	    'TAC' => 'Y', # Tyrosine
	    'TAT' => 'Y', # Tyrosine
	    'TAY' => 'Y', # Tyrosine ambiguous C or T
	    'TAA' => '*', # Stop
	    'TAG' => '*', # Stop
	    'TGC' => 'C', # Cysteine
	    'TGT' => 'C', # Cysteine
	    'TGY' => 'C', # Cysteine ambiguous C or T
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
	    'CAY' => 'H', # Histidine ambiguous C or T
	    'CAA' => 'Q', # Glutamine
	    'CAG' => 'Q', # Glutamine
	    'CAR' => 'Q', # Glutamine ambiguous A or G
	    'CGA' => 'R', # Arginine
	    'CGC' => 'R', # Arginine
	    'CGG' => 'R', # Arginine
	    'CGT' => 'R', # Arginine
	    'CGN' => 'R', # Arginine wobble
	    'ATA' => 'I', # Isoleucine
	    'ATC' => 'I', # Isoleucine
	    'ATT' => 'I', # Isoleucine
	    'ATH' => 'I', # Isoleucine ambiguous A or C or T
	    'ATM' => 'I', # Isoleucine ambiguous A or C
	    'ATW' => 'I', # Isoleucine ambiguous A or T
	    'ATY' => 'I', # Isoleucine ambiguous C or T
	    'ATG' => 'M', # Methionine
	    'ACA' => 'T', # Threonine
	    'ACC' => 'T', # Threonine
	    'ACG' => 'T', # Threonine
	    'ACT' => 'T', # Threonine
	    'ACN' => 'T', # Threonine wobble
	    'AAC' => 'N', # Asparagine
	    'AAT' => 'N', # Asparagine
	    'AAY' => 'N', # Asparagine ambiguous C or T
	    'AAA' => 'K', # Lysine
	    'AAG' => 'K', # Lysine
	    'AAR' => 'K', # Lysine ambiguous A or G
	    'AGC' => 'S', # Serine
	    'AGT' => 'S', # Serine
	    'AGY' => 'S', # Serine ambiguous C or T
	    'AGA' => 'R', # Arginine
	    'AGG' => 'R', # Arginine
	    'AGR' => 'R', # Arginine ambiguous A or G
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
	    'GAY' => 'D', # Aspartic Acid ambiguous C or T
	    'GAA' => 'E', # Glutamic Acid
	    'GAG' => 'E', # Glutamic Acid
	    'GAR' => 'E', # Glutamic Acid ambigous A or G
	    'GGA' => 'G', # Glycine
	    'GGC' => 'G', # Glycine
	    'GGG' => 'G', # Glycine
	    'GGT' => 'G', # Glycine
	    'GGN' => 'G', # Glycine wobble
	    );
	# Add constituent ambiguous characters in place of N as well.  
	my @iupac_ambig_nuc = (qw(R Y S W K M B D H V));	# http://www.bioinformatics.org/sms/iupac.html
	foreach my $codon (keys %converter){
		if ($codon =~ m/N$/){
			foreach my $ambig (@iupac_ambig_nuc){
				my $new_codon = $codon;
				$new_codon =~ s/N$/$ambig/;	# Replace N for the ambiguous nucleotide
				$converter{$new_codon} = $converter{$codon};	# Assign the amino acid value to the same aa as the codon ending in N
			}
		}
	}

	return %converter;
}
#-----------------------------------------------------------------------------
sub print_reports {
	my ($nuc_aa_codon_tally, $label, $gene, $NUC, $CODON, $AA, $cds, $converter, $prefix, $start_time, $verbose, $debug) = @_;
	# Print Reports to filehandles $nucleotide_report_fh, $codon_report_fh, and $amino_acid_report_fh
	# TO DO: Add more documentation about the input data format ($cds, $nuc_aa_codon_tally, etc.)

	#name(chr:gene:aminoAcidPosition:consensusAminoAcid)	chr	gene	aminoAcidPosition	refAminoAcid	consensusAminoAcid	Sample	coverageDepth	numA	numC	numG	numT	numN	[del	SumACGT	numNonzeroACGT	Ref	Nonref]

	# To make things a little less redundant, I could put the medianUnambigFreq and median95%ConfIntHigh in the header (rather than having it in every single row...), e.g, in a line starting with ##, but then I would need to do more modifications in the other scripts to skip lines starting with ##, especially scripts that use the get_header subroutine...
	print STDERR "Printing out frequency tables. ";
	&elapsed($start_time, 'Elapsed', $verbose);


	my $nucleotide_report 	= $prefix . ".nuc.tally.xls";
	my $codon_report 		= $prefix . ".codon.tally.xls";
	my $amino_acid_report 	= $prefix . ".aa.tally.xls";

	my $nucleotide_report_fh 	= open_to_write("$nucleotide_report");
	my $codon_report_fh 		= open_to_write("$codon_report");
	my $amino_acid_report_fh 	= open_to_write("$amino_acid_report");

	print $nucleotide_report_fh "#name\t";
	print $nucleotide_report_fh join "\t", qw(gene nucleotidePosition refNucleotide consensusNucleotide refDiffConsensus Sample coverageDepth unambigCoverageDepth unambigConsensus numUnambigConsensus numUnambigNonConsensus majorAltAllele numMajorAltAllele cMAF cMAF_95%_CI_low cMAF_95%_CI_high medianUnambigFreq median95%ConfIntHigh passThreshold);
	print $nucleotide_report_fh "\tnum"; 
	print $nucleotide_report_fh join "\tnum", @$NUC, "Other";
	# print $nucleotide_report_fh "\t"; 
	# print $nucleotide_report_fh join "\t", qw(del SumACGT numNonzeroACGT Ref Nonref percentNonref);
	print $nucleotide_report_fh "\n";

	print $codon_report_fh "#name\t";
	print $codon_report_fh join "\t", qw(gene codonPosition refCodon consensusCodon refDiffConsensus Sample coverageDepth unambigCoverageDepth unambigConsensus numUnambigConsensus numUnambigNonConsensus majorAltAllele numMajorAltAllele cMAF cMAF_95%_CI_low cMAF_95%_CI_high medianUnambigFreq median95%ConfIntHigh passThreshold);		#  unambigCoverageDepth numConsensus numNonConsensus numMajorAltAllele majorAltAllele
	print $codon_report_fh "\tnum"; 
	print $codon_report_fh join "\tnum", @$CODON, "Other";
	# print $codon_report_fh "\t"; 
	# print $codon_report_fh join "\t", qw(del SumACGT numNonzeroACGT Ref Nonref percentNonref);		# Maybe syn, nonsyn; or maybe save that for the summary report.  
	print $codon_report_fh "\n";

	print $amino_acid_report_fh "#name\t";
	print $amino_acid_report_fh join "\t", qw(gene aminoAcidPosition refAminoAcid consensusAminoAcid refDiffConsensus Sample coverageDepth unambigCoverageDepth unambigConsensus numUnambigConsensus numUnambigNonConsensus majorAltAllele numMajorAltAllele cMAF cMAF_95%_CI_low cMAF_95%_CI_high medianUnambigFreq median95%ConfIntHigh passThreshold); 
	print $amino_acid_report_fh "\tnum"; 
	print $amino_acid_report_fh join "\tnum", @$AA, "Other";  		# Needed to add "num" before amino acids, since R converts _ to "X_" and * to "X."
	# print $amino_acid_report_fh "\t"; 
	# print $amino_acid_report_fh join "\t", qw(SumKnownAminoAcids numNonzeroKnownAminoAcids Ref Nonref percentNonref);
	print $amino_acid_report_fh "\n"; 



	my @all_codons = @$CODON;
	for (@all_codons){
		s/:.+//;		# Was s/:\w+//;  changed to s/:.+//; for numTAA:* and others with *
	}
	
	my $R = Statistics::R->new();		# Open a bridge to R, to pass to get_median_unambig_freq or get_confidence_int sub. 
	
	foreach my $type (qw(nuc codon aa)){		#e.g., $nuc_aa_codon_tally->{'aa'}->{$pos}->{$aa}++;
		my ($median_unambig_freq,$median_high_conf_int) = get_median_unambig_freq($R, $nuc_aa_codon_tally, $type, $converter, $verbose);		# Get the median unambig_nonconsensus_frequency and the high end of the 95% confidence interval for this amplicon or dataset (for this type).  This can be used as a threshold of confidence that a variant is real
		foreach my $pos (sort {$a <=> $b} keys %{$nuc_aa_codon_tally->{$type}}){
			my $ref_sequence = $cds->{$type}->{$pos};
			unless ($ref_sequence){
				print STDERR "Error in print_reports: no reference $type: $gene $pos\n";	
				if ($debug){	print Dumper($cds->{$type});		}
				if ($debug){	print Dumper($nuc_aa_codon_tally->{$type}->{$pos});	}
			}
			my $cons_with_num = find_key_with_biggest_value($nuc_aa_codon_tally->{$type}->{$pos}, 2);		# Assuming there won't be a tie for most abundant.  
			my ($consensus,$numConsensus) = split(/, /, $cons_with_num);
			my $diff = "";
			$diff = "DIFF" if ($consensus ne $ref_sequence);
			my $name = $gene.":".$pos.":".$consensus;
			my $coverage = total(values(%{$nuc_aa_codon_tally->{$type}->{$pos}})); 

			# Also compute values for columns: unambigCoverageDepth	unambigConsensus numUnambigConsensus	numUnambigNonConsensus	numMajorAltAllele	majorAltAllele
			my ($unambig_coverage, $unambig_consensus, $num_unambig_consensus, $num_unambig_nonconsensus, $major_alt_allele, $num_major_alt_allele) = get_unambig_values($nuc_aa_codon_tally, $type, $pos, $converter, $verbose);
			my ($cMAF,$low,$high,$pass_threshold) = ("","","","no");
			if ($unambig_coverage > 0){
				$cMAF = sprintf("%.7f", $num_unambig_nonconsensus / $unambig_coverage);		# Assign cumulative Minor Allele Frequency (if there are some unambiguous residues)
				($low,$high) = get_confidence_interval($R, $num_unambig_nonconsensus, $unambig_coverage);
				$pass_threshold = "yes" if ($low > $median_high_conf_int);
			} 

			# Now print the counts, number of amino acids represented, ref aa counts, and nonref aa counts.
			if ($type eq 'nuc'){
				print $nucleotide_report_fh join "\t", $name, $gene, $pos, $ref_sequence, $consensus, $diff, $label, $coverage, $unambig_coverage, $unambig_consensus, $num_unambig_consensus, $num_unambig_nonconsensus, $major_alt_allele, $num_major_alt_allele, $cMAF, $low, $high, $median_unambig_freq, $median_high_conf_int, $pass_threshold;
				print_tally_line($nucleotide_report_fh, 	$nuc_aa_codon_tally->{$type}->{$pos}, $NUC);
			}
			elsif($type eq 'codon'){
				print $codon_report_fh 		join "\t", $name, $gene, $pos, $ref_sequence, $consensus, $diff, $label, $coverage, $unambig_coverage, $unambig_consensus, $num_unambig_consensus, $num_unambig_nonconsensus, $major_alt_allele, $num_major_alt_allele, $cMAF, $low, $high, $median_unambig_freq, $median_high_conf_int, $pass_threshold;
				print_tally_line($codon_report_fh, 		$nuc_aa_codon_tally->{$type}->{$pos}, \@all_codons);
			}
			else {
				print $amino_acid_report_fh join "\t", $name, $gene, $pos, $ref_sequence, $consensus, $diff, $label, $coverage, $unambig_coverage, $unambig_consensus, $num_unambig_consensus, $num_unambig_nonconsensus, $major_alt_allele, $num_major_alt_allele, $cMAF, $low, $high, $median_unambig_freq, $median_high_conf_int, $pass_threshold;
				print_tally_line($amino_acid_report_fh, 	$nuc_aa_codon_tally->{$type}->{$pos}, $AA);
			}
		}
	}

	close($nucleotide_report_fh);
	close($codon_report_fh);
	close($amino_acid_report_fh);
}
#-----------------------------------------------------------------------------
sub get_unambig_values {
	# This subroutine computes values for a single position: unambigCoverageDepth	unambigConsensus numUnambigConsensus	numUnambigNonConsensus	numMajorAltAllele	majorAltAllele
		# unambigCoverageDepth.  Compute coverage depth for non-ambiguous residues.  For amino acid coverage, only consider regular amino acids (not *, X, or -).  For nuc, only consider A, C, T, G (not N).  For codon, only consider codons with A, C, T, and G (not codons with N).
		# unambigConsensus.  If Consensus residue is an ambiguous residue, report the next most abundant non-ambiguous residue, otherwise, report the consensus.  
		# numUnambigConsensus.  Report the tally number for the consensus allele.  Report '0' if the consensus is an ambiguous residue.
		# numUnambigNonConsensus.  Sum up all residues other than the consensus allele, not including ambiguous residues.  Should be equal to unambigCoverageDepth - numUnambigConsensus.  
		# majorAltAllele.  This is the major alternate allele (i.e., second highest count after the consensus residue, not including ambiguous characters)
		# numMajorAltAllele.  This is the tally number for the majorAltAllele. 
	my ($nuc_aa_codon_tally, $type, $pos, $converter, $verbose) = @_;

	my $unambig_residues = get_unambig_tally($nuc_aa_codon_tally->{$type}->{$pos}, $type, $converter);	 # Returns hashref with the ambiguous residues removed
	my $unambig_coverage = total(values(%$unambig_residues));
	my $unambig_consensus_with_num = find_key_with_biggest_value($unambig_residues, 2);
	my ($unambig_consensus,$num_unambig_consensus) = split(/, /, $unambig_consensus_with_num);
	my $num_unambig_nonconsensus = $unambig_coverage - $num_unambig_consensus;
	my $tally_without_consensus = remove_one_from_tally($unambig_residues,$unambig_consensus);
	my ($major_alt_allele,$num_major_alt_allele) = ("", 0);	# Default empty, in case there is no alternate allele.
	if (%$tally_without_consensus){
		# Then there are alternate allele(s)
		my $major_alt_alleles = find_key_with_biggest_value($tally_without_consensus, 2, 1);	# returns AoA of all max alleles (will usually be an AoA with just one entry, but will have multiple entries if there is a tie for the max number)
		my @alt_alleles;
		my @alt_counts;
		foreach (sort @$major_alt_alleles){
			my ($alt,$count) = @$_;
			#print STDERR "alt: $alt, $count\n";
			push @alt_alleles, $alt;
			push @alt_counts, $count;
		}
		$major_alt_allele = join ",", @alt_alleles;
		$num_major_alt_allele  = join ",", @alt_counts;
	}
	
	if ($verbose){		# $num_unambig_nonconsensus > 0 && 
		print STDERR "\n$pos, $type\n";
		print STDERR Dumper($nuc_aa_codon_tally->{$type}->{$pos});
		print STDERR "unambig\n";
		print STDERR Dumper($unambig_residues);
		print STDERR "unambig_cov: $unambig_coverage\n";
		print STDERR "tally without consensus\n";
		print STDERR Dumper($tally_without_consensus);
		print STDERR "major alt: $major_alt_allele, $num_major_alt_allele\n\n";
		#exit;
	}
	return ($unambig_coverage, $unambig_consensus, $num_unambig_consensus, $num_unambig_nonconsensus, $major_alt_allele, $num_major_alt_allele);
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
sub get_unambig_tally {
	my ($tally_hashref, $nuc_codon_aa, $converter) = @_;
	# Takes a hashref tally for nucleotides, codons, or amino acids (for a particular position) and removes all but the standard nucleotides, codons, and amino acids
	# For nucleotides, only report A, C, T, G
	# For codons, only report codons that unambiguously resolve to a known amino acid (using the %converter hash, which takes into account wobble nucleotides)	
		# Or, should we only report codons that are made up of A, C, T, and G?
	# For amino acids, only report qw(A C D E F G H I K L M N P Q R S T V W Y)
	my %converter = %$converter;

	my %aa; 
	foreach ( qw(A C D E F G H I K L M N P Q R S T V W Y *) ){
		$aa{$_} = 1;
	}
	my %nuc;
	foreach ( qw(A C T G) ){
		$nuc{$_} = 1;
	}

	my $unambig_tally_hashref;

	foreach my $res (keys %$tally_hashref){
		my $unambig = 0;	# Becomes 1 if found to be an unambiguous residue
		if ($nuc_codon_aa eq 'nuc'){
			$unambig++ if (exists($nuc{$res}));
		}
		elsif($nuc_codon_aa eq 'aa'){
			$unambig++ if (exists($aa{$res}));
		}
		elsif($nuc_codon_aa eq 'codon'){
			my $codon_to_aa = "";
			if (exists($converter{$res})){
				$codon_to_aa = $converter{$res};
			}
			$unambig++ if ($codon_to_aa && exists($aa{$codon_to_aa}));		
		}
		$unambig_tally_hashref->{$res} = $tally_hashref->{$res} if $unambig;
	}

	return $unambig_tally_hashref;
}
#-----------------------------------------------------------------------------
sub remove_one_from_tally {
	# Takes a hashref tally for nucleotides, codons, or amino acids (for a particular position) and removes the most abundant residue
	my ($tally_hashref,$key_to_remove) = @_;
	my %tally_without_consensus = %$tally_hashref;  # Make a copy to modify

	delete($tally_without_consensus{$key_to_remove});

	return \%tally_without_consensus;
}
#-----------------------------------------------------------------------------
sub print_merged_report {
	my ($nuc_aa_codon_tally, $label, $gene, $CODON, $cds, $prefix) = @_;
	# Print Reports to filehandles $nucleotide_report_fh, $codon_report_fh, and $amino_acid_report_fh
	#name(gene:aminoAcidPosition:consensusAminoAcid)	gene nucleotidePosition refNucleotide consensusNucleotide aminoAcidPosition refCodon	consensusCodon	refAminoAcid	consensusAminoAcid	Sample	coverageDepth	topNucleotide	secondNucleotide	thirdNucleotide	otherNucleotide	countTopNucleotide	countSecondNucleotide	countThirdNucleotide	countOtherNucleotide	topCodon	secondCodon	thirdCodon	otherCodon	countTopCodon	countSecondCodon	countThirdCodon	countOtherCodon	topAminoAcid	secondAminoAcid	thirdAminoAcid	otherAminoAcid	countTopAminoAcid	countSecondAminoAcid	countThirdAminoAcid	countOtherAminoAcid
	# $nuc_aa_codon_tally format example: $nuc_aa_codon_tally->{'aa'}->{$pos}->{$aa}++;
	
	my $merged_report 		= $prefix . ".merged.tally.xls";

	my $merged_report_fh		= open_to_write("$merged_report");

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



	my @all_codons = @$CODON;	# Make a copy to modify
	for (@all_codons){
		s/:.+//;		# Was s/:\w+//;  changed to s/:.+//; for numTAA:* and others with *
	}
		
	
	my $first_nuc_pos = (sort {$a <=> $b} keys %{ $nuc_aa_codon_tally->{'nuc'} })[0]; 
	#						print "first nuc position: $first_nuc_pos\n"; exit; 
	my $first_aa_pos =  (sort {$a <=> $b} keys %{ $nuc_aa_codon_tally->{'aa'} })[0]; 
	#						print "first aa position: $first_aa_pos\n"; exit;
	my $last_aa_codon_line = "\t"x16 ."\n";			# For repeating the amino acid and codon information for each nucleotide base of the codon.  Default will be an empty line (keep tabs for proper spreadsheet spacing) in case there are a few nt before the first codon.  
	foreach my $nuc_pos (sort {$a <=> $b} keys %{$nuc_aa_codon_tally->{'nuc'}}) {		
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
		if ($ref_nuc && $ref_aa && $ref_codon){
			print $merged_report_fh join "\t", $name, $gene, $nuc_pos, $position_in_codon, $ref_nuc, $consensus_nuc, $diff, $aa_pos, $ref_codon, $consensus_codon, $ref_aa, $consensus_aa, $label, $coverage;
		}
		else {
			print STDERR "Reference missing for position $nuc_pos, skipping... nuc_pos: $nuc_pos, aa_pos: $aa_pos, cons_nuc: $consensus_nuc, cons_codon: $consensus_codon, cons_aa: $consensus_aa\n";
			next;
		}

		# Print Nucleotide top hits, then Codon, then Amino Acid, e.g., topNucleotide	secondNucleotide	thirdNucleotide	otherNucleotide	countTopNucleotide	countSecondNucleotide	countThirdNucleotide	countOtherNucleotide	
		my $arrayRef = [ ['nuc', $nuc_pos], ['codon', $aa_pos], ['aa', $aa_pos] ];
		foreach my $pairs (@$arrayRef){
			my ($type, $pos) = @$pairs;
	#								print "$type $pos\n";
			if (exists($nuc_aa_codon_tally->{$type}->{$pos})){
				print_top_hits($nuc_aa_codon_tally->{$type}->{$pos}, $merged_report_fh);
			}
			else {
				print $merged_report_fh "\t"x8; 	# print blank fields.  Shouldn't have any because the nucleotide sequences are trimmed to the first codon before tallying up nucleotides and amino acids.  Although maybe not trimmed at the 3' end of reads...
			}
		}
		
		# Now print the counts, number of amino acids represented, ref aa counts, and nonref aa counts. 
		
		print $merged_report_fh "\n"; 
	}
	close($merged_report_fh);
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
	my ($tally, $merged_report_fh) = @_;
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
	my ($tally, $merged_report_fh) = @_;
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
	my ($nuc_aa_codon_tally, $label, $gene, $variant_threshold, $prefix, $start_time, $verbose) = @_;
	# Filters the variants (nuc, aa, codon) for those that are above a certain threshold in frequency 
	
	# First get the variants that are above a frequency threshold
		#Type\tgene\tposition\tconsensus\tvariant\tcount\tcoverage\tfrequency\n";
	 	# type is nucleotide, codon or aa.  
	 	# position is either aa or nuc position
	 	# consensus is the one that is highest
	 	# variant is the one with frequency above threshold
	 	# Count is the count for this particular variant
	 	# Coverage is the number of reads at this position
	 	# Frequency is count/coverage
	 # Compute interval of 95% confidence for variant frequency

	my $variants			= $prefix . ".variants.minfreq".$variant_threshold.".xls";		# Maybe add threshold to the name.  Or in the header of the file.

	my $variants_fh 		= open_to_write("$variants");


	print $variants_fh "## Variants above threshold frequency $variant_threshold\n";
	print $variants_fh "#Type\tgene\tposition\tconsensus\tvariant\tcount\tcoverage\tfrequency\tfilter\t95%_CI_low\t95%_CI_high\n";

	print STDERR "Printing out variants report. ";
	&elapsed($start_time, 'Elapsed', $verbose);

	my $R = Statistics::R->new();		# Open a bridge to R, to pass to get_confidence_interval sub

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
					my $frequency = sprintf("%.7f",  $nuc_aa_codon_tally->{$type}->{$pos}->{$variant} / $coverage);
					my ($min,$max) = get_confidence_interval($R, $count,$coverage);
					print $variants_fh join "\t", $type, $gene, $pos, $consensus, $variant, $count, $coverage, $frequency, $variant_filter, $min, $max;
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

	close($variants_fh);
}
#-------------------------------------------------------------------------------
sub get_confidence_interval {
	# This subroutine takes the $count and $coverage at a position, sends to R to run binom.test and returns the 95% confidence interval
	# get_conf_int <- function(x, n){  t <- binom.test(as.integer(x), as.integer(n)); return(max(t$conf.int)); }

	my ($R, $count, $coverage) = @_;

#	my $R = Statistics::R->new();		# http://search.cpan.org/~fangly/Statistics-R-0.31/lib/Statistics/R.pm
	$R->set('x', $count);
	$R->set('n', $coverage);
	$R->run(q`t = binom.test(as.integer(x), as.integer(n))`);
	$R->run(q`min = min(t$conf.int)`);
	$R->run(q`max = max(t$conf.int)`);
	my $min = sprintf("%.7f", $R->get('min') );
	my $max = sprintf("%.7f", $R->get('max') );

	return ($min,$max);
}
#-------------------------------------------------------------------------------
sub get_median_unambig_freq {
	# This subroutine takes an individual "type" tally, e.g., aa, nuc, codon (first level in from $nuc_aa_codon_tally), 
	# computes the unambig_coverage and unambigNonConsensus for each position of that type,
	# determines the 95% confidence interval (CI) for each position,
	# then computes the median frequency and high end of the 95% CI

	my ($R, $nuc_aa_codon_tally, $type, $converter, $verbose) = @_;

	my @unambig_freq;		# Array to store all unambiguous frequencies.  This includes all ACTG nonconsensus bases (cumulative Minor Allele Frequency or cMAF)
	my @unambig_freq_high_conf_int;		# Array to store high end of 95% CI for all frequencies

	foreach my $pos (sort {$a <=> $b} keys %{$nuc_aa_codon_tally->{$type}}){
		my ($unambig_coverage, $unambig_consensus, $num_unambig_consensus, $num_unambig_nonconsensus, $major_alt_allele, $num_major_alt_allele) = get_unambig_values($nuc_aa_codon_tally, $type, $pos, $converter, $verbose);
		if ($unambig_coverage > 0){
			push @unambig_freq, sprintf("%.7f", $num_unambig_nonconsensus / $unambig_coverage );
			my ($low,$high) = get_confidence_interval($R, $num_unambig_nonconsensus, $unambig_coverage);
			push @unambig_freq_high_conf_int, $high;
		}
	}
	my $median_unambig_freq = median(@unambig_freq);
	my $median_high_conf_int = median(@unambig_freq_high_conf_int);
	return ($median_unambig_freq,$median_high_conf_int);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
1;