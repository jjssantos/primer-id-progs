#-------------------------------------------------------------------------------
#----                                MODULE NAME                            ----
#-------------------------------------------------------------------------------
package gff_utils;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use IO::File;
use IO::Zlib; 		#qw(:gzip_external 1);
use Carp;
use Bio::SeqIO;
use Data::Dumper;
use Cwd;
use File::Basename;

# Custom modules
use Fasta_utils;

our @ISA = qw(Exporter);
our @EXPORT = qw(
	parse_gff
	get_gff_attribute
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
# this module since they are deprecated in that script but they might be useful 
# for another script.


#-------------------------------------------------------------------------------
#----------------------------------- FUNCTIONS ---------------------------------
#-------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
sub parse_gff {
	# Takes a gff file and a genome multi-fasta file of chromosomes and outputs a hashref containing the CDS starts, stops, strand, nucleotide sequence and protein sequence for each gene.  Also, nucleotide, amino acid and codons tables for lookup by position for each gene.  (Is it possible to have multiple CDS per gene in a gff file?  How would that affect this?)
	# See primer-id-progs v0.1.1 for deprecated usage.
	my ($gff, $ref_fasta) = @_;
	my $genes;	# hashref to store the coordinates for the CDS of genes on the chromosomes.  first keys, chr; second keys, gene id; third keys start, stop, starts, stops, strand, nucleotide, protein, nuc->pos, aa->pos, codon->pos.
	my $readfh = open_to_read($gff);
	
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
	my $in  = Bio::SeqIO->new(-file => "$ref_fasta" ,
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
				my $seqstr_nuc = Fasta_utils::rev_comp($seqstr_nuc);	# Not tested.  
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
#--------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
1;