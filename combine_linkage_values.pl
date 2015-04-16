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
use File::Temp;
use File::Spec;
#use Statistics::R;
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

my $save = Cwd::cwd();
my $files;
my $verbose;
my $output;
my $gzip;
my $var = "aa,codon,nuc";
my $mu = 0.5;
my $gamma = 0.01;
my $FDR = 0.05;
my $quasi;
my $keeptmp;
my $method = 'fisher';
my $prefix = "Linkage_plot";
GetOptions('save|s=s' => \$save, 'output=s' => \$output, 'verbose' => \$verbose, 'files=s' => \$files, 'gzip' => \$gzip, 'mu=s' => \$mu, 'gamma=s' => \$gamma, 'FDR|F=s' => \$FDR, 'quasi|q' => \$quasi, 'var|v=s' => \$var, 'keeptmp' => \$keeptmp, 'method|m=s' => \$method, 'prefix=s' => \$prefix);

my $Rscript_loc = '/usr/local/bio_apps/R-3.1.0/bin/Rscript';
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $exe = basename($0);

my $usage = "$exe <variants_linkage_file_1> [ <variants_linkage_file_2> ... ]
variants_linkage files are the output of calculate_linkage_disequilibrium.  Each input file
is a separate replicate.  Single replicate is fine too.  
$exe takes variants linkage file from calculate_linkage_disequilibrium.pl and plots them in R.
Requires R on the path and requires the 'network' library to be installed in R.

Example input variants file:
#group	sample	type	gene	comparison	pos1	c1	v1	v1freq	pos2	c2	v2	v2freq	c1c2	c1v2	v1c2	v1v2	p-value	OR	FDR
Group1	Sample1	nuc	Seq12_093009_HA_cds	437:A:C:713:G:A	437	A	C	0.04494	713	G	A	0.00924	25448	255	1255	0	8.18E-06	0	1.17E-04
Group1	Sample1	nuc	Seq12_093009_HA_cds	437:A:C:471:A:C	437	A	C	0.04494	471	A	C	0.02352	25045	656	1255	0	2.83E-14	0	1.22E-12
Group1	Sample1	nuc	Seq12_093009_HA_cds	437:A:G:713:G:A	437	A	G	0.03477	713	G	A	0.00924	25448	255	970	0	0.000139402	0	1.50E-03
Group1	Sample1	nuc	Seq12_093009_HA_cds	437:A:G:471:A:C	437	A	G	0.03477	471	A	C	0.02352	25045	656	970	0	4.40E-11	0	9.45E-10
Group1	Sample1	nuc	Seq12_093009_HA_cds	471:A:C:713:G:A	471	A	C	0.02352	713	G	A	0.00924	27015	255	656	0	0.005284684	0	3.25E-02
Group1	Sample1	codon	Seq12_093009_HA_cds	146:AAC:AGC:157:AAA:AAC	146	AAC	AGC	0.0347	157	AAA	AAC	0.02349	24993	655	969	0	4.37E-11	0	8.73E-10
Group1	Sample1	codon	Seq12_093009_HA_cds	146:AAC:AGC:238:GGT:GAT	146	AAC	AGC	0.0347	238	GGT	GAT	0.00913	25429	255	969	0	0.00013933	0	1.39E-03
Group1	Sample1	codon	Seq12_093009_HA_cds	146:AAC:AGC:157:AAA:GAA	146	AAC	AGC	0.0347	157	AAA	GAA	0.00122	24993	32	969	0	0.6324233	0	1.00E+00
Group1	Sample1	codon	Seq12_093009_HA_cds	146:AAC:ACC:238:GGT:GAT	146	AAC	ACC	0.0449	238	GGT	GAT	0.00913	25429	255	1254	0	8.18E-06	0	1.09E-04
Group1	Sample1	codon	Seq12_093009_HA_cds	146:AAC:ACC:157:AAA:GAA	146	AAC	ACC	0.0449	157	AAA	GAA	0.00122	24993	32	1252	2	0.6767571	1.24764	1.00E+00
Group1	Sample1	codon	Seq12_093009_HA_cds	146:AAC:ACC:157:AAA:AAC	146	AAC	ACC	0.0449	157	AAA	AAC	0.02349	24993	655	1252	0	2.83E-14	0	1.13E-12
Group1	Sample1	codon	Seq12_093009_HA_cds	157:AAA:AAC:238:GGT:GAT	157	AAA	AAC	0.02349	238	GGT	GAT	0.00913	26970	255	655	0	0.005285926	0	3.02E-02
Group1	Sample1	codon	Seq12_093009_HA_cds	157:AAA:GAA:238:GGT:GAT	157	AAA	GAA	0.00122	238	GGT	GAT	0.00913	26970	255	34	0	1	0	1.00E+00
Group1	Sample1	aa	Seq12_093009_HA_cds	146:N:S:157:K:N	146	N	S	0.0347	157	K	N	0.02352	24997	656	969	0	4.43E-11	0	7.31E-10
Group1	Sample1	aa	Seq12_093009_HA_cds	146:N:S:157:K:E	146	N	S	0.0347	157	K	E	0.00122	24997	32	969	0	0.632412	0	1.00E+00
Group1	Sample1	aa	Seq12_093009_HA_cds	146:N:S:238:G:D	146	N	S	0.0347	238	G	D	0.00913	25433	255	969	0	0.000139298	0	1.15E-03
Group1	Sample1	aa	Seq12_093009_HA_cds	146:N:T:157:K:E	146	N	T	0.04494	157	K	E	0.00122	24997	32	1253	2	0.6768856	1.246844	1.00E+00
Group1	Sample1	aa	Seq12_093009_HA_cds	146:N:T:157:K:N	146	N	T	0.04494	157	K	N	0.02352	24997	656	1253	0	2.83E-14	0	9.33E-13
Group1	Sample1	aa	Seq12_093009_HA_cds	146:N:T:238:G:D	146	N	T	0.04494	238	G	D	0.00913	25433	255	1255	0	8.18E-06	0	9.00E-05
Group1	Sample1	aa	Seq12_093009_HA_cds	157:K:N:238:G:D	157	K	N	0.02352	238	G	D	0.00913	26971	255	656	0	0.005305044	0	1.95E-02
Group1	Sample1	aa	Seq12_093009_HA_cds	157:K:E:238:G:D	157	K	E	0.00122	238	G	D	0.00913	26971	255	34	0	1	0	1.00E+00

OPTIONS:
	-F/--FDR	FDR threshold for inclusion of variants in the analysis.  Default = 0.05.
	-v/--var	Type of variant(s) to consider, comma-delimited list.  Possible values: 
			aa codon nuc.  Default is \"aa,codon,nuc\" (i.e., all variants).
	-m/--method	Method to use for merging linkage values if multiple files are provided.
			Possible values are 'fisher' and 'counts'.  'fisher' will use Fisher's method
			to combine the p-values.  'counts' will combine the counts and recompute a 
			p-value and fdr.  Default is 'fisher'.  
	--prefix	Prefix for the output files.  Default = 'Linkage_plot'.
	--keeptmp	Keep temporary files
	-q/--quasi	Calculate quasi-cliques in network using mgqce.  Default is just to make
			a network plot.
	--mu		Mu value for mgqce.  Default 0.5.
	--gamma		Gamma value for mgqce.  Default 0.01.
	-s/--save	Directory in which to save files. Default = pwd.  If folder doesn't exist, 
			it will be created.
";

# To do:
# Maybe add filter for only nonsyn codons?
# Work on optimizing the circle sizes and edge thickness...  see /Users/oleraj/Dropbox/NIAID_Data/WInce/NIH_Research_Festival_Poster/Some_analysis_for_figures/linkage for examples of where it is needed!
# Clean up a little (subroutine names, etc.)
# Merge the input files and output a merged output file!

# Change log
# 2014-04-17
# Created script, starting from calculate_quasi_cliques.pl
# 2014-09-09
# Added merge_files sub.
# Modified $edges variable in calculate_quasi_cliques sub to be a HoA instead of an AoA.  
# Modified expected input format -- added columns for v1freq and v2freq.
# Added (or made functional) options --FDR, --var, --method, --keeptmp --quasi
# Add filter for aa, nuc, codon (default to do all, but can select only aa, or perhaps nonsyn codons or nuc?)
# Added R graphs
# Added File::Spec->rel2abs for output directory
# 2015-03-23
# Modified expected input file format to match output of calculate_linkage_disequilibrium.pl
# 2015-04-16
# Fixed a bug where $save is passed to File::Spec even if it doesn't exist.  Made default for $save to be Cwd::cwd().

unless($ARGV[0]){
	print STDERR $usage;
	exit;
}

my $save_dir = File::Spec->rel2abs( $save) || Cwd::cwd();
unless (-d $save_dir){	mkdir($save_dir) or warn "$!\n"	}

my $start_time = time;
my @suffixes = (qw(.bed .bed.gz .bed12 .bed12.gz .txt .txt.gz .BED .BED.gz .BED12 .BED12.gz .fasta .fa .FA .FASTA .FAS .fas), @SUFFIXES);	#for fileparse.  Feel free to add more accepted extensions.  @SUFFIXES comes from aomisc.pm.  

# Make a hash of the types of variants allowed (nuc, aa, and/or codon)
my @types = split(/,/, $var);
my $types; #Hashref
foreach my $type (@types){
	$types->{$type}++;
}

# Set up global temporary directory
my $tempdir;
my $template = "comb_link_value_tempXXXXX";
if ($keeptmp){
	$tempdir = File::Temp->newdir( $template, DIR => $save_dir, CLEANUP => 0 );		# CLEANUP => 0 so it will not be deleted 
}
else {
	$tempdir = File::Temp->newdir( 	$template  );		# CLEANUP => 1 by default so it will be deleted when out of scope
}

# Merge files if multiple files are provided
my $merged_file;
if (@ARGV > 1){
	$merged_file = merge_files(\@ARGV, $method);	# Merge linkage values based on designated method. 
}
else {
	$merged_file = $ARGV[0];
}


foreach my $file (@ARGV){
#	my ($filename,$dir,$ext) = fileparse($file,@suffixes);		# fileparse($file,qr/\.[^.]*/);
#	my $newfilename = $save_dir.'/'.$filename.'__'.$variable.'.bed';

	my $stats = read_variants($file);
	calculate_quasi_cliques($stats);
	
}




&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub merge_files {
	my ($files,$method) = @_;
	my @files = @$files;
	my $merged_filename;		# Name for merged file.  
	
	
	
	return $merged_filename;
}
#-----------------------------------------------------------------------------
sub read_variants {
	my $file = shift;
	my $stats;	# AoA containing all lines of the variants file
	
	my $readfh = open_to_read($file);
	while (<$readfh>){
		chomp;
		my @line = split(/\t/);
		next if (/^#/);
		push @$stats, \@line;
	}
	
	close($readfh);
	return $stats;
}
#-----------------------------------------------------------------------------
sub calculate_quasi_cliques {
	my $stats = shift;
	# $stats is an AoA where each line represents a comparison of two variants
    #      [
    #        'aa',
    #        'Seq12_093009_HA_cds',
    #        '136',
    #        'K',
    #        'E',
    #        '238',
    #        'G',
    #        'D',
    #        8569,
    #        199,
    #        101,
    #        0,
    #        '0.176207',
    #        '0',
    #        '2.265519e-01'
    #      ],
    
    # New format output of calculate_linkage_disequilibrium.pl :
	#	#type	gene	pos1	c1	v1	v1freq	pos2	c2	v2	v2freq	c1c2	c1v2	v1c2	v1v2	p-value	OR	FDR
	#	nuc	Seq12_093009_HA_cds	406	A	G	0.00684	469	A	G	0.00125	27701	35	191	0	1	0	1.000000e+00
	#	nuc	Seq12_093009_HA_cds	406	A	G	0.00684	559	A	G	0.00269	27662	75	191	0	1	0	1.000000e+00
    # Yet newer format out of calculate_linkage_disequilibrium.pl (as of 2014-09-06) :
	#	#group	sample	type	gene	comparison	pos1	c1	v1	v1freq	pos2	c2	v2	v2freq	c1c2	c1v2	v1c2	v1v2	p-value	OR	FDR
	#	Group1	Sample1	nuc	Seq12_093009_HA_cds	437:A:C:713:G:A	437	A	C	0.04494	713	G	A	0.00924	25448	255	1255	0	8.18E-06	0	1.17E-04
	#	Group1	Sample1	nuc	Seq12_093009_HA_cds	437:A:C:471:A:C	437	A	C	0.04494	471	A	C	0.02352	25045	656	1255	0	2.83E-14	0	1.22E-12
	#	Group1	Sample1	nuc	Seq12_093009_HA_cds	437:A:G:713:G:A	437	A	G	0.03477	713	G	A	0.00924	25448	255	970	0	0.000139402	0	1.50E-03



    # I want to plot the network in R using 'network' package.
    # To do so, I need a graph file of edges and node attributes file
    # It appears that the attributes file row names (numbers starting at 1) are what is to be used to identify the nodes in the edges file.  
    # Read the $stats struct
    	# Get a unique list of the alleles and classify as "var" or "cons".  Store as HoH $alleles, where first key is the allele id and second key is var_cons.   
    	# Store edges in $edges AoA
    # Read the $alleles HoH and assign $alleles->{$allele_id}->{'node_num'}
    # Read the $edges AoA and print out the node_num for each pair (from $alleles); two pairs per row.  
    # After much trial and error, I found that network() command for vertex.attr requires only numeric values for attributes.  vertex.attrnames can be strings though.  The edges can have strings.  
	# I still want to print out the names and attributes file separately from the edges because that is the only way to be sure that the nodes are given the correct order within the network object.  
	
	# Files for plotting in R
	my $temp_attrib	= $tempdir."/attributes.txt";
	my $temp_edges	= $tempdir."/edges.txt";
	my $temp_names	= $tempdir."/names.txt";
	my $temp_edge_lb= $tempdir."/edge_label.txt";
	my $attr_fh 	= open_to_write($temp_attrib);
	my $edges_fh 	= open_to_write($temp_edges);
	my $names_fh	= open_to_write($temp_names);
	my $edge_lb_fh	= open_to_write($temp_edge_lb);
	
	
	my $alleles;		# HoH to contain a non-redundant list of variants and consensus
	my $edges;			# HoA.  This needs to be a hash, since I see that some edges are found multiple times.  Make the pair of residues the key.  Store an array of FDR values for each pair of residues and then use the best for the plot (?). 
						# When does this occur that an edge will be found multiple times?  When there are multiple variant alleles and the variant alleles are linked with OR status is the same for both (>1 or <1) so that the consensus alleles are linked with each other.  It doesn't cause a duplication of the variant pairing, but it does duplicate the consensus pairing. e.g., 'codon-146-AAC:codon-238-GAT'  146 consensus codon is AAC, variant alleles AGC, ACC.  238 consensus codon GGT and variant allele GAT.  146 is compared to 238 twice (once for each allele in 146).  The OR = 0 in both, so codon-146-AAC is linked to codon-238-GAT in both.  p-values are '1.393300e-03', '1.090292e-04'.  Which to use?  They are kind of dependent p-values, so probably shouldn't combine them with Fisher's method.  Most of the time they are very similar (within 1 order of magnitude).  I'll just plot the stronger of the two.  
						
    foreach (@$stats){
 	if (scalar(@{$_}) < 20){
		print STDERR "Check input file to see if it has the correct number of columns.\n";
		exit 1;
	} 
	elsif (scalar(@{$_}) > 1){
    		my ( $group, $sample, $type, $gene, $comparison, $first_pos, $first_cons, $first_var, $first_var_freq, $second_pos, $second_cons, $second_var, $second_var_freq, $AB, $Ab, $aB, $ab, $p, $OR, $fdr ) = @{$_};
    		if ($fdr <= $FDR && exists($types->{$type})){		# Make sure the comparison passes the FDR threshold and is of the correct type.
				$first_cons 	= join "-", $type, $first_pos, $first_cons;		# Consensus residue at the first position in the comparison
				$first_var 	= join "-", $type, $first_pos, $first_var;		# Alternate residue at the first position
				$second_cons	= join "-", $type, $second_pos, $second_cons;	# Consensus residue at the second position
				$second_var	= join "-", $type, $second_pos, $second_var;	# Alternate residue at the second position
				
				my $total = $AB + $Ab + $aB + $ab; 
				
				my %first_cons = (	'var_cons' => "cons", 	'allele' => $first_cons,	'freq' => ($AB + $Ab)/$total );
				$alleles->{$first_cons} = \%first_cons;
				
				my %second_cons = (	'var_cons' => "cons", 	'allele' => $second_cons,	'freq' => ($AB + $aB)/$total );
				$alleles->{$second_cons} = \%second_cons;   

				my %first_var = ( 	'var_cons' => "var", 	'allele' => $first_var,		'freq' => ($aB + $ab)/$total );
				$alleles->{$first_var} = \%first_var;
				
				my %second_var = (	'var_cons' => "var", 	'allele' => $second_var,	'freq' => ($Ab + $ab)/$total );
				$alleles->{$second_var} = \%second_var;
				

				# Use OR to connect individual variants to consensus instead of just position. 
    			# If OR > 1, then variant is linked to variant and consensus linked to consensus; if OR < 1, then variant of one is linked to consensus of another.  Two lines printed for the graph per line in the original file. 
				if ($OR >= 1){
					# variant is linked to variant and consensus linked to consensus
					
					push @{$edges->{$first_var  . ":" . $second_var }}, $fdr;
					push @{$edges->{$first_cons . ":" . $second_cons}}, $fdr;
				}
				else {
					# variant of one is linked to consensus of another
					push @{$edges->{$first_var  . ":" . $second_cons}}, $fdr;
					push @{$edges->{$first_cons . ":" . $second_var }}, $fdr;
				}		
			}
    	}
    }
#	print Dumper($alleles); exit;

#    print Dumper($edges); exit;
    
    # Now sort the alleles (position in the @sorted_alleles AoH will be the node_num)
    my @sorted_alleles;  
    foreach my $allele (sort keys %$alleles){
    	push @sorted_alleles, $allele;
    }    
#   print Dumper(\@sorted_alleles);

	# Print out attributes file and names file for plotting in R (e.g., var or cons)
	# Also add 'node_num' to $alleles
   	print $attr_fh "consensus4.variant2\tfreq\n";
   	print $names_fh "allele\n";
   	my $i = 1;
	foreach my $allele (@sorted_alleles){
		print $names_fh "$allele\n";
		if ($alleles->{$allele}->{'var_cons'} eq 'cons'){
			print $attr_fh "4";		# cons = 4 (blue in plot.network)
		}
		else {
			print $attr_fh "2";		# var = 2 (red in plot.network)
		}
		print $attr_fh "\t$alleles->{$allele}->{'freq'}\n";
		
		$alleles->{$allele}->{'node_num'} = $i; 
		$i++; 
	}
	close ($attr_fh);

#	print Dumper($alleles); exit;
    
    # Read $edges and print out edges print out 'node_num' representing each node instead of the node name.     
	foreach my $edge (keys %$edges){
		my ($first, $second) = split(/:/, $edge);
		my @fdr_array = sort {$a <=> $b} @{$edges->{$edge}};
		my $fdr = $fdr_array[0];	# Use the smaller of the values if there is more than one.  
		print $edges_fh "$alleles->{$first}->{'node_num'}\t$alleles->{$second}->{'node_num'}\n";
		print $edge_lb_fh "$fdr\n";
	}
	close($edge_lb_fh);
	close($edges_fh);

	# Example graphing in R
	# names <- read.table("~/temp/names.txt", header=T)
	# edges <- read.table("~/temp/edges.txt")
	# attr <- read.table("~/temp/attributes.txt", header=T)
	# edglab <- read.table("~/temp/edge_label.txt")
	# library(network)
	# nw <- network(edges, vertex.attr=attr, directed=F)
	# plot.network(nw, displaylabels=T, vertex.col=attr$consensus4.variant2, label=names$allele, label.cex=0.5, vertex.cex=attr$freq+0.5, edge.label=format(edglab$V1, scientific=T, digits=2), edge.label.cex=0.5, edge.label.col=6, edge.lwd=-0.5*log(edglab$V1))
		# This plots red = variant, blue = consensus; thickness of edge represents the fdr value (thicker = lower fdr value); the size of the ball represents the frequency of the allele (bigger = higher frequency)
	
	
	# Or, write R script and run it with Rscript from the system.  
	my $temp_R_script = $tempdir . "/temp_script.R";
	my $temp_R_fh = open_to_write($temp_R_script);
	print $temp_R_fh	"setwd(\"$tempdir\")\n";		# Set working directory to $tempdir to read the files
	print $temp_R_fh	"names <- read.table(\"names.txt\", header=T)\n";
	print $temp_R_fh	"edges <- read.table(\"edges.txt\")\n";
	print $temp_R_fh	"attr <- read.table(\"attributes.txt\", header=T)\n";
	print $temp_R_fh	"edglab <- read.table(\"edge_label.txt\")\n";
	print $temp_R_fh	"library(network)\n";
	print $temp_R_fh	"nw <- network(edges, vertex.attr=attr, directed=F)\n";
	print $temp_R_fh	"setwd(\"$save_dir\")\n";	# Change working directory now to $save_dir for saving files
	my $types_string = join "_", sort keys %$types;
	my $outfile = $prefix . "_$types_string" . "_network.pdf"; 
	print $temp_R_fh	"pdf(\"$outfile\", height=8, width=8)\n";
	print $temp_R_fh	"plot.network(nw, displaylabels=T, vertex.col=attr\$consensus4.variant2, label=names\$allele, label.cex=0.5, edge.label=format(edglab\$V1, scientific=T, digits=1), edge.label.cex=0.5, edge.label.col=6, edge.lwd=-0.5*log(edglab\$V1))\n";
	print $temp_R_fh	"dev.off()\n";
	# These next few steps are testing an alternative method of plotting, scaling frequency separately based on residue designation as either variant or consensus.
	print $temp_R_fh	'meanvar = mean(subset(attr, attr$consensus4.variant2 == 2)$freq)' . "\n";		# Get mean value for all variant frequencies in the table.
	print $temp_R_fh	'meancons = mean(subset(attr, attr$consensus4.variant2 == 4)$freq)'. "\n";		# Get mean value for all consensus frequencies in the table.
	print $temp_R_fh	'for (i in row.names(attr)){if (attr[i, "consensus4.variant2"] == 2){attr[i, "centfreq"] = attr[i, "freq"] / meanvar} else {attr[i, "centfreq"] = attr[i, "freq"] / meancons} }' . "\n";	# transform the frequency by dividing the frequency by the mean of either the consensus frequencies or variant frequencies, depending on what type you are.  Store this in "centfreq" column, for "centered frequency"
	$outfile = $prefix . "_$types_string" . "_network_scaled_freqs.pdf"; 
	print $temp_R_fh	"pdf(\"$outfile\", height=8, width=8)\n";
	print $temp_R_fh	'plot.network(nw, displaylabels=T, vertex.col=attr$consensus4.variant2, label=names$allele, label.cex=0.5, vertex.cex=attr$centfreq+0.5, edge.label=format(edglab$V1, scientific=T, digits=1), edge.label.cex=0.5, edge.label.col=6, edge.lwd=-0.5*log(edglab$V1))' . "\n";
	print $temp_R_fh	"dev.off()\n";

	close($temp_R_fh);
	
	# if (check_for_Rscript()){
	# 	# Then we didn't find Rscript on the path.
	# 	my $copy_R_script_file_dir = $save_dir . "/$prefix" . "_network_plot_files";
	# 	print STDERR "Didn't find Rscript on the PATH!\nCopying the R script and files to $copy_R_script_file_dir.  You can run it manually like this:\ncd $copy_R_script_file_dir; Rscript temp_script.R\n";
	# 	system("cp -r $tempdir $copy_R_script_file_dir");
	# 	exit;
	# }
	# else {
		system("$Rscript_loc $temp_R_script");
	# }
	
	
	if ($quasi){
		# This is for making quasi-cliques with mgqce

		# Files for quasi-cliques in mgqce
		my $temp_graph	= $tempdir."/graph.mgqce.txt";
		my $temp_rgc 	= $tempdir."/graph.rgc.txt";
		my $temp_key	= $tempdir."/key.rgc.txt";
		my $temp_query	= $tempdir."/query.txt";
		my $temp_clique	= $tempdir."/cliques.txt";
		my $fh = open_to_write($temp_graph);    
	
		foreach (@$stats){
			if (scalar(@{$_}) > 1){
				my ( $group, $sample, $type, $gene, $comparison, $first_pos, $first_cons, $first_var, $second_pos, $second_cons, $second_var, $AB, $Ab, $aB, $ab, $p, $OR, $fdr ) = @{$_};
				if ($fdr < $FDR){
					my $first  = join "-", $type, $first_pos, $first_cons . "/" . $first_var; 
					my $second = join "-", $type, $second_pos, $second_cons . "/" . $second_var; 
					print $fh "$first\t$second\n";
				}
			}
		}
		close($fh);

		# Make mgqce/DENSE compatible graph and key files with rgc.
		my $cmd = "rgc $temp_graph $temp_rgc $temp_key";
						if ($verbose){	print STDERR "Running commdand: $cmd\n";	}
		system($cmd);
	
		# Make query file from $temp_key. 
		my $count = `cat $temp_key | wc -l`;
		my $query_fh = open_to_write($temp_query);
		my $key_fh = open_to_read($temp_key);
		print $query_fh "$count\n";
		while(<$key_fh>){
			chomp;
			my @line = split(/\s/);
			print $query_fh "$line[0]\n";
		}
		close($query_fh);
		close($key_fh);
	
		# Run mgqce to get quasi-cliques
		# E.g., mgqce graph.rgc.txt 0.75 0.01 query.txt out.txt
		# mgqce <graph> <mu> <gamma> <query> <cliques>
		# "mgqce enumerates all subgraphs of the graph where at least mu percent of subgraph is made up of "query" vertices and each vertices is adjacent to at least gamma percent of the other vertices."
		# "Note, gamma must be at least 0.5.  Also, the algorithm may enumerate the enriched cliques of a graph more "intelligently" (i.e., quickly) if you use a value of gamma = 0.999..., rather than gamma = 1."		
		my $mu = 0.5; 		# at least mu percent of subgraph is made up of "query" vertices .  My default 0.5
		my $gamma = 0.01;	# each vertex is adjacent to at least gamma percent of the other vertices.  My default 0.01.
		$cmd = "mgqce $temp_rgc $mu $gamma $temp_query $temp_clique > stdout.txt";
		system($cmd);
		`cat $temp_clique`;
	}
		
	# Now I need to convert back to variant information using the key
	# Graph in R.
    
	
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
