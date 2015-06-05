#!/usr/local/bio_apps/perl-5.16.2/bin/perl

use warnings;
$| = 1;

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


# Andrew J. Oler, PhD
# Computational Biology Section
# Bioinformatics and Computational Biosciences Branch (BCBB)
# OCICB/OSMO/OD/NIAID/NIH
# Bethesda, MD 20892
#
# andrew.oler@gmail.com
# 
# This package is free software; you can redistribute it and/or modify
# it under the terms of the GPL (either version 1, or at your option,
# any later version) or the Artistic License 2.0.


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
my $label;
my $debug;
my $keeptmp;
GetOptions('save=s' => \$save, 'output=s' => \$output, 'verbose' => \$verbose, 'files=s' => \$files, 'gzip' => \$gzip, 'label=s' => \$label, 'debug' => \$debug, 'keeptmp' => \$keeptmp, );

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $exe = basename($0);

my $usage = "$exe takes frequency tables (nucleotide, amino acid or codon) from 
convert_reads_to_amino_acid.pl and compares frequency distribution between two samples or 
two groups of samples using Fisher's Exact test with FDR correction.

Usage: $exe <Sample 1> <Sample 2>
First argument is a frequency table or a comma-delimited list of frequency tables from 
replicates for Sample 1.
Second argument is a frequency table or a comma-delimited list of frequency tables from
replicates for Sample 2.

Nucleotide input file has columns:
#name	gene	nucleotidePosition	refNucleotide	consensusNucleotide	refDiffConsensus	Sample	coverageDepth	numA	numC	numG	numT	numN	numOther
Basically, there should be 7 columns before the coverageDepth column.

Output file has columns
name	chr	coordinate	refBase	Sample	coverageDepth	numA	numC	numG	numT	numN	del	SumACGT	numNonzeroACGT	Ref	Nonref
One output file per sample
Note, it only prints rows if the coverage is at least 1% of the maximum coverage for a sample.
OPTIONS:
	--label		Comma-delimited list of labels for the two samples or groups of samples.  
			Defaults to the Sample column if there is a single replicate per sample, or 
			\"A,B\" if multiple replicates are used per sample.
	--keeptmp	Keep temporary files.  Default is to delete them (e.g., R script, individual pdf files)

	--nostats	Just make frequency table, don't run fisher's exact test.  Default is to do both.
	--printall	Print all rows, even if lower than 1% coverage.  Default is to not print any rows with 
			< 1% of the maximum coverage.
	--stranded	Bases from reads aligning to forward strand will be tallied separately from bases of
			reads aligning to the reverse strand.  Forces --nostats.  
	--save		Directory in which to save files. Default = pwd.  If folder doesn't exist, it 
			will be created.
";

# Change Log
# 2014-09-06
# Created script.  Borrowed code from MpileupFrequency.pl for statistical analysis 
# 2014-09-08
# Changed File::Temp::tempdir call to File::Temp::newdir to be compatible with File::Temp 0.23.  See http://www.dagolden.com/index.php/2109/why-the-latest-filetemp-might-surprise-you/
# 2015-06-05
# Fixed merge_samples to use the labels in case no replicates were provided.  The main issue is if both samples were run through convert_reads_to_amino_acid without supplying the --label option, both samples will have Sample_1 as Sample id, so no comparisons can be made in Fishers test. Fix that by substituting Sample_1 for label provided in --label.  

# To do
# When merging replicates, check for a change in the consensus based on the new frequencies and change the consensus accordingly.  
# Maybe allow different methods for merging replicates.  Default is to sum up the counts from each of the replicates.  Makes a sort of average; good for 2-3 replicates.  For 4 replicates, could possibly use log-linear model.  
# Set up parallelization for statistical analysis?  Maybe run via Statistics::R instead of Rscript?
# Make the merged file not temporary, saved in output.  Better name.  
# Add add_consensus_columns_to_frequency_tables.pl functionality to get these columns.  Also make sure all programs are compatible with the extra columns.  

unless ($ARGV[1]){	print STDERR "$usage\n";	exit;	}

my $save_dir = $save || Cwd::cwd();
unless (-d $save_dir){	mkdir($save_dir) or warn "$!\n"	}
#warn "save directory: $save_dir\n";

my $start_time = time;
my @suffixes = (qw(.bed .bed.gz .bed12 .bed12.gz .txt .txt.gz .BED .BED.gz .BED12 .BED12.gz .fasta .fa .FA .FASTA .FAS .fas), @SUFFIXES);	#for fileparse.  Feel free to add more accepted extensions.  @SUFFIXES comes from aomisc.pm.  


# Get labels
my @labels;
if ($label){
	@labels = split (/,/, $label);
}
else {
	my $lab = "A";
	foreach my $file ($ARGV[0], $ARGV[1]){
		if ($file =~ m/,/){
			push @labels, $lab;
		} 
		else {
			# Just a single replicate for this sample so get label from Sample column 
			($lab) = get_label($file);
			push @labels, $lab;
		}
		$lab = "B";	# second time through, $lab will be "B"
	}
}


# Get header and remove #
my @files = split(/,/, $ARGV[0]);
my @header = get_header($files[0]); 	# In case multiple files are provided
$header[0] =~ s/^#//; 


# Determine type of file (either nuc aa or codon)
my $type = get_residue_type($ARGV[0]);


# Set up global temporary directory
my $tempdir;
my $template = "comp_var_freq_tempXXXXX";
if ($keeptmp){
	# $tempdir = File::Temp->tempdir( $template, DIR => $save_dir );		# CLEANUP => 0 by default so it will not be deleted 		# Doesn't work in File::Temp 0.23
	$tempdir = File::Temp->newdir( $template, DIR => $save_dir, CLEANUP => 0 );		# CLEANUP => 0 so it will not be deleted 
}
else {
	$tempdir = File::Temp->newdir( 	$template, DIR => $save_dir  );		# CLEANUP => 1 by default so it will be deleted when out of scope
}


## Merge samples
# If there are replicates for either sample, merge them into one table of the same format as the input file
my $merged_sample1_file;
my $merged_sample2_file;

# Sample 1
if ($ARGV[0] =~ m/,/){
	# Replicates to merge for first Sample
	$merged_sample1_file = merge_replicates($ARGV[0], $labels[0]);
}
else {
	# No replicates to merge.
	$merged_sample1_file = $ARGV[0];
}

# Sample 2
if ($ARGV[1] =~ m/,/){
	# Replicates to merge for first Sample
	$merged_sample2_file = merge_replicates($ARGV[1], $labels[1]);
}
else {
	# No replicates to merge.
	$merged_sample2_file = $ARGV[1];
}


## Now take the merged replicate files and make one merged file for the statistical analysis in R.
my $merged_both_samples = merge_samples($merged_sample1_file,$merged_sample2_file, $labels[0], $labels[1]);
print STDERR "Merged file with both samples: $merged_both_samples\n";


## Now run the statistical analysis on the merged file
run_fisher_test($merged_both_samples);

#foreach my $file (@files){
#	read_and_write($file);
#}





&elapsed($start_time, 'Total', $verbose);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub get_label {
	# Takes a file, reads the header line and the line following the header and returns the
	# value for the first column that has a header with "Sample" as part of the name.
	# Also returns an arrayref of the headers.

	my $file = shift;
	my $label = "none";
	
	# Get header line and first data line from the file.
	my $readfh = open_to_read($file);
	my @header;
	my @first_line;
	LINE: while(<$readfh>){
		chomp;
		my @F = split(/\t/);
		if (m/^#/){
			@header = @F;
			next;
		}
		@first_line = @F;
		last LINE;
	}
	close($readfh);
	
	# Determine index for "Sample" column and then obtain name for label
	my $sample_index;
	COLUMN: for (my $i = 0; $i < @header; $i++){
		if ($header[$i] =~ m/Sample/i){
			$sample_index = $i;
			last COLUMN;
		}
	}
	$label = $first_line[$sample_index] if ($sample_index);
	return ($label,\@header);
}
#-----------------------------------------------------------------------------
sub get_residue_type {
	# Takes a file and determines whether it is a frequency table of nucleotides, amino
	# acids, or codons.  For now, it is simply looking for .nuc .aa or .codon in the file
	# name.  It depends on the fact that convert_reads_to_amino_acids.pl adds this to the 
	# file name. Perhaps in the future, I'll do something more sophisticated for this sub.
	# 
	my $files = shift;
	my @files = split(/,/, $files);	# in case multiple files are passed, take the first
	my $file = $files[0];
	my $type;		# Will be 'nuc', 'aa', or 'codon'
	foreach my $test (qw(nuc aa codon)){
		$type = $test if $file =~ m/$test.tally.xls/;
	}
	
	print STDERR "type: $type\n";	
	
	return $type;
}
#-----------------------------------------------------------------------------
sub merge_replicates {
	# Takes a comma-delimited list of frequency tables and merges the counts into a single
	# file of the same format 
	my ($files,$label) = @_;
	my @files = split(/,/, $files);
	my $hash;	# hashref to store the counts as I'm adding them up from the different files.  First keys are the "Position" column from the frequency table (column 3).  Second level keys are 'info' (string of first 7 columns), 'coverageDepth' and all 'num' columns in the file.
	my @header;			# All headers
	my @data_headers;	# order of headers for data columns, including coverageDepth and all "num" columns.
	
	foreach my $file (@files){
		my $readfh = open_to_read($file); 
		while(<$readfh>){
			chomp;
			# Nucleotide tally.xls file:
			# #name	gene	nucleotidePosition	refNucleotide	consensusNucleotide	refDiffConsensus	Sample	coverageDepth	numA	numC	numG	numT	numN	numOther
			# Seq12_093009_HA_cds:4:A	Seq12_093009_HA_cds	4	A	A		Sample_1	5495	5470	1	23	1	0	0
			
			# Amino acid tally.xls file:
			# #name	gene	aminoAcidPosition	refAminoAcid	consensusAminoAcid	refDiffConsensus	Sample	coverageDepth	numA	numC	numD	numE	numF	numG	numH	numI	numK	numL	numM	numN	numP	numQ	numR	numS	numT	numV	numW	numY	num*	numX	num-	numOther
			# Seq12_093009_HA_cds:2:K	Seq12_093009_HA_cds	2	K	K		Sample_1	5495	0	0	0	23	0	0	0	0	5469	0	0	0	0	1	1	0	0	0	0	0	1	0	0	0
			
			# Codon tally.xls file similar to amino acid tally.xls frequency file.  
			
			# Basically, the first few columns are the same for all types and these columns can be used as a key to the hash.  
			# The column "coverageDepth" and all columns starting with "num" should be saved for adding up between the different files.
			# I can think of one potential problem in adding these up where the consensus might not be the same among all replicates, so the key will not work.  Maybe I'll just use the reference base as the key.  I should calculate a new consensus if the tally says it has changed in the merged situation... for now, I'm not going to do that because the likelihood is very low.
			# For now, 1) Save the information for the columns before "coverageDepth" from the first file, and 2) Use the reference position as a key as it will be unchanged between samples, and will be good for sorting the data for printing.  
			
			my @F = split(/\t/); 	# Split by tabs, not whitespace because sometimes (most of the time) the refDiffConsensus column is empty
			if (m/^#/){
				@header = @F unless @header; 	# Just save from first file	
				@data_headers = @F[7..$#F] unless @data_headers;	# Just save from first file	
				next;
			}
			my $pos = $F[2]; 
			my $info_columns = join "\t", @F[0..5], $label; 	# All columns before Sample, plus the label determined for this sample.  Stored info about the position, reference residue, consensus residue, etc.  
			
			$hash->{$pos}->{info} = $info_columns unless (exists($hash->{$pos}->{info}));  	# Only store it if it isn't there yet, thus giving preference to first file in the list.
			
			for (my $i = 7; $i < @F; $i++){
				my $col = $header[$i]; 
				if ($F[$i] =~ m/^\d+$/){
					$hash->{$pos}->{$col} += $F[$i]; 	# Add the value to the current value in the hash. Only add if an integer (skip the entries in numOther such as this: "1(-:1)")
				}
				else {
					my $num = 0;
					if ($F[$i] =~ m/^(\d+)\D/){		# e.g., 5493(NNN:5464,TGN:4,TNN:25)
						$num = $1;
					}
					$hash->{$pos}->{$col} += $num;			# Added this so there is at least a value in that cell. 
				}
			}
		}
		close($readfh);
	}
	
	
#	print Dumper($hash); exit;

	# Prepare output file for merged 
	my $merged_file = $tempdir . "/$label".".mergedreps.tally.xls";
	my $writefh = open_to_write($merged_file);
	$header[0] =~ s/^#//;	# R doesn't like column headers to begin with #
	print $writefh join "\t", @header; 
	print $writefh "\n";
	
	foreach my $pos (sort {$a <=> $b} keys %$hash){
		my $info = $hash->{$pos}->{info};
		print $writefh $info;
		foreach my $col (@data_headers){
			if (exists($hash->{$pos}->{$col})){
				print $writefh "\t".$hash->{$pos}->{$col};
			}
			else {
				print STDERR "Doesn't exist in Hash: Label $label Position $pos Column: $col\n";
			}
		}
		print $writefh "\n";
	}
	 
	
	close($writefh);
	
	return $merged_file;
}
#-----------------------------------------------------------------------------
sub merge_samples {
	# Takes two frequency tables and merges them into one file of the format required
	# for the analysis in R.  Basically, we just need the name column and the residue
	# columns (e.g., A C G T).  We will also want any other useful information columns
	# we want to retain in the output.
	# For our current purposes, we can just concatenate the files, actually.
	# In the future, I'll probably want to compute a few extra columns, e.g., Ref, Nonref, SumResidues (like SumACGT), numNonzeroResidues (like numNonzeroACGT).
	# Or, we could compute these extra columns up above.
	my ($file1, $file2, $label1, $label2) = @_;
	
	my $labels_for_filename = join "_", $label1, $label2;
	my $merged_filename = $save_dir . "/Merged_" . $labels_for_filename . ".tally.xls";
	
	# Initialize output file and print out header
	my $merged_file_fh = open_to_write($merged_filename);
	my @header = get_header($file1);
	$header[0] =~ s/^#//;	# remove '#' character for R header
	print $merged_file_fh join "\t", @header;
	print $merged_file_fh "\n";
	
	my @files = ($file1, $file2);
	my @labels = ($label1, $label2);
	for (my $i = 0; $i < @files; $i++){
		my $header = column_header_lookup_hash($file1);
		my $sample_index = $header->{"Sample"};		# column index for Sample
		my $readfh = open_to_read($files[$i]);
		while(<$readfh>){
			my @F = split(/\t/);
			next if (m/coverageDepth/);	# skip header row
			unless ($F[$sample_index] && $F[$sample_index] eq $labels[$i]){
				# Sample label doesn't exist or doesn't match the input sample label, replace it with the input sample label
				$F[$sample_index] = $labels[$i];
			}
			print $merged_file_fh join "\t", @F;	
		}
		close ($readfh);
	}
	close ($merged_file_fh);
		
	return $merged_filename;
}
#-----------------------------------------------------------------------------
sub run_fisher_test {
	my $merged_file = shift;
	
	# Now write an R script to run the statistical analysis
	# Input file format (column headers):
	#	name	chr	coordinate	refBase	Sample	coverageDepth	numA	numC	numG	numT	numN	del	SumACGT	numNonzeroACGT	Ref	Nonref

	# name	gene	aminoAcidPosition	refAminoAcid	consensusAminoAcid	refDiffConsensus	Sample	coverageDepth	numA	numC	numD	numE	numF	numG	numH	numI	numK	numL	numM	numN	numP	numQ	numR	numS	numT	numV	numW	numY	num*	numX	num-	numOther


	# Check for R
	if (check_for_Rscript() ){
		die; 
	}

	# Write the R script, then run it.  
	my $temp_script = $tempdir."/temp_script.R";
	
	my $scriptfh = open_to_write($temp_script, 0, 0, 1);

	my $labels_for_filename = join "_", @labels;
	my $all_output_file  = $save_dir.'/'.$labels_for_filename.".$type.freq.pvalue.all.xls";
	my $filt_output_file = $save_dir.'/'.$labels_for_filename.".$type.freq.pvalue.filt.xls";

	# Need to allow user to control "min" value below
	
	my $col_names;
	my @array;
	if ($type eq "nuc"){
		@array = qw(A C G T);
	}
	elsif($type eq "codon"){
		#$col_names = join "\",\"num", qw(GCA:A GCC:A GCG:A GCT:A TGC:C TGT:C GAC:D GAT:D GAA:E GAG:E TTC:F TTT:F GGA:G GGC:G GGG:G GGT:G CAC:H CAT:H ATA:I ATC:I ATT:I AAA:K AAG:K CTA:L CTC:L CTG:L CTT:L TTA:L TTG:L ATG:M AAC:N AAT:N CCA:P CCC:P CCG:P CCT:P CAA:Q CAG:Q AGA:R AGG:R CGA:R CGC:R CGG:R CGT:R AGC:S AGT:S TCA:S TCC:S TCG:S TCT:S ACA:T ACC:T ACG:T ACT:T GTA:V GTC:V GTG:V GTT:V TGG:W TAC:Y TAT:Y);  # Copied from an existing codon frequency file header.  Removed codons with Ns:  GCN:A GGN:G CTN:L CCN:P CGN:R TCN:S ACN:T GTN:V		# Colons and * are converted to '.' by read.delim in R.  Stop codons will look like "numTGA.."
		@array = qw(GCA.A GCC.A GCG.A GCT.A TGC.C TGT.C GAC.D GAT.D GAA.E GAG.E TTC.F TTT.F GGA.G GGC.G GGG.G GGT.G CAC.H CAT.H ATA.I ATC.I ATT.I AAA.K AAG.K CTA.L CTC.L CTG.L CTT.L TTA.L TTG.L ATG.M AAC.N AAT.N CCA.P CCC.P CCG.P CCT.P CAA.Q CAG.Q AGA.R AGG.R CGA.R CGC.R CGG.R CGT.R AGC.S AGT.S TCA.S TCC.S TCG.S TCT.S ACA.T ACC.T ACG.T ACT.T GTA.V GTC.V GTG.V GTT.V TGG.W TAC.Y TAT.Y TAA.. TAG.. TGA..);
	}
	elsif($type eq "aa"){
		@array = qw(A C D E F G H I K L M N P Q R S T V W Y .);	# stop codon symbol * is converted to '.' by read.delim in R.  'num*' becomes 'num.'.  'num-' becomes 'num..1' in R.
	}
	else {
		print STDERR "Type not found: $type!\n";
		exit;
	}
	$col_names = join "\",\"num", @array;
	$col_names = "\"num" . $col_names . "\""; 
	
	
#	print "colnames:\n$col_names\n"; exit;
	
	print $scriptfh '
fisherbyresidue <- function(d, min=50){
	p.value = NA;
';
	print $scriptfh "col.names = c($col_names);\n";		# Need to make applicable to amino acid too. All num columns except num*	numX	num-	numOther for aa.  For nuc, all num except numN	numOther.  Problem is that numN is valid in aa.  
	print $scriptfh 'check=apply(d[,col.names],2,max)
	num=length(check[check!=0])
	if (all(sum(d[,col.names]) >= min) && (num > 1) && (nrow(d[,col.names]) > 1)) {		
		temp <- fisher.test(d[,col.names], simulate.p.value=T, B=100000);
		p.value = temp$p.value;
	}
	return(p.value);	
}
'."\n";

	print $scriptfh "data = read.delim(\"$merged_file\",fill=TRUE,header=TRUE,sep=\"\\t\",stringsAsFactors=FALSE)\n";
	print $scriptfh "print(\"Sample table:\")\n";
	print $scriptfh "print(table(data\$Sample))\n";
	print $scriptfh "nsamples <- nrow(table(data\$Sample))\n";
	print $scriptfh "if (nsamples < 2){\n";
	print $scriptfh "  stop(\"Only one sample found so no comparisons can be made. Quitting...\")\n";
	print $scriptfh "}\n\n";	
	
	print $scriptfh 'data.list = split(data,data[,1])'."\n";
	print $scriptfh 'y=unlist(lapply(data.list,fisherbyresidue))'."\n";
	print $scriptfh 'data.sample.list = split(data,data$Sample)'."\n";
	
	my $out_col_names = join "\",\"", @header[0..7]; # qw(name gene	aminoAcidPosition	refAminoAcid	consensusAminoAcid	refDiffConsensus	Sample	coverageDepth);
	$out_col_names = "\"" . $out_col_names . "\""; 
	$out_col_names = $out_col_names . "," . $col_names;
	# Need to maybe add more columns here...
#	print "colnames:\n$out_col_names\n"; exit;
	print $scriptfh "out.col.names = c($out_col_names)\n";
	
	print $scriptfh 'z=merge(data.sample.list[[1]],data.sample.list[[2]][,out.col.names],by=c(1,2,3),sort=FALSE)'."\n";
	if (scalar (@labels) > 2){
		print $scriptfh 'z=merge(z,data.sample.list[[3]][,out.col.names],by=c(1,2,3),sort=FALSE)'."\n";
	}
	if (scalar (@labels) > 3){
		print $scriptfh 'z=merge(z,data.sample.list[[4]][,out.col.names],by=c(1,2,3),sort=FALSE)'."\n";
			# Haven't tested this.  If it works, could be put into a loop based on number of samples in @labels	
	}
	print $scriptfh 'z=merge(z,cbind(names(y),p.value=y),by.x=1,by.y=1, sort=FALSE)'."\n";
	print $scriptfh 'z$p.value = as.numeric(as.vector(z$p.value))'."\n";
	print $scriptfh 'z$adj.p.value = p.adjust(z$p.value, method="BH",n = sum(!is.na(z$p.value)))'."\n";
	print $scriptfh 'table(z$adj.p.value < 0.05)'."\n"; 	# Not necessary
	print $scriptfh 'filt = z[z$adj.p.value < 0.05 & !is.na(z$adj.p.value),]'."\n";
	print $scriptfh "write.table(filt, file = \"$filt_output_file\", sep=\"\\t\", quote= FALSE, row.names = FALSE)\n";
	print $scriptfh "write.table(z, file = \"$all_output_file\", sep=\"\\t\", quote= FALSE, row.names = FALSE)\n";

	#Make some graphs
#	print $scriptfh 'matplot(z$p.value,z$adj.p.value)'."\n";	# Not necessary
	print $scriptfh 'filt$transf.p.value = -10 * log(filt$adj.p.value, base=10)'."\n";
	print $scriptfh "positions = filt[,c(\"$header[1]\", \"$header[2]\", \"transf.p.value\")]\n";
	print $scriptfh 'pos.list = split(positions,positions[,1])'."\n";
	print $scriptfh "for (name in names(pos.list)) { print(name); pdf(paste(\"$tempdir\", \"/\", name, \"_dot.pdf\", sep=\"\")); dotchart(pos.list[[name]]\$$header[2], labels = pos.list[[name]]\$$header[2], main=paste(name, \"base\", \"changes\"), xlab=\"Position\", cex=0.5); dev.off()}\n";		# Dot plot
	print $scriptfh "for (name in names(pos.list)) { print(name); pdf(paste(\"$tempdir\", \"/\", name, \"_xy.pdf\", sep=\"\")); plot(pos.list[[name]]\$$header[2], pos.list[[name]]\$transf.p.value, main=paste(name, \"base\", \"changes\"), xlab=\"Position\", ylab=\"-10*log(adj.p.value)\", cex=0.5); text(pos.list[[name]]\$$header[2], pos.list[[name]]\$transf.p.value, pos.list[[name]]\$$header[2], pos=4, cex=0.4); dev.off()}\n";
	close ($scriptfh);

	#my $output = `Rscript $temp_script`;
	my $output = system("Rscript $temp_script");
	
	#print STDERR "Rscript output: \n$output\n";
	
	# Merge graphs with gs if no exit code from Rscript
	if($output){
		print STDERR "R script $temp_script exited with code " , $output / 256, "\n";
	}	
	else {
			#One potential problem with this is that if gs is not on the system, individual pdfs will be deleted because they are in the temp directory.   To avoid deleting them, use --keeptmp
		system("gs -q -dBATCH -dNOPAUSE -dAutoRotatePages#/None -sDEVICE=pdfwrite -sOutputFile=$save_dir/$labels_for_filename"."_genes_dot.pdf $tempdir/*_dot.pdf");
		system("gs -q -dBATCH -dNOPAUSE -dAutoRotatePages#/None -sDEVICE=pdfwrite -sOutputFile=$save_dir/$labels_for_filename"."_genes_xy.pdf $tempdir/*_xy.pdf");
	}

	
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
