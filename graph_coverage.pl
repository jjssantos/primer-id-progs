#!/usr/bin/env perl
$| = 1;
use strict;
use warnings;
use aomisc;
use File::Temp;
use File::Basename;


my $file = $ARGV[0];

my $usage = basename($0)." <coverage file or BAM>
Coverage file is the output of genomeCoverageBed with option -d.
Or, if input file is BAM, genomeCoverageBed will be run first with options -split -d.
If input filename is /path/file.bam, output will be /path/file.pdf.
";

# Change log
# 2013-11-04 Modified to plot each chromosome separately

die "$usage" unless $ARGV[0];


# Check for R
if (check_for_Rscript() ){
	die; 
}

my ($filename,$dir,$ext) = fileparse($ARGV[0],@SUFFIXES);		

my $coverage_file = $ARGV[0];

if ($ext =~ m/bam/i){
	# Then run genomeCoverageBed first
	# E.g.,	for i in {1..18}; do b=$(ls Sample*_${i}_*max10.bam); echo $b; genomeCoverageBed -ibam $b -split -d > ${b/.bam/.cov}; done
	$coverage_file = $dir . "/" . $filename . '.cov'; 
	print STDERR "Coverage file: $coverage_file\n";
	my $cmd = "genomeCoverageBed -ibam $ARGV[0] -split -d > $coverage_file";
	system($cmd);
}



# Write the R script, then run it.  
my $tempdir = File::Temp->newdir(  "graph_coverage_tempXXXXX" );

my $temp_script = $tempdir."/temp_script_for_plot.R";
my $scriptfh = open_to_write($temp_script, 0, 0, 1);

my $title = $filename;
$title =~ s/\.cov\w+$//i;
my $pdf = $dir . "/" . $filename . '.pdf'; 

print $scriptfh "tally = read.delim(\"$coverage_file\", header=FALSE)\n";
print $scriptfh "names = as.character(unique(tally\$V1))\n";
print $scriptfh "pdf(\"$pdf\")\n";
print $scriptfh "for (i in names) {\n";
print $scriptfh "a <- subset(tally, tally\$V1 == i)\n";
print $scriptfh "plot(a\$V2, a\$V3, main=paste(i, \"$title Coverage Plot\"), xlab=\"Reference Position\", ylab=\"Coverage\", pch=19, type=\"l\", lab=c(20,20,10), cex.axis=0.4)\n";		# xlim=c(1,$length)
print $scriptfh "}\n";
print $scriptfh "dev.off()\n";

close ($scriptfh);

my $output = `Rscript $temp_script`;


