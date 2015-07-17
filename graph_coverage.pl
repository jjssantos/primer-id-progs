#!/usr/local/bio_apps/perl-5.16.2/bin/perl
$| = 1;
use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use aomisc;
use Getopt::Long;

use File::Temp;
use File::Basename;


#my $file = $ARGV[0];
my $PWD = pwd_for_hpc();
my $usage = basename($0)." <coverage file or BAM>
Coverage file is the output of genomeCoverageBed with option -d.
Or, if input file is BAM, genomeCoverageBed will be run first with options -split -d.
If input filename is /path/file.bam, output will be /path/file.pdf.
";

# Change log
# 2013-11-04 Modified to plot each chromosome separately

# die "$usage" unless $ARGV[0];
my %options;
GetOptions(\%options,
	   'bam_file=s',
	   'output_dir=s'
	   );

die "I need a bam file to graph \n". $! unless ($options{bam_file} && (-e $options{bam_file}));
die "I need an output_dir \n". $! unless ($options{output_dir});
$options{output_dir} =~ s/\/$//;
system ("mkdir $options{output_dir}") unless -e $options{output_dir};
die "failed to make output dir $options{output_dir}\n".$! unless -d $options{output_dir};

# Check for R
if (check_for_Rscript() ){
	die; 
}

my ($filename,$dir,$ext) = fileparse($options{bam_file},@SUFFIXES);		

my $coverage_file = $options{bam_file};

if ($ext =~ m/bam/i){
	# Then run genomeCoverageBed first
	# E.g.,	for i in {1..18}; do b=$(ls Sample*_${i}_*max10.bam); echo $b; genomeCoverageBed -ibam $b -split -d > ${b/.bam/.cov}; done
	$coverage_file = $options{output_dir} . "/" . $filename . '.cov'; 
	#print STDERR "Coverage file: $coverage_file\n";
	my $cmd = pwd_for_hpc()."/genomeCoverageBed -ibam $options{bam_file} -split -d > $coverage_file";
	system($cmd);
}


# Write the R script, then run it.  
my $tempdir = File::Temp->newdir(  "graph_coverage_tempXXXXX" );
my $temp_script = $tempdir."/temp_script_for_plot.R";

my $scriptfh = open_to_write($temp_script, 0, 0, 1);

my $title = $filename;
$title =~ s/\.cov\w+$//i;
my $pdf = $options{output_dir} . "/" . $filename . '.pdf'; 

print $scriptfh "tally = read.delim(\"$coverage_file\", header=FALSE, na.strings=\"NULL\")\n";
print $scriptfh "names = as.character(unique(tally\$V1))\n";
print $scriptfh "pdf(\"$pdf\")\n";
print $scriptfh "for (i in names) {\n";
print $scriptfh "a <- subset(tally, tally\$V1 == i)\n";
print $scriptfh "plot(a\$V2, a\$V3, main=paste(i, \"$title Coverage Plot\"), xlab=\"Reference Position\", ylab=\"Coverage\", pch=19, type=\"l\", lab=c(20,20,10), cex.axis=0.4)\n";		# xlim=c(1,$length)
print $scriptfh "}\n";
print $scriptfh "dev.off()\n";

close ($scriptfh);

my $output = `Rscript $temp_script`;


