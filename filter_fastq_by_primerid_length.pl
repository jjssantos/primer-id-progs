#!/usr/local/bio_apps/perl-5.16.2/bin/perl
use strict;
use warnings;
$| = 1;

use Bio::SeqIO; 
use Bio::Seq::Quality;
use File::Basename;
use lib dirname (__FILE__);	# # philip macmenamin
use aomisc;
use Getopt::Long;
use Data::Dumper;
use Statistics::R;


if (@ARGV){		
	print STDERR "Arguments: ", join " ", @ARGV, "\n";	
}


my $pre = "";
my $post;
my $stdout;
my $primerid_min = 1;
my $removepost;
my $n; 				# philip
my $file ='';			# philip
my $output_dir = './';	# default.
my $r2 = '';

GetOptions("n=i"=>\$n,
	   "output_dir=s"=>\$output_dir,
	   'pre=s' => \$pre, 
	   'post=s' => \$post, 
	   'stdout' => \$stdout, 
	   'primerid_min=s' => \$primerid_min, 
	   'removepost' => \$removepost, 
	   "file_in=s" => \$file,
	   'r2=s' => \$r2,
	   );

my $usage= basename($0)." looks for a primerID or barcode of a defined length in 
the 5' end of a read.  If found, the primerID is trimmed and placed in the seqid in the 
output fastq file.  Input fastq file(s) must be unzipped.  It processes about 750 reads per 
second.
filter_fastq_by_primerid_length.pl [OPTIONS] --file_in <fastq> --n <size> --post <seq>
OPTIONS:
--file_in	Input fastq file.  Single-end fastq file, or R1 of paired-end fastq file. 
	Required.  
--r2	Input fastq file.  R2 of paired-end fastq file.  Optional.  If provided, the 
	primerID sequence will be added to the header of this file as well.
--n     Size of primerID.  Required.
--post	Sequence following the primerID.  Default, any sequence. 
	Required!  Otherwise, the length of the primerID is meaningless.
--output_dir	Output directory to save the files. Default current directory.
--pre	Sequence preceding the primerID.  Default, no sequence before primerID (i.e., 
	primerID starts at extreme 5' end).
--removepost	Remove the \"post\" sequence in addition to the primerID.  Default is to keep the 
	post sequence.
e.g., filter_fastq_by_primerid_length.pl --pre GC --post AAGCAGT --file_in myfile.fastq --n 8 2> out.err 
	Outputs the files myfile.pid.fastq and myfile.pid.primerids.gz
";
# --stdout	Print to STDOUT.

# To do
# Allow input of a second fastq file (R2).  If supplied, the primerid will be placed in the header for that one too.
# Printout graph of primerID counts (not just table)
# Make option for --barcode_max 

# Change Log
# 2013-09-12 
# Added progress tracker 
# Trim quality scores for the same bases that are trimmed for barcode
# Remove post sequence option --removepost.
# Removed this option (put the functionality in a separate script, filter_fastq_by_primerid_count.pl)
#--barcode_min 	Minimum number of barcoded molecules.  Default = 1 (all retained).  It is 
#	highly suggested to set this to at least 2 or 3.  
# Renamed script from filter_fastq_by_barcode_length.pl to filter_fastq_by_primerid_length.pl
# 2014-01-09
# Changed output file name suffix from primerid10bp to pid10bp
# Updated the usage.  
# Make reference to the barcode as "primerID" instead of "barcode"
# Removed section at the end for filtering by min/max # reads per group (after "exit;").  This is now in merger_primerid_read_groups.pl
# 2014-06-03
# Print out plot of read counts based on the primerID group size they are assigned to. 

# philip = 
# either we use getopts or not. using argv and get opts is confusing. For now I'm getting rid of
# argv[1] and calling it n (because it's a number?)

# 2015-04-16
# Updated usage to match Philip's changes.
# 2015-06-19
# Updated plot_counts so that we are plotting numbers and fractions of reads instead of primerID groups. Also made the y-axis fit the data better.
# 2015-11-06
# Added --r2 option.  



die $usage unless ($n);

unless ($post){
	warn "Please enter value for --post.\n";
	exit 1;
}
unless ($output_dir){
	warn "Please enter an output directory.\n";
	exit 1;
}

mkdir $output_dir unless(-e $output_dir);


if ($r2 && $stdout){
	print STDERR "Can't use --stdout with paired-end fastq files.  Will write output for R1 and R2 to separate files.\n";
}

my $start_time = time;

$output_dir .= "/" unless ($output_dir =~ m/\/$/);	# Add '/' to the end of the directory name if it doesn't exist already.

my $newfilename = get_new_filename($file);

my $in=Bio::SeqIO->new(-format => 'fastq', 
	-variant => 'sanger', 
#	-file => $ARGV[0],
	-file => $file,
	); 
	
my $out;
if ($stdout){
	$out = Bio::SeqIO->new(-format => 'fastq',
	-fh => \*STDOUT,
	);
} else {
	$out = Bio::SeqIO->new(-format => 'fastq',
	-file => ">$newfilename",
	);
}

# Set up FASTQ reader and output file for R2
my $r2_in;
my $r2_out;
if ($r2){
	my $r2_newfilename = get_new_filename($r2);
	$r2_in=Bio::SeqIO->new(-format => 'fastq', 
		-variant => 'sanger', 
		-file => $r2,
		); 
	$r2_out = Bio::SeqIO->new(-format => 'fastq',
		-file => ">$r2_newfilename",
		);
}

my $count = 0;
my $primerid_tally;
my $i = 0;
my $with_primerID = 0;
my $no_primerID = 0;
while(my $seq = $in->next_seq){
	my $r2_seq = $r2_in->next_seq if $r2;
	if ($seq->seq =~ m/^$pre([ACTG]{$n})($post)([ACTGN]+)$/){ 
		my $length = length($seq->seq);
		$seq->id($seq->id .':'."$1");	# Add the primerid to the id.  This is important for the script merge_primerid_read_groups.pl
		$r2_seq->id($r2_seq->id .':'."$1") if $r2;	# Add the primerid to the id for the R2 read as well, if R2 file provided.  No other modifications required for R2 read (i.e., no trimming)
		my $new_seq;
		my $start = 1;
		if ($removepost){
			$new_seq = $3;
			$start += length($pre) + $n + length($post);
		}
		else {
			$new_seq = $2 . $3;
			$start += length($pre) + $n;
		}
#		$new_seq = $seq->trunc($start,$length)
		$seq->seq($new_seq);		# replace the sequence with just the sequence after the primerid.
		
		my $qual = $seq->qual;	# Array ref of quality scores (integers)
#					print Dumper($qual);
#					print scalar(@$qual)."\n";
		my @subqual = @$qual[-length($seq->seq)..-1];	# Get subset of quality scores to match the sequence
#					print join "|", @subqual; 
#					print "\n".scalar(@subqual)."\n";
#		my $deleted_length = length($pre) + $n;
#		$deleted_length += length($post) if ($removepost);
#		my $start = $deleted_length + 1;
#		print "start: $start original_length: $length new_length: ".length($seq->seq)."\n";
#		$seq->qual( $seq->subqual( $start, ) );	# Set new quality score 
		$seq->qual( \@subqual );		# Set new quality score 
		$seq->force_flush("true");
		print STDERR "not flush\n" unless ($seq->quality_is_flush);
#					print "length of seq: ".length($seq->seq)." supposed length of qual: ".scalar(@subqual)." length of qual: ".scalar(@{$seq->qual})."\n";
		$out->write_fastq($seq);		# Print out the fastq record
		$r2_out->write_fastq($r2_seq) if $r2;	# Print out the fastq record for R2 as well (if R2 file provided)
		$primerid_tally->{$1}->{tally}++;
		push @{$primerid_tally->{$1}->{list}}, $seq->id;
		$with_primerID++;
	} 
	else {
		$no_primerID++;
	}
	$i++;
	if (($i % 10000) == 0){	
		print STDERR "Done processing $i sequences. ";
		&elapsed($start_time, ' Elapsed');
	}
}

#print Dumper($primerid_tally);

# Print out the table of counts
my $counts; 
print STDERR "primerID_counts\tNumber_primerid_groups_for_each_level\n";
foreach my $primerid (keys %$primerid_tally){
	$counts->{$primerid_tally->{$primerid}->{tally}}++;
} 
$counts->{0} = $no_primerID;		# reads with no primerID, or "zero" 

foreach my $num (sort {$a <=> $b} keys %$counts){ 
	print STDERR "$num\t$counts->{$num}\n";
}

my ($filename,$dir,$ext) = fileparse($newfilename,@SUFFIXES);
plot_counts($filename, $counts);		

printf STDERR "\n%i sequences with primerID\t\t%.3f", $with_primerID, $with_primerID/$i;
printf STDERR "\n%i sequences without primerID\t\t%.3f\n", $no_primerID, $no_primerID/$i;

# Now print out the primer IDs with the individual counts and the read ids for each.  
my $primerID_filename = $output_dir.$filename.'.pid.primerids';
my $writefh = open_to_write($primerID_filename, 1);			# gzip compression for second argument
foreach my $ID (keys %$primerid_tally){
	print $writefh "$ID\t";
	print $writefh "$primerid_tally->{$ID}->{tally}\t";
	my $ids = "";
	foreach (@{$primerid_tally->{$ID}->{list}}){
		$ids .= $_ . ",";
	}
	$ids =~ s/,$//;
	print $writefh "$ids\n";
}


&elapsed($start_time, 'Total');

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub plot_counts {
	my ($newfilebase,$counts) = @_;
	## Graph Primer ID group size distribution

	# Before filtering the $read_tally data structure, make a graph of the primerID groups by number of reads per group, and then how many groups have that number.  
	# For now just make a simple barplot.  Single sample.  X-axis is "Counts in each PrimerID group" and Y-axis is fraction of Primer IDs
	# It would be nice to make a summary plot with all samples in the workflow as well.  
	# So, I should make a temporary file that can be read into R and parsed later as well.
	# Does it make sense to have the "zero" category?  The y-axis is is number of primerID groups, which doesn't really fit with the group of reads with no primerID.  For the "zero" group, the y-axis should be # of reads.
	# It would make more sense to have *reads* for the y-axis
	
	my $group_count_file = 	$newfilebase . ".read.counts.by.group.txt";
	my $group_count_graph = $newfilebase . ".read.counts.by.group.barplot.pdf";
	my $counts_reads; 	# $counts is the tally of primerID groups, but $counts_reads will be the number of reads. 
	foreach my $group_size (keys %$counts){
		if ($group_size == 0){
			$counts_reads->{$group_size} = $counts->{$group_size};
		}
		else {
			$counts_reads->{$group_size} = $counts->{$group_size} * $group_size;
		}
	}
	
	my $total_reads = total(values %$counts_reads);
	my $count_fh = open_to_write($group_count_file);
	print $count_fh "PrimerID_group_size\tNumber_of_reads\tFraction_of_reads\n";
	foreach my $group_size (sort {$a <=> $b} keys %$counts){
		my $number_of_reads = $counts->{$group_size};
		$number_of_reads = $number_of_reads * $group_size if ($group_size > 0);
		printf $count_fh "%2d\t%2d\t%.3f\n", $group_size, $number_of_reads, $number_of_reads / $total_reads;
	}
	close ($count_fh);

	# Now graph this in R
	my $R = Statistics::R->new();		# http://search.cpan.org/~fangly/Statistics-R-0.31/lib/Statistics/R.pm
	my $out1 = $R->run(
		"table <- read.delim(\"$group_count_file\")",
		"ymax <- min(max(table\$Fraction_of_reads) * 1.25,1)",		# In plot, change ylim to lesser value of (max Fraction_of_reads * 1.25 or 1) in case all values are ~0.2 or 0.1 (ymax of 1 is too high in that case so the plot doesn't look good)
		"pdf(\"$group_count_graph\")",
		"barplot(table\$Fraction_of_reads, names.arg=table\$PrimerID_group_size, ylim=c(0,ymax), ylab=\"Fraction of Reads\", xlab=\"PrimerID Group Size (Number of Reads)\", main=\"Read Distribution by PrimerID Group for $newfilebase\", cex.main=0.9)",
#		"library(ggplot2)",
#		"ggplot(table, aes(x = Reads_in_PrimerID_group, y = Fraction_of_groups)) + geom_bar(stat = \"identity\")",		# Testing, so that I can add labels above the bars for the actual numbers
		"dev.off()"
	);
	
}
#-----------------------------------------------------------------------------
sub get_new_filename {
	# This subroutine takes a file as input.
	# 1. Parses out the basename, extension
	# 2. Adds '.pid'
	# 3. Shifts R1/R2 to the end of the file name (if it exists)
	# 4. Adds new output directory and adds back extension
	# 5. Returns the new filename
	my $file = shift;

	my ($filename,$dir,$ext) = fileparse($file,@SUFFIXES);		# fileparse($file,qr/\.[^.]*/);
	
	$filename .= '.pid';
	
	if ($filename =~ m/(_R[12])\.?/){		# -f 24_S1_R1.fastq -r 24_S1_R2.trim.fastq	# fileparse removes the .fastq part
		# For PE FASTQ files
		# Put R1/R2 just before the extension.  e.g., 24_S1_R2.trim.fastq output file will be 24_S1.trim.pid_R2.fastq; 24_S1_R1.fastq output file will be 24_S1.pid_R1.fastq
		my $r = $1;
		$filename =~ s/$r//;
	#	$r =~ s/_/\./;
		$filename .= $r;	
							
	}
#	else {
#		print STDERR "no match in file! $filename\n";
#	}
	
	my $newfilename = $output_dir.$filename.$ext;

	print STDERR "new file: $newfilename\n"; 
	
	return $newfilename;

}
#-----------------------------------------------------------------------------
