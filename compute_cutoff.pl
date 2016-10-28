#!/usr/bin/env perl

$| = 1;
use warnings;
use aomisc;


my $usage = "
compute_cutoff.pl <max primerID group size in dataset>
This script computes the minimum reliable primerID group size.
The formula is based on a simulation using the consensus_cutoff.rb 
script from Zhou et al, J Virol, 2015 for a 12bp primerID.
";

die $usage unless $ARGV[0];
my $m = $ARGV[0];
my $cutoff = -6E-22 * $m ** 6 + 0.00000000000000002 * $m ** 5 - 0.0000000000003 * $m ** 4 + 0.000000002 * $m ** 3 - 0.000006 * $m ** 2 + 0.0183 * $m + 2.2794;

my $min = 5;		# Modify this to your personal minimum acceptable threshold.  Suggested 3 to 5.  Even if $cutoff is below this, the script won't report a value below this minimum.

$cutoff = $min if ($cutoff < $min);
$cutoff = round($cutoff);

print "$cutoff\n";
