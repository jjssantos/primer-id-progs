#!/usr/bin/env perl

$| = 1;
use warnings;
use aomisc;

# This script computes the minimum reliable primerID group size.
# Formula is based on simulation using consensus_cutoff.rb script from Zhou et al, J Virol, 2015 with 12bp primerID.

my $usage = "
compute_cutoff.pl <max primerID group size in dataset>
";

die $usage unless $ARGV[0];
my $m = $ARGV[0];
my $cutoff = -6E-22 * $m ** 6 + 0.00000000000000002 * $m ** 5 - 0.0000000000003 * $m ** 4 + 0.000000002 * $m ** 3 - 0.000006 * $m ** 2 + 0.0183 * $m + 2.2794;

my $min = 5;
$cutoff = $min if ($cutoff < $min);
$cutoff = round($cutoff);

print "$cutoff\n";
