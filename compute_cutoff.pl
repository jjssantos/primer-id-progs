#!/usr/bin/env perl

use lib "/nethome/oleraj/lib";
$| = 1;
use warnings;
use aomisc;

my $usage = "
compute_cutoff.pl <max primerID size>
";

die $usage unless $ARGV[0];
my $m = $ARGV[0];
my $cutoff = -6E-22 * $m ** 6 + 0.00000000000000002 * $m ** 5 - 0.0000000000003 * $m ** 4 + 0.000000002 * $m ** 3 - 0.000006 * $m ** 2 + 0.0183 * $m + 2.2794;

my $min = 5;
$cutoff = $min if ($cutoff < $min);
$cutoff = round($cutoff);

print "$cutoff\n";
