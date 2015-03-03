#!/usr/local/bio_apps/perl-5.16.2/bin/perl

##### #!/usr/bin/env perl

# Takes a bed file as STDIN and outputs a single bed region for the majority sequence.
# E.g., bamToBed -i $bam | get_majority_start_stop.pl > majority.bed
use File::Basename;
use lib dirname (__FILE__);	# # philip macmenamin
use aomisc;

# use lib './specific_progs';		# philip macmenamin
# use aomisc; 
my $hash; 
while(<>){
	my @F = split(/\t/); 
	$hash->{start}->{$F[1]}++; 
	$hash->{end}->{$F[2]}++; 
	$chr = $F[0];
} 
my $start = find_key_with_biggest_value($hash->{start}); 
my $end = find_key_with_biggest_value($hash->{end}); 
print "$chr\t$start\t$end\n";

