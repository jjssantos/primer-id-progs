#!/usr/bin/env perl
#	#!/usr/local/bio_apps/perl-5.16.2/bin/perl
#Add use lib statement to assume there is a directory at the same level as bin in which the script is run, called 'lib'
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin";

use strict;
use warnings;
use Getopt::Long;
use aomisc;

#my $BWA='/usr/local/bio_apps/bwa/bwa';
my $PWD = pwd_for_hpc();
my $BWA = $PWD . "/bwa";

my $p=8;						# Number of threads to use for BWA. 

my %options;
$options{ref} = '';		# eg data_default/Seq12_093009_HA_cds.fa

GetOptions(\%options,
           'ref=s',
	   'output_dir=s'
	   );

die "I need a ref file to index (-ref)\n". $! if (($options{ref} eq '') || (! -e $options{ref}));
die "I need an output_dir (-output_dir)\n". $! if ($options{output_dir} eq '');
$options{output_dir} =~ s/\/$//;

system ("mkdir $options{output_dir}") unless -e $options{output_dir};
die "failed to make output dir $options{output_dir}\n".$! unless -d $options{output_dir};

system ("cp $options{ref} $options{output_dir}");
my $ref = $options{output_dir}.'/'.extract_file_name($options{ref});
die "failed to copy files properly $ref\n".$! unless (-e $ref);

# Check to see that the reference is indexed.  If not, attempt to index it
system ("$BWA index $ref") unless -e $ref.'.bwt';
die "I tried to index $ref but failed.".$! unless  -e $ref.'.bwt';

sub extract_file_name{
    my $file = shift;
    my @a = split /\//,$file;
    return $a[$#a];
}
