#!/usr/local/bio_apps/perl-5.16.2/bin/perl

use strict;
use warnings;
use Getopt::Long;
#use lib::helper_funcs;

my $output_dir = '';
GetOptions('output_dir=s'=>\$output_dir) or die ("Error in command line arguments\n");

die "Please spec an output dir: --output_dir [dir_name]\nDying\n" if $output_dir eq '';
mkdir $output_dir unless (-d $output_dir);

# Mapping and Gap size stats for R1 - R2 spacing (before merging)
# for i in logs/concat_nonoverlap_fastq_*.o*; do echo $i; cat $i | egrep "(mapped to the reference sequence)" ; cat $i | grep -A 10 "^Gap_size" $i; echo "..."; echo; done  > mapping_and_gap_size_stats.txt
my $mapping_and_gap_size_stats_file = dirify($output_dir,'mapping_and_gap_size_stats.txt');
rm_if_exists($mapping_and_gap_size_stats_file);
map {chomp $_; get_num_mappers($_,$mapping_and_gap_size_stats_file)} `ls logs/concatenate_fastq.pl*`;


# PrimerID match summary
my $find_primerID_summary_stats_file = dirify($output_dir,'find_primerID_summary_stats.txt');
rm_if_exists($find_primerID_summary_stats_file);
# for i in logs/filter_fastq_by_primerid_length_*.o*; do echo $i; cat $i | grep "sequences with"; echo; done > find_primerID_summary_stats.txt
map {chomp $_; grep_file_for('sequences with', $_, $find_primerID_summary_stats_file)} `ls logs/filter_fastq_by_primerid_length.pl*`;


# Primer trimming and splitting into regions
my $primer_trimming_split_into_region_stats_file = dirify($output_dir,'primer_trimming_split_into_region_stats.txt');
rm_if_exists($primer_trimming_split_into_region_stats_file);
map {chomp $_; extract_btrim_stats($_,$primer_trimming_split_into_region_stats_file)} `ls logs/Btrim64*`;


# Auto-detected Gap start and size stats 
my $gap_start_length_during_merge_stats_file = dirify($output_dir,'gap_start_length_during_merge_stats.txt');
rm_if_exists($gap_start_length_during_merge_stats_file);
# for i in logs/merge_primerid_groups_*.o*; do echo $i; cat $i | egrep "^(Many reads|Fraction of reads|Most common|Using)"; echo; done > gap_start_length_during_merge_stats.txt
map {chomp $_; grep_file_for('^(Many reads|Fraction of reads|Most common|Using)', $_, $gap_start_length_during_merge_stats_file)} `ls logs/merge_primerid_read_groups.pl.*`;


# Retention based on primerID group size thresholds (also see group plots)
my $merge_primerID_retention_stats_file = dirify($output_dir,'merge_primerID_retention_stats.txt');
rm_if_exists($merge_primerID_retention_stats_file);
# for i in logs/merge_primerid_groups_*.o*; do echo $i; head -50 $i | grep "primerID groups "; done > merge_primerID_retention_stats.txt
map {chomp $_; grep_file_for('primerID groups', $_,$merge_primerID_retention_stats_file)} `ls logs/merge_primerid_read_groups.pl.*`;


# Retention in conversion to amino acid (tossed based on early stop codons)
my $convert_reads_good_tossed_stop_codons_stats_file = dirify($output_dir,'convert_reads_good_tossed_stop_codons_stats.txt');
rm_if_exists($convert_reads_good_tossed_stop_codons_stats_file);
# for i in logs/convert_reads_to_amino_acid_*.o*; do echo -e "$i\t$(cat $i | grep "^Good")\t$(cat $i | grep "^Tossed")"; done > convert_reads_good_tossed_stop_codons_stats.txt
map {chomp $_; grep_file_for('^Good|^Tossed', $_,$convert_reads_good_tossed_stop_codons_stats_file,1)} `ls logs/convert_reads_to_amino_acid.pl*`;

# Retention after majority block filtering.
my $majority_block_filtering_stats_file = dirify($output_dir,'majority_block_filtering_stats.txt');
rm_if_exists($majority_block_filtering_stats_file);

sub get_majority_block_filtering_stats{
    my $file = shift;
    my $cat_into = shift;
    
# out=majority_block_filtering_stats.txt;
 # echo -e "Sample\tRegion\tAll\tMajority\tAllBAM\tMajorityBAM" > $out;
 # for n in 0 1 2;
 # do for s in 30_S1 31_S2 32_S3;
 # do allbam=$(ls ${s}.contigs.pid.btrim.fastq*${n}*.bam);
 # majbam=$(ls ${s}.contigs.pid.btrim.*${n}*.majority.bam);
 # echo -e "$s\t$n\t$(samtools view $allbam | wc)\t$(samtools view $majbam | wc)\t${allbam}\t${majbam}" >> $out;
 # done;
 # done 
}


sub extract_btrim_stats{
    # for i in logs/Btrim*.o*; do echo $i; 
    # cat $i | awk '{
    # if ($0 ~ /Total sequences/ || $0 ~ /^Pattern distribution/){p = 1} 
    # else if ($0 ~ /^$/){p = 0} } 
    # {if (p == 1){print}}' ; echo; done > primer_trimming_split_into_region_stats.txt

    my $file = shift;
    my $cat_into = shift;

    open IN, $file or die "at extract, file : $file\n". $!;    
    open OUT, ">>$cat_into" or die "at extract, cat_into : $cat_into\n". $!;
    print OUT $file."\n";

    my $do_print = 0;
    while (<IN>){
	next unless (($_ =~ /(Total sequences)|(Pattern distribution)/) || $do_print);
	$do_print = 1;
	print OUT $_;
	$do_print = 0 if $_ =~/^\s?$/;
    }
    close IN;
    close OUT;
}

sub grep_file_for{
    my $search_for = shift;
    my $file_name = shift;
    my $cat_into = shift;
    my $single_line = shift;

    open IN, $file_name or die "at grep file, file_name : $file_name\n". $!;    
    open OUT, ">>$cat_into" or die "at grep file, cat_into : $cat_into\n". $!;
    print OUT "\n",$file_name;
    $single_line ? print OUT "\t" : print OUT "\n";
    while (<IN>){
	next unless $_ =~ /$search_for/;
	if ($single_line){
	    chomp $_;
	    print OUT $_."\t";
	}
	else{
	    print OUT $_;
	}
    }
    close IN;
    close OUT;
}

sub rm_if_exists{
    my $file = shift;
    system ("rm $file") if -e $file;
}

sub dirify{
    my $dir = shift;
    my $file = shift;
    return slashify($dir).$file;
}

sub slashify{
    my $str = shift;
    return $str if $str eq '';
    chomp $str;
    return ($str =~ /\/$/) ? $str : $str.'/'
}

sub get_num_mappers{
    my $max_lines = 10;		#    only want 10 lines for some reason
    my $file = shift;
    my $write_to = shift;

    my $i = 0;
    open IN, $file or die "at get_num_mappers() $file\n";
    open OUT, ">>$write_to" or die "at get_num_mappers() write_to $write_to\n";
    print OUT $file,"\n";
    while (<IN>){
	next unless (($_ =~ /^Gap_size/) || $i);

	$i++;
	print OUT $_;
	if ($i > $max_lines){
	    print OUT "...\n";
	    last;
	}
    }
    close IN or warn $!;
    close OUT or warn $!;
}
