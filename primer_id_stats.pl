#!/usr/local/bio_apps/perl-5.16.2/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#use lib::helper_funcs;

my $output_dir = '';
my $meta_file = '';
my $SAMTOOLS = '/usr/local/bio_apps/samtools/samtools';

GetOptions('output_dir=s'=>\$output_dir,
	   'meta_file=s'=>\$meta_file
	  ) or die ("Error in command line arguments\n");

die "Please spec an output dir: --output_dir [dir_name]\nDying\n" if $output_dir eq '';
die "Please spec meta file: --meta_file [file_name]\nDying\n" if $meta_file eq '';

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


# Retention based on primerID group size thresholds (also see group plots) and ambiguous characters outside of the gap region
my $merge_primerID_retention_stats_file = dirify($output_dir,'merge_primerID_retention_stats.txt');
rm_if_exists($merge_primerID_retention_stats_file);
# for i in logs/merge_primerid_groups_*.o*; do echo $i; head -50 $i | grep "primerID groups "; done > merge_primerID_retention_stats.txt
map {chomp $_; grep_file_for('(primerID groups|Consensus reads) ', $_,$merge_primerID_retention_stats_file)} `ls logs/merge_primerid_read_groups.pl.*`;	 # Andrew


# Plot position of ambiguous nucleotides in consensus reads
map {chomp $_; plot_pos_ambig_nucs($_)}  `ls logs/merge_primerid_read_groups.pl*`;


# Retention in conversion to amino acid (tossed based on early stop codons)
my $convert_reads_good_tossed_stop_codons_stats_file = dirify($output_dir,'convert_reads_good_tossed_stop_codons_stats.txt');
rm_if_exists($convert_reads_good_tossed_stop_codons_stats_file);
# for i in logs/convert_reads_to_amino_acid_*.o*; do echo -e "$i\t$(cat $i | grep "^Good")\t$(cat $i | grep "^Tossed")"; done > convert_reads_good_tossed_stop_codons_stats.txt
map {chomp $_; grep_file_for('^Good|^Tossed', $_,$convert_reads_good_tossed_stop_codons_stats_file,1)} `ls logs/convert_reads_to_amino_acid.pl*`;


# Retention after majority block filtering.
my $majority_block_filtering_stats_file = dirify($output_dir,'majority_block_filtering_stats.txt');
rm_if_exists($majority_block_filtering_stats_file);
cat_into($majority_block_filtering_stats_file, "Sample\tRegion\tAll\tMajority\tAllBAM\tMajorityBAM\n");
my $samples = hash_file($meta_file,0,3);
get_majority_block_filtering_stats($samples,$majority_block_filtering_stats_file);

# Linked Variant pairs with various threshold FDR values:
my $summary_linked_variants_by_type_and_FDR_stats_file = dirify($output_dir,'summary_linked_variants_by_type_and_FDR_stats.txt');
rm_if_exists($summary_linked_variants_by_type_and_FDR_stats_file);
cat_into($summary_linked_variants_by_type_and_FDR_stats_file, "## Note: each threshold column represents variants *in addition to* the previous columns.\n");
map {chomp $_; do_summary_linked_variants($_,$summary_linked_variants_by_type_and_FDR_stats_file)} `ls outputs/calculate_linkage_disequilibrium.pl/*btrim.*linkage.minfreq0.0*.xls`;


# Enriched variants in treatment vs. control with various FDR values:
my $summary_compare_variants_treatment_control_by_type_and_FDR_stats_file = dirify($output_dir,'summary_compare_variants_treatment_control_by_type_and_FDR_stats.txt');
rm_if_exists($summary_compare_variants_treatment_control_by_type_and_FDR_stats_file);
cat_into($summary_compare_variants_treatment_control_by_type_and_FDR_stats_file, "## Note: each threshold column represents variants *in addition to* the previous columns.\n");
map {chomp $_; by_type_and_FDR_stats($_,$summary_compare_variants_treatment_control_by_type_and_FDR_stats_file)  } `ls outputs/compare_variant_frequency.pl/Passage_Parent.*.freq.pvalue.all.xls`;


sub plot_pos_ambig_nucs{
# for i in logs/merge_primerid_groups_*.o*;
#  do base=$(basename $i);
#  table=${base%.*}.ambigpos.txt;
#  cat $i | egrep -A 10000 "^#?Position" | grep -v "^Total time" | sed 's/^#//g' | awk '{if($3 ~ /GAP/)'{print $1"\t"$2"\tT"} else if($1 ~ /Position/){print $1"\t"$2"\tGap"} else {print $1"\t"$2"\tF" } }' > $table;
# done
# for i in *ambigpos.txt;
#  do ./graph_ambig_pos.R $i;
#  done

    my $file = shift;
    my $bn = basename($file);
    my $i = 0;
    open IN, $file or die $!, "Oops at plot_pos_ambig_nucs\n";
    my @lines;
    while (<IN>){
	$i = 1 if (($_ =~ /#?Position/) && ($i ==0));
	next unless $i;
	last if $i > 10000;
	$i++;
	next unless $_ !~ /^Total time/;
	push @lines, t_or_f($_);       
    }

    my $table_file = dirify($output_dir,$bn.'.ambigpos.txt');
    open OUT, ">$table_file" or die $!;
    print OUT join "\n", @lines;
    close OUT;
}

sub t_or_f{
    my $line = shift;
    chomp $line;
    $line =~ s/^#//g;
    my @l = split /\t/,$line;
    my @l_ret;
    if ($l[2] && ($l[2] =~ /GAP/) ){
     	@l_ret = ($l[0],$l[1],'T');
    }
    elsif($0 =~ /Position/){
    	@l_ret = ($l[0],$l[1],'Gap')
    }
    else {
    	@l_ret = ($l[0],$l[1],'F')
    }
    return join "\t",@l_ret;
}

sub by_type_and_FDR_stats{

# # Enriched variants in treatment vs. control with various FDR values:
# out=summary_compare_variants_treatment_control_by_type_and_FDR_stats.txt;
#  echo "## Note: each threshold column represents variants *in addition to* the previous columns." > $out;
#  for type in nuc aa codon;
#  do for i in compare_variant_freq_region_*/Passage_Parent.${type}.freq.pvalue.all.xls;
#  do if [ -e $i ];
#  then ls $i;
#  cat $i | grep -v "^name" | perl -se'my $hash;
#  print "## Note: each threshold column represents variants in addition to the previous columns.\n";
#   LINE: while(<>){ my @F = split(/\t/);
#  $hash->{$a}->{'total'}++;
#   foreach my $t (0.0001, 0.001, 0.01, 0.05, 0.1){ if ($F[-1] < $t){ $hash->{$a}->{$t}++;
#  next LINE;
#  } } } print "#Type\t<0.0001\t<0.001\t<0.01\t<0.05\t<0.10\ttotal\n";
#  foreach my $type (qw(nuc codon aa)){ next unless (exists($hash->{$type}));
#  print "$type";
#  foreach my $t (0.0001, 0.001, 0.01, 0.05, 0.10, 'total'){ my $c = $hash->{$type}->{$t} ? $hash->{$type}->{$t} : 0;
#  print "\t$c";
#  } print "\n";
# }' -- -a=$type;
#  fi;
#  done;
#  done >> $out

    my $file = shift;
    my $out = shift;
    my %hash;
    open IN,$file or die $!;
    while(<IN>){
	chomp $_;
	next if $_ =~ /^#|^\s?$/;
	next unless $_ =~ /^name/;
	my @F = split(/\t/);
	$hash{'total'}++;
	foreach my $t (0.0001, 0.001, 0.01, 0.05, 0.1){
	    next unless ($F[-1] < $t);
	    $hash{$t}++;
	}
    }
    close IN;
    cat_into($out,$file."\n");

    print join "\t", keys %hash;
    print "\n";
    print join "\t", values %hash;
    print "\n";

    dump_n_n_total($out,\%hash);
}
# ' -- -a=$type;

sub do_summary_linked_variants{

# Linked Variant pairs with various threshold FDR values:
# out=summary_linked_variants_by_type_and_FDR_stats.txt;
#  echo "## Note: each threshold column represents variants *in addition to* the previous columns." > $out;
#  for i in *btrim.*linkage.minfreq0.0*.xls;
#  do if [ -e $i ];
#  then ls $i;
#  cat $i | grep -v "^#group" | perl -e'my $hash;
#  LINE: while(<>){ my @F = split(/\t/);
#  $hash->{'total'}->{$F[2]}++;
#   foreach my $t (0.0001, 0.001, 0.01, 0.05, 0.1){ if ($F[-1] < $t){ $hash->{$t}->{$F[2]}++;
#  next LINE;
#  } } } print "#Type\t<0.0001\t<0.001\t<0.01\t<0.05\t<0.10\ttotal\n";
#  foreach my $type (qw(nuc codon aa)){print "$type";
#  foreach my $t (0.0001, 0.001, 0.01, 0.05, 0.10, 'total'){ my $c = $hash->{$t}->{$type} ? $hash->{$t}->{$type} : 0;
#  print "\t$c";
#  } print "\n";
# }' ;
#  fi;
#  done  >> $out
    # the below is ~ cut and pasted from Andrew's code

    my $file = shift;
    my $out = shift;
    my %hash;
    open IN, $file or die $!;
    while (<IN>){
	next if $_ =~ /^#/;
	my @F = split(/\t/);
	$hash{'total'}->{$F[2]}++;
  	foreach my $t (0.0001, 0.001, 0.01, 0.05, 0.1){
  	    next unless ($F[-1] < $t);
	    $hash{$t}->{$F[2]}++;
	}
    }
    close IN;
    cat_into($out,$file."\n");
    dump_n_n_total($out,\%hash);
}

sub dump_n_n_total{
    # print "#Type\t<0.0001\t<0.001\t<0.01\t<0.05\t<0.10\ttotal\n";
    # foreach my $type (qw(nuc codon aa)){
    # 	next unless (exists($hash->{$type}));
    # 	print "$type";
    # 	foreach my $t (0.0001, 0.001, 0.01, 0.05, 0.10, 'total'){
    # 	    my $c = $hash->{$type}->{$t} ? $hash->{$type}->{$t} : 0;
    # 	    print "\t$c";
    # 	} 
    # 	print "\n";
    # }

    my $file = shift;
    my $hash = shift;
    open OUT, ">>$file" or die $!;

    print OUT "#Type\t<0.0001\t<0.001\t<0.01\t<0.05\t<0.10\ttotal\n";
    foreach my $type (qw(nuc codon aa)){
  	print OUT "$type";
  	foreach my $t (0.0001, 0.001, 0.01, 0.05, 0.10, 'total'){ 
  	    my $c = $hash->{$t}->{$type} ? $hash->{$t}->{$type} : 0;
  	    print OUT "\t$c";
  	} 
  	print OUT "\n";
    }
    close OUT;
}

sub get_bam_file{
    my $n = shift;
    my $sample = shift;
    my $type = shift;
    return "Incorrect type" unless $type =~ /^all$|^maj$/;
    my $fn = $type eq 'all' ? join '.',('outputs/get_majority_block_bam.pl/'.$sample,'contigs.pid.btrim.fastq',$n,'bam') : 
      join '.',('outputs/get_majority_block_bam.pl/'.$sample,'contigs.pid.btrim',$n,'majority.bam');;
    return $fn if (-e $fn);
    warn "No such file $fn\n";
    return 0;
}

sub get_majority_block_filtering_stats{
    my $samples = shift;
    my $file = shift;
 # echo -e "Sample\tRegion\tAll\tMajority\tAllBAM\tMajorityBAM" > $out;
 # for n in 0 1 2;
 # do for s in 30_S1 31_S2 32_S3;
 # do allbam=$(ls ${s}.contigs.pid.btrim.fastq*${n}*.bam);
 # majbam=$(ls ${s}.contigs.pid.btrim.*${n}*.majority.bam);
 # echo -e "$s\t$n\t$(samtools view $allbam | wc)\t$(samtools view $majbam | wc)\t${allbam}\t${majbam}" >> $out;
 # done;
 # done 

    for my $sample (keys %$samples){
	for my $n (0, 1, 2){
	    my $all_bam_file = get_bam_file($n,$sample,'all');
	    next unless $all_bam_file;
	    my $maj_bam_file = get_bam_file($n,$sample,'maj');
	    my $phrase = join "\t",($sample, $n, samtools_view_wc($all_bam_file), samtools_view_wc($maj_bam_file), $all_bam_file, $maj_bam_file);
	    cat_into($file, $phrase."\n");
	}
    }
}

sub samtools_view_wc{
    my $file = shift;
    my $str = `$SAMTOOLS view $file | wc`;
    chomp $str;
    return $str;
}
    
sub cat_into{
    my $file = shift;
    my $phrase = shift;
    open IN, ">>$file" or die $!;
    print IN $phrase;
    close IN;
}

sub hash_file{
    my $file = shift;
    my $key_col = shift;
    my $val_col = shift;
    my %h;
    open FILE, $file or die "Dying; at hash_file, file name : $file\n";
    while (<FILE>){
	next if $_ =~ /^#|^\s?$/;
	chomp $_;
	my @cols = split /\t/, $_;
	$h{$cols[$key_col]} = $cols[$val_col];
    }
    close FILE or warn $!;
    return \%h;
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
