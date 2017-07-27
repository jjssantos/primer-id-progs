#-------------------------------------------------------------------------------
#----                                MODULE NAME                            ----
#-------------------------------------------------------------------------------
package aomisc;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use IO::File;
use IO::Zlib; 		#qw(:gzip_external 1);
use Carp;
use Data::Dumper;
use Cwd;
use File::Basename;
my @suffixes = qw(.bed .bed12 .bed6 .txt .fasta .fastq .fq .fa .fas .fna .png .pdf .gtf .gff .gff3 .sam .bam .xls .tab .ucsc .csv .names .qual .diff .gct .bg .bedGraph .g.vcf .vcf); 	#for fileparse.  Feel free to add more accepted extensions.		
our @SUFFIXES = _get_suffixes(\@suffixes);	#Needs to be "our" instead of "my" in order to export 
our @ISA = qw(Exporter);
our @EXPORT = qw(
		    pwd_for_hpc
	round
	average
	total
	hours_and_minutes
	median
	range
	simple_intersection
	stranded_intersection
	open_to_read
	open_to_write
	index
	get_files
	find_biggest_key
	find_key_with_biggest_value
	simple_hash
	elapsed
	get_file_num
	get_header
	column_header_lookup_hash
	bin_hash_of_values
	trim
	@SUFFIXES
	scalar_HoA
	check_for_Rscript
	stats
);


# Andrew J. Oler, PhD
# Created at
# Howard Hughes Medical Institute
# Dept of Oncological Sciences
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#
# With additions at
# Computational Biology Section
# Bioinformatics and Computational Biosciences Branch (BCBB)
# OCICB/OSMO/OD/NIAID/NIH
# Bethesda, MD 20892
#
# andrew.oler@gmail.com
# 
#This package is free software; you can redistribute it and/or modify
#it under the terms of the GPL (either version 1, or at your option,
#any later version) or the Artistic License 2.0.


#-------------------------------------------------------------------------------
#----------------------------------- CHANGELOG ---------------------------------
#-------------------------------------------------------------------------------
#111228
#Added scalar_HoA.
#111229
#Added total and changed average so it uses total
#120114
#Added step to ignore non-numerical keys in bin_hash_of_values. 
#Assign 0 as bottom of range before calculating bin width in bin_hash_of_values.
#120128
#Changed bin_hash_of_values to use Scalar::Util to check if value is numeric or not so that decimals are allowed. 
#120215
#Changed get_header to trim the headers of spaces on the ends.  
#120426
#Changed simple_hash to count every million lines read.
#121015
#Changed simple_hash to print out file name when keys not unique
#121019
#Added check_for_Rscript
#130114
# Allow passing of array either by reference or by its values.  Add this to each sub that takes an arrayref or array.  From here: http://ods.com.ua/win/eng/program/Perl5Unleashed/ch29.phtml
#	e.g., my $a = ref($_[0]) ? $_[0] : \@_;
#130222
# Modified find_key_with_biggest_value to return an arrayref if there are multiple keys with the same max value (optional boolean third argument).
# 2013-10-17
# Added if (defined($_)){ to trim subroutine.
# 2013-11-12
# Modified simple_hash subroutine to use //= for first and second column values, to allow zero to be a valid value.  
# 2017-04-05
# Removed subroutines to make venn diagram

# philip
sub pwd_for_hpc{
    my $prog_loc = Cwd::abs_path($0);	  # philip macmenamin
    my @a = split /\//,$prog_loc;	  # philip macmenamin
    return join '/', @a[0..$#a-1]; # philip macmenamin 
}
# To Do


#-------------------------------------------------------------------------------
#----------------------------------- FUNCTIONS ---------------------------------
#-------------------------------------------------------------------------------
sub round {
    my($number) = shift;
    return int($number + .5 * ($number <=> 0));
}
#-------------------------------------------------------------------------------
sub _get_suffixes {
	#Takes a basic list of suffixes and makes all combinations, including capitalized and adding zip extensions
	my $suffixes = shift;
	my @new_suffixes;
	my @zip = qw(.gz .zip .bz2);	#different type of zip extensions.
	foreach my $suffix (@{$suffixes}){
		$suffix = lc($suffix);
		push @new_suffixes, $suffix;
		push @new_suffixes, uc($suffix);
		foreach my $zip (@zip){
			push @new_suffixes, $suffix."$zip";
			push @new_suffixes, uc($suffix)."$zip";	
		}
	}
	return @new_suffixes;
}
#-------------------------------------------------------------------------------
sub average {
	my $array = ref($_[0]) ? $_[0] : \@_;		# Takes an array or arrayref Was my $array = $_[0];
	my $num_elements = @{$array};
	my $total = total($array);
	my $average = 0;
	if ($num_elements > 0 ){
		$average = $total / $num_elements;
	}
	return $average;
}
#-------------------------------------------------------------------------------
sub total {
	my $array = ref($_[0]) ? $_[0] : \@_;		# Takes an array or arrayref Was my $array = shift;
	my $total = 0;
	foreach my $element (@{$array}){
		$total = $total + $element;
	}
	return $total;
}
#-------------------------------------------------------------------------------
sub hours_and_minutes{
	my $time = $_[0];
	my $total_minutes = $time / 60;
	my $hours = $total_minutes / 60;
	my $hours_int = int($hours);
	my $remaining_minutes_fraction = $hours - $hours_int;
	my $remaining_minutes = $remaining_minutes_fraction * 60;
	my $remaining_minutes_int = int($remaining_minutes);
	my $hours_and_minutes_string;
	if ($hours_int > 0){
		$hours_and_minutes_string = "$hours_int hours $remaining_minutes_int minutes";
	}	
	elsif ($remaining_minutes_int > 1) {
		$hours_and_minutes_string = "$remaining_minutes_int minutes";
	}
	else {
		$hours_and_minutes_string = round($time) . " seconds";
	}
	return $hours_and_minutes_string;	
}
#-------------------------------------------------------------------------------
sub median {
	#this could be in a module (it is!)
	my $array = ref($_[0]) ? $_[0] : \@_;		# Takes an array or arrayref Was my $array = $_[0];
	my @array = sort { $a <=> $b } @{$array}; 
	my $count = scalar(@array); 
	my $median;
	if ($count % 2) { 	# % returns the remainder.  If there is a remainder, then it is odd-numbered.
		$median = $array[int($count/2)]; 
	} else { 	#even number of items in the array
		$median = ($array[int($count/2)] + $array[int($count/2 - 1)]) / 2; 
	}
	return $median;
}
#-------------------------------------------------------------------------------
sub range {
	#Can get range (default), min, or max value in an array.  
	my $array = ref($_[0]) ? $_[0] : \@_;		# Takes an array or arrayref Was my $array = $_[0];
	shift;
	my $max_min = shift;	#optional.  String, either 'min' or 'max'.  If 'max', only return maximum.  If 'min', only return min.  If not specified, return an array min, max.
	my @array = sort { $a <=> $b } @{$array}; 
	my $first = $array[0];
	my $last = $array[-1];
	if ($max_min){
		if ($max_min eq 'max'){
			return $last;
		}
		elsif ($max_min eq 'min'){
			return $first;
		}
		else {
			print STDERR "$max_min not matching required value, i.e., either max or min.  Will return default, both min and max.\n";
			return $first, $last;
		}
	}
	else {
		return $first, $last;
	}
}
#-------------------------------------------------------------------------------
sub scalar_HoA {
	my $HoA = shift;
	my $count = 0;
	foreach my $key (keys %{$HoA}){
		foreach my $element (@{$HoA->{$key}}){
			$count++;
		}
	}
	return $count;
}
#-------------------------------------------------------------------------------
sub simple_intersection {
	#Conditions of overlap: stop 1 [+gap]> start2 && start 1 [-gap]< stop 2
	my ($start1,$stop1,$start2,$stop2,$gap) = @_;
	if((($stop1+$gap)>=$start2)&&(($start1-$gap)<=$stop2)){
		return 1;
	}
	else{
		return 0;
	}	
}
#-------------------------------------------------------------------------------
sub stranded_intersection {
	#Conditions of overlap: stop 1 [+gap]> start2 && start 1 [-gap]< stop 2
	my ($start1,$stop1,$strand1,$start2,$stop2,$strand2,$gap) = @_;
	if((($stop1+$gap)>=$start2)&&(($start1-$gap)<=$stop2)&&($strand1 eq $strand2)){
		return 1;
	}
	else {
		return 0;
	}
}
#-------------------------------------------------------------------------------
sub index {
	#This seems to work.		I could try to shorten the index so it won't have to search as many.
	my $length = length($_[0])-4;
	my $index = 0;
	if ($length >= 0){
		$index = substr($_[0],0,$length) || 0;	#this should get everything except the last 4 digits.	Does this work for the size of region to know how many extra indices to search?
	}
	return $index;
}
#-------------------------------------------------------------------------------
sub open_to_read {
	#Similar to open_to_read_fh in tim_file_helper.pm
	#Can open a regular file or one ending in .gz.  Only opens the file for READING, not writing or read/write.
	#my $readfh = open_to_read($file);
	my $path_to_open = shift;
	if (-e $path_to_open){
		my $readfh;
		if ($path_to_open =~ m/\.gz$/){	
			$readfh = IO::Zlib->new;
		}
		else {
			$readfh = IO::File->new;
		}
		
#		print Dumper($readfh);
		
		if ($readfh->open($path_to_open, "r") ) {
			return $readfh;
		}
		else {
			carp "Can\'t open file '$path_to_open': " . $readfh->error . "\n";		#Not sure if $readfh->error is proper call of method.  This is what Tim Parnell had done.
			return;	
		}
	}
	else {
		warn "aomisc::open_to_read: file $path_to_open does not exist\n";
		exit 1;		
	}
#	return $readfh;
}
#-------------------------------------------------------------------------------
sub open_to_write {
	#Can open a regular file or one ending in .gz.  Only opens the file for WRITING, not reading or read/write.
	#First argument is the path to the file to write to, second argument is whether to compress the output file (1 for yes, 0 for no; default = no).
	#By default, it creates a new file, which deletes any existing file with that name.  Alternatively, you can append with '>>' in the third argument
	#(I may want to add a function that allows you to have the option to append instead of create, like in open_to_write_fh in tim_file_helper.  I would use the sub in his module, but it has features that are deprecated.)
	#111117: Changed to increment $gzip if .gz is found at the end.  Then in checking for .gz to add, only check if $gzip.
	#use PerlIO::gzip;		#does this work to have it here instead?  Looks like it works okay.
	my ($path_to_open,$gzip,$write_symbol,$quiet) = @_;
	$write_symbol ||= '>';		#default is to not append
	if ($path_to_open =~ m/\.gz$/){	#if the extension is .gz, then we should zip it.
		$gzip++;
	}
	my $writefh;
	my $mode;
	if ($write_symbol =~ m/^>>$/){	#append
		$mode = 'a';
	}
	else {	#write a new file; this is the default
		$mode = 'w';
	}
	if($gzip){	#binary write mode, make it ab or wb.  NOTE that ab mode is not legal with external gzip in IO::Zlib
		$mode .= 'b';
	}
	
	#Add .gz if necessary
	if ($gzip){	
		unless ($path_to_open =~ m/\.gz$/){	#add .gz if it's not there already.
			$path_to_open .= '.gz';
		}
	}
	
	#create an anonymous filehandle
	if ($gzip){
		$writefh = IO::Zlib->new;
	}
	else {
		$writefh = IO::File->new;
	}
	
	if ($writefh->open($path_to_open, $mode) ){
		unless ( ($write_symbol =~ m/>>/)||($quiet) ){	warn "out: $path_to_open\n";	}		#print to STDERR when creating a new file.
		return $writefh;
	}
	else {
		carp "Can\'t open file '$path_to_open' in mode $mode: " . $writefh->error . "\n";		#Not sure if $writefh->error is proper call of method.  This is what Tim Parnell had done.
		return;
	}

}
#-------------------------------------------------------------------------------
sub get_files{
	#Takes either a directory of files or a comma-delimited list of files
	#First argument, list of files or a directory; second argument, extension to use in grep to grab if first argument is a directory. **extension can be a regular expression.** e.g., 'bed' or '(bed|bed.gz)'
	my @files;
	my $ext = '[^\.]+';		#default extension: [anything after the last period, see grep below]
	if ($_[1]){
		$ext = $_[1];
	}
	if ($_[0] =~ m/,/){	#comma-delimited list of files
		@files = split(/,/, $_[0]);
	}
	else {				#should be a directory, or a single file, or string 'stdin'
		if (-d $_[0]){		#directory exists
			opendir(DIR, $_[0]);
			@files = grep(/\.$ext$/,readdir(DIR));	#Find the extension regex at the end of the line 
			for (my $i = 0; $i<@files; $i++){
				$files[$i] = $_[0]."/".$files[$i];
			}
			closedir(DIR);
		}
		elsif(-e $_[0]){
			$files[0] = $_[0];
		}
		elsif($_[0] =~ m/stdin/i){
			$files[0] = 'stdin';
		}
		else {
			warn "File/Directory $_[0] does not exist.\n";
			exit 1;
		}
	}	
	#warn join "\n", @files, "\n";
	return @files;
}
#-------------------------------------------------------------------------------
sub find_biggest_key {
	#could be replaced by "range" function, right?
	my $hash = shift;
	my @sizes = sort { $b <=> $a } (keys (%{$hash}));	#sort (keys (%{$data}))
	my $max = $sizes[0];
	return $max;	
}
#-------------------------------------------------------------------------------
sub find_key_with_biggest_value {
	# Finds the key-value pair with the highest value, then reports the key, the value, or both
	# What if two keys have the same highest value?  Specify $all_highest to find them all
	my $hash = shift;
	my $with_value = shift || 0;	#optional. If not specified, it will return just the key. If value is 1, it will return the value.  If value is 2, it will return key, value. 
	my $all_highest = shift || 0;	#optional.  If not specified, it will return just one key-value pair (either key, value, or both, as specified by $with_value).  If value is 1, it will return an arrayref of the results.
	my @keys = sort { $hash->{$b} <=> $hash->{$a} } (keys (%{$hash}));
	my $max_key = shift(@keys);
	if ($all_highest == 1){
		if ($with_value == 1){		#values
			my @max = $hash->{$max_key};
			while($keys[0] && ($hash->{$keys[0]} == $hash->{$max_key}) ){
				my $key = shift @keys;
				push @max, $hash->{$key};
			}
			return \@max;
		}
		elsif($with_value == 2){	#keys and values
			my $max->[0] = [$max_key, $hash->{$max_key}];
			while($keys[0] && ($hash->{$keys[0]} == $hash->{$max_key}) ){
				my $key = shift (@keys);
				push @$max, [$key, $hash->{$key}];
			}
			return $max;
		}
		else {						#keys
			my @max = ("$max_key");
			while($keys[0] && ($hash->{$keys[0]} == $hash->{$max_key}) ){
				my $key = shift (@keys);
				push @max, $key;
			}
				return \@max;
		}
	}
	else {
		if ($with_value == 1){
			return $hash->{$max_key};
		}
		elsif($with_value == 2){
			unless ($max_key && exists($hash->{$max_key})){
				carp "ERROR: Max key or value in hash not defined.  max_key: $max_key\n" . Dumper($hash) . "\n";
			}
			return $max_key.", ".$hash->{$max_key};
		}
		else {
			return $max_key;
		}
	}
}
#-------------------------------------------------------------------------------
sub simple_hash {
	#Arguments: $fh, $first_column, $second_column, [$unique, $skip, $separator, $quiet]
	#e.g., my $fh = open_to_read($file); my $hash = simple_hash($fh,0,4); close ($fh);
	#e.g., my $fh = open_to_read($file); my $hash = simple_hash($fh,0,'count'); close($fh);
	#New 110607 (counts by default, default first_column = 0, can take file instead of fh):
	#Arguments: $file or $fh, [$first_column, $second_column, $unique, $skip, $separator, $quiet]
	#e.g., my $hash = simple_hash($file);
	#returns hashref
	my $hash;
	my ($fh,$first,$second,$make_unique,$skip,$separator,$quiet) = @_;
	#	$fh	#a previously opened filehandle
	$first //= 0;	#first column, e.g., 0.  Will be key, therefore should be unique. (zero is a valid value)
	$second //= 'count';	#second column, will be value.  If 'count', then it will count instead of using another column as value.   (zero is a valid value)
	$make_unique ||= 0;	#optional, 1 if yes, 0 if no.  Default = 0;  no will assume they are unique.  Prints an error if at least one is found that is not unique...
	$skip ||= '^#';			#optional, regexp for line to skip.  If nothing, then it will skip commented lines.  To force retention of commented lines, send 0 to this argument.
	$separator ||= '\t';	#optional, regexp for separator, e.g., '\t'.  Default is \t
	$quiet ||= 0;	#optional.  default is to be verbose, mentioning when keys are not unique, when skipping keys, etc.  If 1, then will not give these warnings.
	
	#If passed a file instead of a filehandle, open an anonymous filehandle
	my $opened_in_sub = 0;
	my $file = "";
	if (defined(IO::Handle::fileno($fh))){	#if file handle, this will be defined as a number.
		#print "file handle: " . IO::Handle::fileno($fh) . "\n";
	}
	else {	#then it is just a file path
		#print "no file in fh.\n";
		$file = $fh;
		$fh = open_to_read($fh);
		$opened_in_sub++;
	}
	my $warned_unique = 0;
	my $warned_count = 0;
	my $lines_done = 0;
	my $start_time = time;
	while (<$fh>){
		chomp;
		next if ($_ =~ m/$skip/i);
		my @line = split(/$separator/, $_);
		next unless(@line);	#to skip empty lines
		if ($second =~ m/count/i){
			$hash->{$line[$first]}++;
		}
		else {
			#Check if there really is a value in the second column
			if (defined($line[$second])){	#use "defined" in case the value actually is zero, which is a valid value
				#There is a value in second column, so it will become the value in the hash
				if (exists($hash->{$line[$first]})){	#to stop hash entries from being erased 
					if ($make_unique){
						UNIQUE: for (my $i=1;$i<10000;$i++){		#Make them unique
							if (exists($hash->{$line[$first].'.'.$i})){
								#Then don't use this name.
							}
							else {	#This name is available; take it.
								$hash->{$line[$first].'.'.$i}=$line[$second];	
								last UNIQUE;
							}
							if ($i == 10000){
								warn "aomisc::simple_hash warning: $line[$first] already exists in hash > 10000 times, skipping...\n";
							}
						}
					}
					else {	#just warn
						unless ($warned_unique){	
							unless ($quiet){
								warn "aomisc::simple_hash warning: Simple hash keys are not unique, e.g., $line[$first].\n";
								if ($opened_in_sub){	warn "\tfile: $file\n";	}	
								$warned_unique++;	
							}
						}
					}
				}
				else {	#not in hash, so use it.
					$hash->{$line[$first]}=$line[$second];
				}
			}
			else {
				#There really isn't a value in the second column, so default to counting
				unless ($warned_count){
					unless($quiet){
						warn "aomisc::simple_hash warning: no value in column of index $second, will return a hash of count values instead\n";
						$warned_count++;
					}
				}
				$hash->{$line[$first]}++;
			}
		}
		unless($quiet){
			$lines_done++;
			if (($lines_done % 1000000) == 0){	
				print STDERR "Done processing $lines_done lines.";
				&elapsed($start_time, ' Elapsed', 1);
			}			
		}
	}
	if ($opened_in_sub){	#then close in sub
		close ($fh);
	}
	return $hash;
	#remember to close your $fh in the script (unless opened in subroutine)
}
#-------------------------------------------------------------------------------
sub elapsed {
	#If using this module, you should have this at the beginning of the script:
	#my $start_time = time;
	#e.g., 	&elapsed($start_time, ' Elapsed', $verbose);
	my $start_time = shift;		
	my $pretext = shift;	#e.g., 'Elapsed',or 'Total'
	my $verbose = shift;
	my $elapsed_time;
	if($verbose){	
		$elapsed_time = (time - $start_time);
		$elapsed_time .= ' seconds';
	}
	else{
		$elapsed_time = &hours_and_minutes(time - $start_time);
	}	
	warn "$pretext time: $elapsed_time\n";
	
}
#-------------------------------------------------------------------------------
sub get_file_num {
	#Used for determining a good number to append to a file.
	#e.g., my $outfile = get_file_num($save_dir,$file_text, 'fa', $verbose);
	#$save_dir is the directory in which to look for an existing file
	#$file_text is the name of the file.  $ext is the extension that will be used.  
	#$verbose is 1 if verbose output is desired, 0 or nothing if not desired.
	#returns the full path to the filename with the appropriate number.
	my ($save_dir,$text,$ext,$verbose) = @_;
	my $good = 0;
	for (my $num = 1; $good<1; $num++){
		my $filename = $save_dir."/$text"."_".$num.".".$ext;
		if (-e $filename){
			if ($verbose){	print STDERR "file exists $filename ... trying another\n";	}
		}
		else {
			$good++;
			#return $num;
			return $filename;
		}
	}
}
#-------------------------------------------------------------------------------
sub bin_hash_of_values{
	#Generic 10 bins within range.  Takes a hashref of key-value pairs to be binned, also, second argument is the number of bins (default 10).  Third argument is whether verbose output is desired (1) or ignored (0).
	#This could be in a module
	#This comes from BedAlignStats.pl and was adapted from bin_values. 
	#example of calling:
		#		print STDERR "Alignment scores, binned:\n";
		#		bin_hash_of_values($alignment_scores, $bins, $verbose, $bin_width, $min, $max);
	#To do:  
	#1. Add features: 1) bin width, 2) max range (all values above this will be included in the same bin)
	#2. Change input options to have a hash of key/values instead, maybe, to simplify things.  
	#110827
	#Attempting to allow $quotient to be passed to the sub (i.e., bin width).
	#Logic: if $bin_num and bin width are chosen, they could conflict.  E.g., if range is 1-100, 
	#bin_num is 10 and bin_width is 20, you cannot satisfy both.  To satisfy bin_num 10, you need 
	#bin_width 10; to satisfy bin_width 20, you need bin_num 5.  $quotient takes precedence.  If 
	#$quotient is defined, $bin_num is not used.  
	#Also, allowing $min to be passed to sub.  Default mode is to use the lower end of the range.  
	#Any values below $min will not be included.  
	#Also, allowing $max to be passed to sub.  Default mode is to use the higher end of the range.  
	#If max is defined, values above that will still be binned, but all in one bin.  
	use Scalar::Util qw(looks_like_number);
	
	my $scores = $_[0];		#hash where keys are numbers and values are counts for the number of time that value is present in the data. (e.g., quality scores of reads, if there are 30 reads with quality score of 60, then the key '60' will have value 30.)
	my $bin_num = $_[1] || 10;
	my $verbose = $_[2] || 0;
	my $quotient = $_[3] || 0; 
	my $min = defined $_[4] ? $_[4] : -2**31;	#zero is a valid value for this variable, so use // operator.  default is -2**31, which will be unused.   $_[4] // -1;
	my $max = defined $_[5] ? $_[5] : 2**31;
	if ($verbose){	print STDERR "finding range...\n";	}
	my @keys = keys(%{$scores});
	my $skipped_keys = 0; my @new_keys; foreach (@keys){ if (looks_like_number($_)){	push @new_keys, $_; }	else {$skipped_keys++;}	}		#Added this @new_keys array to get rid of non-numerical keys.  Sometimes we get these in HashByColumn.pl using option --count.  (Alternatively, I could fix the hash so that it is all numerical before sending to this subroutine...)  120114.  Fixed this 120128 to use Scalar::Util qw(looks_like_number); instead of 'if ($_ =~ m/^\d+$/)' to check if numeric.
	@keys = @new_keys;		#Also added 120114
	if ($skipped_keys){	print STDERR "skipped $skipped_keys non-numerical keys\n";}		#Also added 120114.
	my @range = aomisc::range(\@keys);
	my %bins;
#	print Dumper(\@range);
#	print Dumper(\@keys);
	if ($min>-2**31){		#including if it == 0.  
		#Then assign this to the bottom range.  #changed 120114 so that this was done before calculating $quotient.  Before, this was done after defining quotient, which caused the number of bins to be incorrect.
		$range[0] = $min;
	}
	if ($max<2**31){
		#Then assign this to the top range.
		$range[-1] = $max;
	}
	
	unless ($quotient){	
		if ($verbose){	print STDERR "range: ", join " ", @range, "\nfinding bin size...\n";	}
		$quotient = aomisc::round(($range[1]-$range[0])/($bin_num));		
	}
	unless ($quotient){	$quotient = 1;}	#so that it doesn't become zero to make an infinite loop
	
	my @bins;
	my $bin = 0;
	if ($verbose){	print STDERR "bin size: $quotient\nmaking bins...\n";	}
	for (my $i=1; $bin<$range[1]; $i++){
		$bin = $i*$quotient-1+$range[0];
		push @bins, $bin;
	}
	pop @bins; 
	push @bins, $range[1];	#so the histogram doesn't go beyond the actual range of values.  (we take care of the beginning value below in printing: $lastbin = $range[0]; starting with range[0])
	if ($verbose){	print STDERR "bins: ", join " ", @bins, "\n";	}
	foreach $bin (@bins){
		$bins{$bin} = 0;
	}
	SCORE: foreach my $score (sort {$a <=> $b} keys %{$scores}){
		my $last_bin = -1;
		
		BIN: foreach $bin (@bins){	#Find the right bin for the score
			if (($score > $last_bin)&&($score<=$bin)){
				$bins{$bin}+= $scores->{$score};
				next SCORE;
			}
			else {
				$last_bin = $bin;
				next BIN;
			}
		}
	}
	#print Dumper(\%bins);
	my $lastbin = $range[0];
	foreach $bin(sort {$a <=> $b} keys %bins){
		if ($lastbin eq $bin){
			print STDERR "$bin\t$bins{$bin}\n";
		}
		else {
			print STDERR "$lastbin\t$bin\t$bins{$bin}\n";
		}
		$lastbin = $bin+1;
	}
}
#-------------------------------------------------------------------------------
sub get_header {
	#Retrieves the headers of a file (first line) and puts into an array, returns 
	#the array.  Column headers should be tab-delimited.
	my $file = shift;
	my $readfh = open_to_read($file);
	my @array;
	LINE: while(<$readfh>){
		chomp;
		my @line = split("\t", $_);
		@array = trim(@line);		#Added 120215.  was @array = @line.
		last LINE;
	}	
	close($readfh);
	return @array;
}
#-------------------------------------------------------------------------------
sub column_header_lookup_hash {
	# Takes a file or an arrayref as argument
	# Returns a hashref for you to be able to lookup the index of a column header 
	# For example, headers (chr, start, end), hash would look like this:
	# %hash = (
	#	'chr' 	=> 0,
	#	'start' => 1,
	#	'end' 	=> 2,
	# );
	# Takes a file, returns a hashref
	# Each column header should be unique
	# Column headers should be tab-delimited
	
	#my $file = shift;
	#my @header = get_header($file);
	my $input = shift;	# Either a file or an array ref.  
	my @header = 	ref($input) ? @$input : get_header($input);
	
	my $hash;
	for (my $i = 0; $i < @header; $i++){
		my $col_header = $header[$i];
		$hash->{$col_header} = $i;		# Maybe add something here to check if that already exists so not to overwrite an earlier stored value.  
	}	
	return $hash;	
}
#-------------------------------------------------------------------------------
sub trim {
	#Modified from Perl Cookbook, recipe 1.19
	#Trims whitespace before and after text
	#Can take an array or a single value
	#e.g., @line = trim(@line);
	#e.g., $number = trim($number);
	my @out = @_;
	for (@out){
		if (defined($_)){
			s/^\s+//;	#trim left
			s/\s+$//;	#trim right
		}
	}
	return @out == 1
			? $out[0]	#only one to return
			: @out;		#or many
}
#-------------------------------------------------------------------------------
sub check_for_Rscript {
	#Check for Rscript on the PATH.
	my $output = `Rscript 2>&1`;	#captures stdout and stderr.  If present, value will start: Usage: /path/to/Rscript [--options] [-e expr] file [args]".  If not present, value will be undef.
	
	unless ($output){
		print STDERR "No R installation detected.  Please install R or add to PATH environment variable to run this program.\n"; 
		return 1;
	}
	return 0;	#successfully found Rscript
}
#-------------------------------------------------------------------------------
sub stats {
	# Modified from http://ods.com.ua/win/eng/program/Perl5Unleashed/ch29.phtml
	# E.g.,  my ($ave,$max,$min,$std,$count) = stats(\@array);
	# Allow passing of array either by reference or by its values.

	my $a = ref($_[0]) ? $_[0] : \@_;
	my $count = $#{$a} + 1;
	
	# Bail out in case of erroneous data.
	return(-1,-1,-1,-1) if ($count < 2);
#	print "$count items \n";

	my $i;
	# Initialize local variables. The assignment to 0 
	# is unnecessary for all scalars except $max and $min
	# since Perl will initialize them to zero for you.
	#
	my $min = $$a[0];
	my $max = $$a[0];
	my $sum = 0;
	my $sum2 = 0;
	my $ave = 0;
	my $std = 0;

	# Get the required statistics 
	for $i (@$a) {
		$sum += $i;
		$max = $i if ($max < $i);
		$min = $i if ($min > $i);
	}
	$ave = $sum/$count;
	for $i (@$a) {		# Added this new for st dev calculation.  The formula at the website listed was incorrect
		$sum2 += ($i - $ave)**2;
	}
	$std = sqrt($sum2/($count - 1));			# Was $std = (($sum2 - $sum * $ave)/($count - 1));

	# Return the list of values back from function.
	#
	return ($ave,$max,$min,$std,$count); 
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
1;


