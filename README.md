# primer-id-progs
Programs needed to run the PrimerID pipeline

# Tutorial

### Test files and downloading output files.
The input files for the test are in the repository in the `primer-id-progs/data/test/input` directory.  All of the commands can be run in this directory.

### Preparing the environment
Before running the commands, make sure that you have a version of Perl on your PATH that has the necessary modules installed.  This has been tested with Perl 5.22.1.  Some of the required libraries are 

* BioPerl, including
  * Bio::Tools::Run::Alignment::Clustalw
  * Bio::DB::Sam
* Statistics::R
* Statistics::Distributions (included in the repo)
* Parallel::Loops
* File::Which
* Array::Utils
* File::Slurp
* IO::Zlib
* aomisc

To use the modules in the repository, modify your PERL5LIB to include the `primer-id-progs` directory.

The primerID workflow also depends on a number of external dependencies.  Binary executables for most of these dependencies have been provided.  In order to use these, modify your PATH variable to include the `primer-id-progs` directory and the directory to the `mafft` binary, e.g.,
```
export PATH=/path/to/primer-id-progs:/path/to/primer-id-progs/mafft-7.221:$PATH
```
For mafft to work, you also need to set this environment variable, pointing to the directory containing the mafft binary files (modify to suit your environment):
```
export MAFFT_BINARIES=/path/to/primer-id-progs/mafft-7.221/libexec/mafft
```

Some of the scripts use R as well, so be sure to have R and Rscript on your PATH.   This has been tested with R 3.2.3.  One required library is the 'network' package.


### Running the PrimerID workflow
The commands below are suggestions.  You are welcome to put them together in a shell script, or submit them on a cluster as you please.  In fact, it is recommended to run many of these on a cluster.  This small dataset has about 1000 primerID groups per amplicon (1000 primerID groups * 4 amplicons per dataset = 4000 primerID groups total approximately), so these jobs run relatively quickly compared to a full dataset, which may have 10s or 100s of thousands of primerID groups.  For full datasets, it is recommended to run jobs in parallel on a cluster (one amplicon per job).  Many of the scripts are parallelized to use multiple threads.  The usage statement may give a recommended number of threads to try, usually 8-12.

1. Trim first 4bp of R2 reads.  
  In our amplicons, we put 4bp of random nucleotide for the first 4 cycles for cluster generation on the MiSeq.  Your library design may be different, so modify as appropriate.
  ```
  for i in *_R2.fastq; do fastx_trimmer -i $i -o ${i/.fastq/}.trim.fastq -f 5 -Q 33; done
  ```
  This takes about 30 sec. to run.

2. Run fastqc (optional; not provided)
  ```
  for i in *R1.fastq *R2.trim.fastq; do fastqc --nogroup $i; done
  ```
  This takes about 1 min. 5 sec. to run.

3. Align with bwa mem and get coverage to see how many amplicons to expect and whether a gap or not. 
  ```
  bwa index HA_orf.fasta
  for i in 151_62_S1 141_64_S1 141_65_S2; do echo $i; bwa mem -t 8 HA_orf.fasta ${i}_R1.fastq ${i}_R2.trim.fastq | samtools view -bS - > ${i}.bam; java -Xmx2G -jar ~/Apps/primer-id-progs/picard.jar SortSam I=${i}.bam CREATE_INDEX=true O=${i}.sort.bam SO=coordinate TMP_DIR=/hpcdata/scratch/; graph_coverage.pl --bam_file ${i}.sort.bam --output_dir ./ ; done
  ```
  Takes 41 sec.

4. Make contigs from overlapping reads
  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do pandaseq -f ${i}_R1.fastq -r ${i}_R2.trim.fastq -F -T 8 -u ${i}_pandaseq_unaligned.fastq > ${i}.contigs.fastq; done
  ```
  Takes 5 sec.  Note that the binary of pandaseq provided will probably not work on your system.  You will need to compile a version for your system.  I am working on making a portable version and may include this at a later date.
  If your reads are not overlapping, you can use the `concatenate_reads.pl` script, which will reverse complement R2 and insert an appropriate number of Ns between the reads in a pair based on their alignment with the reference.

5. Extract primerID from sequences and place in the first line of the fastq record
  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do f=${i}.contigs.fastq; filter_fastq_by_primerid_length.pl --removepost --post CA --file_in $f --n 12; done
  ```
  Takes ~1 min. each. (3 min. 30 sec. total)

6. Split into amplicon regions and strip off primer sequences
  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do Btrim64 -p primers.txt -t ${i}.contigs.pid.fastq -o ${i}.contigs.pid.btrim.fastq -u 2 -v 2 -S -B -e 300; done
  ```
  Takes less than 1 sec.

7. Align with bwa mem, convert to BAM and get majority start/stop
  The purpose of this step is to remove off-target sequences.
  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do for fastq in ${i}.contigs.pid.btrim.fastq.*; do get_majority_block_bam.pl --ref HA_orf.fasta --fastq $fastq --output_dir ./; done; done
  ```
  Takes 45 sec.

8. Graph coverage, before and after cleaning up (optional)
  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do for sample in ${i}.contigs.pid.btrim.*.majority.bam ${i}.contigs.pid.btrim.fastq.*bam; do graph_coverage.pl --bam_file $sample --output_dir ./ ; done; done
  ```
  Takes 20 sec.

9. Merge the reads by primerid
  First run the command with `--plot_counts --plot_only` options to make the \*majority.group.counts.txt files, then run as normal.
  ``` 
  for i in 151_62_S1 141_64_S1 141_65_S2; do for sample in ${i}.contigs.pid.btrim.*.majority.bam; do a=${sample/.majority.bam/}; merge_primerid_read_groups.pl --plot_only --plot_counts $sample; done; done
  ```
  Takes 29 sec.
  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do for sample in ${i}.contigs.pid.btrim.*.majority.bam; do a=${sample/.majority.bam/}; groups=${a}.majority.group.counts.txt; cutoff=$(compute_cutoff.pl $(cat $groups | tail -n 1 | cut -f 1)); echo merge_primerid_read_groups.pl -m $cutoff --ambig 600 --min_freq 0.75 -p 8 $sample; done; done
  ```
  The \*majority.group.counts.txt files are required input to the `compute_cutoff.pl` script, in order to determine the minimum PrimerID group size.  There is a default minimum of 5 (this can be modified; suggested 3 to 5), and the script will calculate a statistical minimum based on the maximum group size.  In this case, you should see a minimum of 5 for all except 141_65_S2.contigs.pid.btrim.2.majority.bam will have a computed minimum of 6.

  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do for sample in ${i}.contigs.pid.btrim.*.majority.bam; do a=${sample/.majority.bam/}; groups=${a}.majority.group.counts.txt; cutoff=$(compute_cutoff.pl $(cat $groups | tail -n 1 | cut -f 1)); merge_primerid_read_groups.pl -m $cutoff --ambig 600 --min_freq 0.75 -p 8 $sample; done; done
  ```
  Takes 29 min. total, or about 2 min. each.

10. Convert merged reads to codons and amino acids and get frequency tables and clean read/peptide alignment files.
  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do for file in ${i}.contigs.pid.btrim.*.majority.cons.fasta; do a=${file/.majority.cons.fasta/}; convert_reads_to_amino_acid.pl --files $file --ref HA_orf.fasta --prefix ${file/.fasta/} -p 8; done; done 
  ```
  Takes 5 min. 21 sec.

11. Calculate linkage disequilibrium for all variants found above a certain frequency level 
  You can choose whatever threshold you would like.  It will first calculate the number of comparisons that will be performed, and then will give an estimate of the time it will take.  Right now the estimates are a little off, so it will take *longer* than the estimated time most likely.
  The command below is the suggested format for datasets where you have replicates. Using the `--group` option, the output will quantify all variants that are above the specified threshold in *either* of the replicates.  If this option is not used, then it will only look at variants that are above the threshold in that particular replicate so you won't have allele counts to compare between replicates for some variants. If your dataset doesn't have replicates, this is not as important.
  Note that lately I have noticed some instability in the Statistics::R package as implemented in this script, particularly with large datasets, so if you experience some problems this may be a known issue.  This script may get an update at some point in the near future.
  ```
  for v in 0.005; do for i in "141_64_S1 141_65_S2"; do set $i; for n in 0 1 2 3; do sample1=${1}.contigs.pid.btrim.${n}.majority.cons.variants.minfreq0.xls; sample2=${2}.contigs.pid.btrim.${n}.majority.cons.variants.minfreq0.xls; calculate_linkage_disequilibrium.pl $sample1,$sample2 ${sample1/.variants.minfreq0.xls/}.uniq.cleanreads.txt,${sample2/.variants.minfreq0.xls/}.uniq.cleanreads.txt ${sample1/.variants.minfreq0.xls/}.uniq.cleanpeptides.txt,${sample2/.variants.minfreq0.xls/}.uniq.cleanpeptides.txt --variant_threshold $v --group_id WT_Passage.${n} --label Passage1,Passage2; done; done; done
  ```
  Takes 42 sec.  Only a few comparisons are performed with this dataset, at this frequency threshold.
  ```
  for v in 0.005; do for i in 151_62_S1; do for n in 0 1 2 3; do sample=${i}.contigs.pid.btrim.${n}.majority.cons.variants.minfreq0.xls; calculate_linkage_disequilibrium.pl $sample ${sample/.variants.minfreq0.xls/}.uniq.cleanreads.txt ${sample/.variants.minfreq0.xls/}.uniq.cleanpeptides.txt --variant_threshold $v --prefix WT_Parent.${n}  --group_id WT_Parent.${n} --label Parent; done; done; done
  ```
  Takes 5 sec.  Only one comparison.

12. There are a few other steps that I haven't tested with this dataset at the moment, including `merge_overlapping_tally_regions.pl`, `combine_linkage_values.pl`, `compare_variant_frequencies.pl`, `primer_id_stats.pl`, and `graph_ambig_pos.R`.  I'll add some documentation for them at some point.  In the mean time, feel free to test them out. There is a usage statement for most that should describe how they are used (or you can inspect the script).

### Downloading the output files
You can download the output files expected for running these commands here:

https://drive.google.com/open?id=0B_uaeWUQ6aiJNFRkNkZWMl9qMzA

# Software terms of use
Use of this software prior to publication constitutes a collaboration. Please see Chris Brooke who has a copy of the terms of use.


# Notes for running on Locus.
v0.1.0 was tested to work on NIAID Locus.
Modules to load prior to running scripts:
* module load R/3.2.3-goolf-1.7.20-2016-Q1
* module load Perl/5.22.1-goolf-1.7.20-2016-Q1
* module load MAFFT/7.221-goolf-1.7.20-with-extensions





# Todo
* Description of scripts 
* Additional steps in the tutorial
