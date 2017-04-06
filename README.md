# primer-id-progs
Programs needed to run the PrimerID pipeline

# Tutorial

## Test files and downloading output files.
The input files for the test are in the repository in the `primer-id-progs/data/test/input` directory.  All of the commands can be run in this directory.
You can download the output files expected for running these commands here:

https://drive.google.com/open?id=0B_uaeWUQ6aiJNFRkNkZWMl9qMzA

## Preparing the environment
Before running the commands, make sure that you have a version of Perl on your PATH that has the necessary modules installed.  This has been tested with Perl 5.22.1 in a Linux environment (RHEL).  Some of the required libraries are 

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
* aomisc (part of this repository)
* primerid (part of this repository)

To use the modules in the repository, modify your PERL5LIB to include the `primer-id-progs` directory, e.g.,
```
export PERL5LIB=/path/to/primer-id-progs:$PERL5LIB
```

The primerID workflow also depends on a number of external dependencies.  Binary executables for most of these dependencies have been provided.  In order to use these, modify your PATH variable to include the `primer-id-progs` directory and the directory to the `mafft` and `pandaseq` binaries, e.g.,
```
export PATH=/path/to/primer-id-progs:/path/to/primer-id-progs/mafft-7.221:/path/to/primer-id-progs/pandaseq-2.9:$PATH
```
For mafft to work, you also need to set this environment variable, pointing to the directory containing the mafft binary files (modify to suit your environment):
```
export MAFFT_BINARIES=/path/to/primer-id-progs/mafft-7.221/libexec/mafft
```

Some of the scripts use R as well, so be sure to have R and Rscript on your PATH.   This has been tested with R 3.2.3.  One required library is the 'network' package.


## Running the PrimerID workflow
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
  tempdir=temp	# Modify this to suit your environment
  mkdir -p $tempdir
  for i in 151_62_S1 141_64_S1 141_65_S2; do echo $i; bwa mem -t 8 HA_orf.fasta ${i}_R1.fastq ${i}_R2.trim.fastq | samtools view -bS - > ${i}.bam; java -Xmx2G -jar ~/Apps/primer-id-progs/picard.jar SortSam I=${i}.bam CREATE_INDEX=true O=${i}.sort.bam SO=coordinate TMP_DIR=$tempdir; graph_coverage.pl --bam_file ${i}.sort.bam --output_dir ./ ; done
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
  First run the command with `--plot_counts --plot_only` options to make the \*majority.group.counts.txt files.
  ``` 
  for i in 151_62_S1 141_64_S1 141_65_S2; do for sample in ${i}.contigs.pid.btrim.*.majority.bam; do a=${sample/.majority.bam/}; merge_primerid_read_groups.pl --plot_only --plot_counts $sample; done; done
  ```
  Takes 29 sec.
  
  The \*majority.group.counts.txt files are required input to the `compute_cutoff.pl` script, in order to determine the minimum PrimerID group size.  There is a default minimum of 5 (this can be modified; suggested 3 to 5), and the script will calculate a statistical minimum based on the maximum group size.  In this case, you should see a minimum of 5 for all except 141_65_S2.contigs.pid.btrim.2.majority.bam will have a computed minimum of 6.  

  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do for sample in ${i}.contigs.pid.btrim.*.majority.bam; do a=${sample/.majority.bam/}; groups=${a}.majority.group.counts.txt; cutoff=$(compute_cutoff.pl $(cat $groups | tail -n 1 | cut -f 1)); merge_primerid_read_groups.pl -m $cutoff --ambig 600 --min_freq 0.75 -p 8 $sample; done; done
  ```
  Takes 29 min. total, or about 2 min. each.

10. Convert merged reads to codons and amino acids and get frequency tables and clean read/peptide alignment files.
  ```
  for i in 151_62_S1 141_64_S1 141_65_S2; do for file in ${i}.contigs.pid.btrim.*.majority.cons.fasta; do convert_reads_to_amino_acid.pl --files $file --ref HA_orf.fasta --prefix ${file/.fasta/} -p 8; done; done 
  ```
  Takes 5 min. 21 sec.
  
  If you have a large number of unique reads in your dataset (i.e., > 20,000) and you only need frequency tables, it is recommended to follow one of the alternative procedures described in the section below titled "Alternative procedures for convert_reads step".

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

12. There are a few other steps that I haven't tested with this dataset at the moment, including `merge_tally.pl`, `combine_linkage_values.pl`, `compare_variant_frequencies.pl`, `primerid_stats.pl`, and `graph_ambig_pos.R`.  I'll add some documentation for them at some point.  In the mean time, feel free to test them out. There is a usage statement for most that should describe how they are used (or you can inspect the script).

## Alternative procedures for convert_reads step for fasta file with > 20,000 reads
  Fasta files with a large number of reads can take a significant amount of RAM and time to complete the convert_reads step.  For example, it takes about 80 hours and 100GB RAM to process ~120,000 unique reads using 8 threads.  If you have more than about 20,000 unique reads in your dataset, it is recommended to split the input file, processing each file with convert_reads and then merging the output with `merge_tally.pl`.  This process can reduce the processing time to hours/minutes.  For example, a file with ~1.5 million unique reads split into files about 10,000 reads each will take about 30 min. to process if all ~150 convert_reads jobs are running in parallel. Tally files for nuc, codon, and aa are merged separately.  If desired, you can merge tally files for different amplicons of the same chromosome/segment using `merge_tally.pl` as well.  Simply put all of the files to merge in the same directory and pass the directory to the --input parameter for `merge_tally.pl`. The Phylip output (required for the next step, calculating linkage disequilibrium) is only available using the regular workflow.
  To try this procedure with the tutorial files, start with the \*.cons.fasta files after the merge_primerid_read_groups step.  First split the files into ~100 reads each (in reality, you would only need to do this procedure with much larger datasets and you would split them into files of about 10,000 reads each).  We'll use the Unix `split` command to do this. To keep one merged output file for each amplicon, follow procedure 1 below. To merge frequency tables for all amplicons of a segment/chromosome into one file, follow procedure 2 below.

### Procedure 1 (keep amplicons separate): 
  Create a separate directory for each amplicon and split the reads for the amplicon into that directory.  Use -l 200 to get 100 reads per file.
```
for i in *majority.cons.fasta; do mkdir ${i/.majority.cons.fasta/}_split; done
for i in *majority.cons.fasta; do split -a 3 -d --additional-suffix=".fasta" -l 200 $i ${i/.majority.cons.fasta/}_split/${i/majority.cons.fasta/}; done
```
Now run convert_reads for each split fasta file.
```
for dir in *_split; do for file in $dir/*fasta; do convert_reads_to_amino_acid.pl --files $file --ref HA_orf.fasta --prefix ${file/.fasta/} -p 8; done; done 
```
Takes 7 min.
Now merge the reads with `merge_tally.pl`
```
for dir in *_split; do sample=${dir/_split/}; sample=${sample/.contigs.pid.btrim/}; merge_tally.pl -i $dir --prefix Merged --sample $sample; done
```
Takes 8 sec.

### Procedure 2 (merge data for amplicons belonging to the same segment):
Create separate directory for each sample (and segment if you are assessing multiple segments) and split the reads into that directory.  Use -l 200 to get 100 reads per file.
```
for i in 151_62_S1 141_64_S1 141_65_S2; do mkdir ${i}_split; done
for i in 151_62_S1 141_64_S1 141_65_S2; do for file in ${i}*majority.cons.fasta; do split -a 3 -d --additional-suffix=".fasta" -l 200 $file ${i}_split/${file/majority.cons.fasta/}; done; done
```
Now run convert_reads for each split fasta file.
```
for dir in *_split; do for file in $dir/*fasta; do convert_reads_to_amino_acid.pl --files $file --ref HA_orf.fasta --prefix ${file/.fasta/} -p 8; done; done
```
Takes 7 min.
Now merge the reads for each sample with `merge_tally.pl`
```
for dir in *_split; do sample=${dir/_split/}; merge_tally.pl -i $dir --prefix Merged --sample $sample; done
```
Takes 6 sec.


# Software terms of use
Use of this software prior to publication constitutes a collaboration. Please see Chris Brooke who has a copy of the terms of use.  Use and distribution of the external dependencies is subject to their respective licenses (see below).

## Public Domain license

The software in this repository authored by officers or employees of the National Institutes of Health (NIH) is free and unencumbered software released into the public domain. This does not include third-party software also included in the repository; see below for details on third-party software licenses.
</br>*Be kind, and provide attribution when you use this code.*

United States government creative works, including writing, images, and computer code, are usually prepared by officers or employees of the United States government as part of their official duties. A government work is generally not subject to copyright in the United States and there is generally no copyright restriction on reproduction, derivative works, distribution, performance, or display of a government work. Unless the work falls under an exception, anyone may, without restriction under U.S. copyright laws:

* Reproduce the work in print or digital form
* Create derivative works
* Perform the work publicly
* Display the work
* Distribute copies or digitally transfer the work to the public by sale or other transfer of ownership, or by rental, lease, or lending

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Learn more about how copyright applies to U.S. government works at [USA.gov](https://www.usa.gov/government-works)


## Attribution for Third-Party Software

This repository includes open-source third-party applications in object (binary) form. See corresponding links below for license and copyright notices.  Please cite applicable references below if using this software in your research (in addition to citing this repository itself).

#### [Samtools](http://www.htslib.org/)
_Used by `merge_primerid_read_groups.pl` and `primerid_stats.pl` and included for convenience to manipulate SAM/BAM files_

Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. PMID: 19505943.

#### [BWA](https://github.com/lh3/bwa)
_Used by `bwa_index_ref.pl`, `concatenate_fastq.pl`, `get_majority_block_bam.pl` and `merge_primerid_read_groups.pl` and also included for convenience to map reads to a reference_

Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 

### [MAFFT](http://mafft.cbrc.jp/alignment/software/)
_Used by `merge_primerid_read_groups.pl`_

Katoh, Standley 2013 MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution 30:772-780.  PMID: 23329690. 

#### [PANDAseq](https://github.com/neufeld/pandaseq) 
_Included for convenience to merge overlapping paired-end reads_

Andre P Masella, Andrea K Bartram, Jakub M Truszkowski, Daniel G Brown and Josh D Neufeld. PANDAseq: paired-end assembler for illumina sequences. BMC Bioinformatics 2012, 13:31. http://www.biomedcentral.com/1471-2105/13/31 PMID: 22333067.
Source code for v2.9 is included as a git submodule.

#### [BEDTools](https://github.com/arq5x/bedtools2)
_Used by `graph_coverage.pl`, `get_majority_block_bam.pl`, and `get_majority_block_bam.pl`_

Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842. PMID: 20110278.

#### [Btrim](http://graphics.med.yale.edu/trim/)
_Included for convenience to trim primer sequences from reads and to split reads into separate amplicon groups_

Kong, Y (2011) Btrim: A fast, lightweight adapter and quality trimming program for next-generation sequencing technologies. Genomics, 98, 152-153. http://dx.doi.org/10.1016/j.ygeno.2011.05.009 PMID: 21651976.

#### [Picard](https://github.com/broadinstitute/picard)
_Used by `get_majority_block_bam.pl` and included for convenience to manipulate SAM/BAM files_

Cite: https://broadinstitute.github.io/picard/

#### [Seqtk](https://github.com/lh3/seqtk)
_Used by `merge_primerid_read_groups.pl`_

Cite: https://github.com/lh3/seqtk

#### [ClustalW2](ftp://ftp.ebi.ac.uk/pub/software/clustalw2/)
_Used `merge_primerid_read_groups.pl` when option `--clustalw` is selected_

Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA, McWilliam H, Valentin F, Wallace IM, Wilm A, Lopez R, Thompson JD, Gibson TJ and Higgins DG. Bioinformatics 2007 23(21): 2947-2948. doi:10.1093/bioinformatics/btm404 PMID: 17846036 

#### [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
_`fastx_trimmer` included for convenience to trim nucleotides from reads_

Cite: https://github.com/agordon/fastx_toolkit

#### [TopHat2](https://ccb.jhu.edu/software/tophat/manual.shtml)
_`bam2fastx` used by `merge_primerid_read_groups.pl`_

Kim D, Pertea G, Trapnell C, Pimentel H, Kelley R, Salzberg SL. TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. Genome Biology 2013, 14:R36. PMID: 23618408.




# Notes for running on NIAID Locus.
v0.1.0 was tested to work on NIAID Locus.
Modules to load prior to running scripts:
* module load R/3.2.3-goolf-1.7.20-2016-Q1
* module load Perl/5.22.1-goolf-1.7.20-2016-Q1
* module load MAFFT/7.221-goolf-1.7.20-with-extensions





# Todo
* Description of scripts 
* Additional steps in the tutorial
