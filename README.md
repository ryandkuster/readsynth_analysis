# **read simulation**

*note: this file corresponds to the statonlabprivate wiki titled 'Read simulation' as well as the github repository /read_simulation stored on flora*

# **directory structure**
- README.md **(you're reading me)**
- analyses/ **(all specific objectives for processing data)**
  - 1520_genomes_RAD_characteristics
    - 2021_06_28_test_100_genomes
    - 2021_06_28_test_1520_genomes
    - 2021_06_30_quantify_30_genomes
- code/ **(specific software used for analyses)**
- raw_data/ **(unprocessed data)**
  - 1520_genomes/
- misc/
  - visuals/ **(visuals presented in this markdown file)**

# **raw_data**

## **1520_genomes**

publication:
https://doi.org/10.1038/s41587-018-0008-8

bioproject:
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA482748

This dataset is relevant for simulating reads of similar sequencing preparation (human gut microbiota).

## 2021.06.14
---
## **automate download of all accessions**

Installed the 1,520 genome accessions produced in https://doi.org/10.1038/s41587-018-0008-8; wget called on the list of accession numbers. 

PRJNA482748 assembly details list provided as sys.argv[1] in download_genomes.py.

This script downloads (wget) all assembled genomes as 'fna.gz' into the genomes/ directory and creates a collection of md5checksums for each.

## 2021.06.16
---
## **summarize 1520 accessions**

Automatically download all accessions and provide summary statistics on the assembly.

Notes: Summary statistics for these assemblies uses python scripts to pull the data. The scripts are stored in the /raw_data/genomes/1520_genomes directory with documentation.

PRJNA482748 accession sheet provided as sys.argv[1] in download_assembly_stats.py.

This script downloads (wget) all GenBank assembly statistics for each of the 1,520 genomes listed into the assembly_stats/ directory and produces a summary csv of all the accessions.

# **analyses**

## 2021.06.28
---
## **in silico digest each of the 1,520 genomes using MseI and EcoRI motifs using 'readsynth' scripts**
*<sub><sup>analyses/1520_genomes_RAD_characteristics/2021_06_28_test_1520_genomes/1_RE_digest_genomes</sub></sup>*

readsynth.py from commit 1e72bf20d58ef8b7af03edc661cca810a3598b85

For the purposes of exploratory analysis, only RE digestion (complete) will be performed using complete digestion. The MseI and EcoRI have no issue of overlapping recognition sites, but this needs to be addressed for future simulations, as duplicate loci will be represented in these events.

The output of each simulation will consist of 4 files:
- raw_digest...csv (the absolute, raw fragments produced by the simulated RE digest) 
- hist_raw_digest...pdf (histogram of the raw fragments)
- sampled...csv (the absolute, raw fragments produced by the simulated RE digest)
- hist_sampled (histogram of the sampled fragments)

*note: The raw_digest...csv files will contain all possible cut sites considering both template and reverse-complemented template orientations. The originating strand is represented in the 'strand' column using a '+' (template) or '-' (reverse). Therefore, if 10,000 fragments are produced, we can assume that there are truly only 5,0000 fragments, as long as digests are complete and there is no strict orientation of cut sites on the 5' or 3' end. If there is an orientation requirement (e.g. adapter ligation requires that MseI be only on the 5' end and EcoRI can only be on the 3' end) then it makes sense to store reads in this strand-specific manner.*

By default the maximum fragment length is set to mean + 6sd, so for these test, using mean of 400 and standard deviation of 100, we can expect the maximum fragment length to be included to max out at 1000bp.


## 2021.06.29
---
## **get summary statistics of the 1,520 digests**
*<sub><sup>analyses/1520_genomes_RAD_characteristics/2021_06_28_test_1520_genomes/2_pull_digest_stats</sub></sup>*

for every template strand ('+') from the raw_digest...csv file, pull the total number of fragments as well as ranges of fragment sizes (e.g. 1 to 200bp, 201 to 400bp, 401 to 600bp).

The R file 1520_genomes_stats.R creates correlation plots for the read fragment lengths and various characteristics of the genome assemblies and genome size and gc content.

![correlation plot of all variables](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/correlation.png)

After normalizing the fragment counts per genome (count/genome(bp)), the most notable relationship is the negative correlation of gc_ratio and framents in the 1 to 200 bp range:

![scatter plot of gc and 1 - 200 bp fragments](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/gc_ratio_and_fragments_plot.png)


## 2021.06.30
---
*<sub><sup>analyses/1520_genomes_RAD_characteristics/2021_06_30_quantify_30_genomes/1_create_genome_keys</sub></sup>*

## **produce community profile**
select 30 genomes and randomly assign copy number to each, then simulate fastq files for each abundance profile provided in the key

first, get approximate sampling rate:

```
python3 -c "print(str(30/1520))"
0.019736842105263157
```

so we'll use ~0.02 for sampling genomes randomly

```
cat file_list.txt | awk -v var="0.02" 'BEGIN {srand()} !/^$/ { if (rand() <= var) print $0}' > sampled_files.txt
```

run 'genome_copy_number.py' to randomly generate copy numbers between 1 and 100 per sample and save as a tab-separated key for benchmarking the profiling success of RAD libraries.

## **produce fastq files from community profile**
*<sub><sup>analyses/1520_genomes_RAD_characteristics/2021_06_30_quantify_30_genomes/2_simulate_fastq_files</sub></sup>*

copy or link the 30 chosen files to the local 'genomes' directory

use 'run_quantify_sim.sh' script to automate the fastq simulation of each of the 30 genomes (all files written to 'output' directory)

readsynth.py from commit 9b95510ebd1d25a45310d8bd8b60caf4af9e63c5

concatenate the resulting 30 fastq files into a single file:

```
cat ./output/*1.fastq > 2021_06_30_combined_R1.fastq
cat ./output/*2.fastq > 2021_06_30_combined_R2.fastq
```

## 2021.07.06
---
## **recreate community profiling using updated script**

readsynth.py from commit 4df8ce1f90394987637b4a83443ec53e808c1af2

The readsynth.py code was updated to more accurately represent the composition of the original genome size. For example, a complete digest with no R1/R2 orientation will produce an equal number of fragments on the template and non-template strands, so the sampled distribution should represent half this distribution (shown below in orange).

![sampling from complete digest at 1 X](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/hist_sampled_half_chr_1_complete.png)

Similarly, using an incomplete digest will produce an overabundance of fragments from some parts of the genome, and a weighting system must account for this. The figure below details an incomplete digest distribution (blue) with the weight-adjusted sampling distribution approximating 1x representation of the original genome (note the count similarity to the complete digest above).

![sampling from incomplete digest at 1 X](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/hist_sampled_half_chr_1_incomplete.png)

Importantly, simulating complete EcoR1/MseI digests where R1/R2 adapters ligate in a 5' to 3' specific manner requires weight consideration to accurately approximate genome copy number. Complete digests with orientation won't encounter overlap at any point on either template or non-template strands.

An example from today's simulation from the GCA_003468235 accession using 4 as the original genome copy number:

![4 X coverage example](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/GCA_003468235_4X.png)

Finally, a genome copy number of 20 for accession GCA_003433755:

![20 X coverage example](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/GCA_003433755_20X.png)

## 2021.07.06
---
## **community profiling using Kraken2**

*<sub><sup>/pickett_flora/projects/read_simulation/analyses/1520_genomes_RAD_characteristics/2021_06_30_quantify_30_genomes/3_profile_simulated_samples/</sub></sup>*

https://github.com/DerrickWood/kraken2/wiki

```
spack list | grep 'kraken'
kraken
kraken2
```

kraken not installed:

```
spack load kraken2@2.0.8-beta
==> Error: Spec 'kraken2@2.0.8-beta' matches no installed packages.
```

Attempt manual install using spack:

```
spack install kraken2@2.0.8-beta
==> Warning: gcc@8.3.1 cannot build optimized binaries for "zen2". Using best target possible: "zen"
==> Error: Failed to acquire a write lock for berkeley-db-18.1.40-nebwsa4tt25wkohwwlzeikbln2unmipv due to LockROFileError: Can't take write lock on read-only file: /pickett_shared/spack/opt/spack/.spack-db/prefix_lock
==> Error: Can't take write lock on read-only file: /pickett_shared/spack/opt/spack/.spack-db/prefix_lock
```

I installed version 2.0.7-beta using conda:

```
conda install -c bioconda kraken2
kraken2 --version
Kraken version 2.0.7-beta
Copyright 2013-2018, Derrick Wood (dwood@cs.jhu.edu)
```

Then built the kraken2 database using the subsetted bacteria library:

```
mkdir bacteria_db
kraken2-build --download-library bacteria --db ./bacteria_db/
Step 1/2: Performing rsync file transfer of requested files                                                                                                    
Rsync file transfer complete.                                                                                                                                  
Step 2/2: Assigning taxonomic IDs to sequences                                                                                                                 
Processed 25876 projects (59362 sequences, 107.02 Gbp)... done.                                                                                                
All files processed, cleaning up extra sequence files... done, library complete.                                                                               
Masking low-complexity regions of downloaded library... done.
```

version 2.0.7-beta is no longer compatible with downloading from the NCBI database and building the kraken2 database with taxonomy

```
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh /pickett_flora/projects/read_simulation/code/kraken2
```

Attempt building the standard database after uninstalling the conda kraken2 packaged (important as it has global tools in PATH that are outdated). To avoid installing the NCBI dustmasker and segmasker tools, use the --no-masking argument.

```
mkdir standard_db
/pickett_flora/projects/read_simulation/code/kraken2/kraken2-build --standard --no-masking --threads 10 --db standard_db/
```

This failed for rsync reasons documented in current GitHub issues. A simple workaround is to download premade databases stored at https://benlangmead.github.io/aws-indexes/k2.


## 2021.07.07
---
## **community profiling using Kraken2 and Bracken**

using Kraken2 from commit 84b2874e0ba5ffc9abaebe630433a430cd0f69f4
using Bracken from commit b014034d5c0efa7c4cb3eb9c1f0913c09cdce728

```
cd ./standard_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz
tar -xvzf k2_standard_16gb_20210517.tar.gz
rm k2_standard_16gb_20210517.tar.gz
cd ..
```

Using the above library, align the paired-end fastq mock community reads created with readsynth.py.

```
cp ../2_simulate_fastq_files/2021_06_30_combined_R1.fastq ./
cp ../2_simulate_fastq_files/2021_06_30_combined_R2.fastq ./
/pickett_flora/projects/read_simulation/code/kraken2/kraken2 --paired --threads 12 --report 2021_07_07.kreport --classified-out cseqs#.fq --unclassified-out useqs#.fq --db /pickett_flora/projects/read_simulation/raw_data/standard_db/ 2021_06_30_combined_R1.fastq 2021_06_30_combined_R2.fastq > 2021_07_07_run.txt
Loading database information... done.                                                                                                                          
875990 sequences (438.00 Mbp) processed in 1.284s (40919.3 Kseq/m, 20459.63 Mbp/m).
  679276 sequences classified (77.54%)
  196714 sequences unclassified (22.46%)
```

Using the README at https://github.com/jenniferlu717/Bracken/ as a guide.

Run Bracken for abundance estimation using 10 as the default -t threshold level ("any species with <= 10 (or otherwise specified) reads will not receive any additional reads from higher taxonomy levels when distributing reads for abundance estimation"). The -r argument is read length (our mock communities are at 250bp), and the -l level value of S for species.

```
/pickett_flora/projects/read_simulation/code/Bracken/bracken -d /pickett_flora/projects/read_simulation/raw_data/standard_db/ -i 2021_07_07.kreport -o 2021_07_07.bracken -r 250 -l S -t 10
/pickett_flora/projects/read_simulation/code/Bracken/bracken -d /pickett_flora/projects/read_simulation/raw_data/standard_db/ -i 2021_07_07.kreport -o 2021_07_07_G.bracken -r 250 -l G -t 10
```

## **comparing Bracken abundance profiles with fastq files and input abundance keys**
*<sub><sup>/pickett_flora/projects/read_simulation/analyses/1520_genomes_RAD_characteristics/2021_06_30_quantify_30_genomes</sub></sup>*

First, the abundance of reads in the fastq file must be considered as separate from the key file, as GC content of the source genome is expected to skew the number of fragments. Second, the original key file needs to be reflected as a percentage of total reads on a per taxid basis.

Grab the header lines only (R2 file will be the same, so only R1 is needed):

```
cp ../2_simulate_fastq_files/2021_06_30_combined_R1.fastq ./
sed -n '1~4p' 2021_06_30_combined_R1.fastq > R1_headers.txt
```

Then, run 'count_fq_taxids.py' to collect the relevant information from the key file and the previously created digested_genome_stats.csv file from 2021_06_28:

```
cp ../../2021_06_28_test_1520_genomes/2_pull_digest_stats/digested_genome_stats.csv ./
python3 /pickett_flora/projects/read_simulation/code/scripts/count_fq_taxids.py R1_headers.txt sampled_files_key.txt digested_genome_stats.csv
```

The resulting sampled_genome_stats.csv will contain the original abundances from the key ('copies'), the ratio of these copies vs. the total copies produced in the key file ('copy_ratio'), the count of fastq reads per sample ('reads'), and these reads as a ratio to the total reads in the simulated fastq file ('read_ratio').

![correlation plot comparing original key vs. simulated fastq abundance profiles](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/abundance_corr.png)

To pull the Bracken abundances (only at the genus level):

```
python3 /pickett_flora/projects/read_simulation/code/scripts/pull_bracken_levels.py ../3_profile_simulated_samples/2021_07_07_bracken_genuses.kreport G sampled_genome_stats.csv
```

This will produce the profiled_genome_stats.csv file.

After collapsing the 30 reads into the genus level, the following visuals were produced. In general, the Bracken results closely match the expected values in most genera. We can see four groups (Ruminococcus, Lachnospira, Blautia, Clostridium, and an unidentified genome) that were detected at much lower levels than expected.

![correlation plot comparing original key copies with the bracken output](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/bracken_to_reads_plot.png)

The following plot shows the underlying relationship between the input key copies and the resulting fastq read counts:

![correlation plot comparing original key copies with the reads](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/reads_to_copy_plot.png)

After considering GC content, we see there is an underlying issue in the dataset. The original key file, by chance, gave lower abundances to GC-rich organisms.

![correlation plot comparing original key copies with the bracken output](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/gc_ratio_in_key_plot.png)

However, this does not fully explain the outcome of Bracken. Bracken did appear to perform best with genomes containing lower GC content.

![correlation plot comparing original key copies with the bracken output](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/3d_plot.png)

## 2021.07.09
---
## **profiling Liu 2017 datasets using Kraken2 and Bracken**

Download the datasets:

```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5298272/SRR5298272
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5298274/SRR5298274
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5360684/SRR5360684

spack load sratoolkit@2.10.7

fastq-dump SRR5298272
fastq-dump SRR5298274
fastq-dump SRR5360684
```

The files are paired-end, but can't be split using the expected fastq-dump flags, so we'll just manually do it:

```
sed -n '1~8p;2~8p;3~8p;4~8p' SRR5298272.fastq > SRR5298272_R1.fastq                                                           
sed -n '5~8p;6~8p;7~8p;8~8p' SRR5298272.fastq > SRR5298272_R2.fastq
sed -n '1~8p;2~8p;3~8p;4~8p' SRR5298274.fastq > SRR5298274_R1.fastq                                                           
sed -n '5~8p;6~8p;7~8p;8~8p' SRR5298274.fastq > SRR5298274_R2.fastq
sed -n '1~8p;2~8p;3~8p;4~8p' SRR5360684.fastq > SRR5360684_R1.fastq                                                           
sed -n '5~8p;6~8p;7~8p;8~8p' SRR5360684.fastq > SRR5360684_R2.fastq
```

There. But wait, the header naming won't work with bwa mem:

```
for i in *R2.fastq ; do python3 ../../code/scripts/fix_R2_headers2.py $i && mv modified_${i} $i ; done
```

*<sub><sup>/pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/1_profile_reads</sub></sup>*

Run Kraken as before:

```
/pickett_flora/projects/read_simulation/code/kraken2/kraken2 --paired --threads 12 --report 2021_07_09.kreport --classified-out cseqs#.fq --unclassified-out useqs#.fq --db /pickett_flora/projects/read_simulation/raw_data/standard_db/ SRR5298272_R1.fastq SRR5298272_R2.fastq > 2021_07_09_run.txt

/pickett_flora/projects/read_simulation/code/Bracken/bracken -d /pickett_flora/projects/read_simulation/raw_data/standard_db/ -i 2021_07_09.kreport -o 2021_07_09.bracken -r 75 -l S -t 10
```

Find the top hits:

```
python3 /pickett_flora/projects/read_simulation/code/scripts/pull_kraken_taxids.py cseqs_1.fq
python3 /pickett_flora/projects/read_simulation/code/scripts/pull_kraken_taxids.py cseqs_2.fq
```

The three most prevalent taxa were 816 (Bacteroides), 853 (Faecalibacterium prausnitzii), and 74426 (Collinsella aerofaciens).

*<sub><sup>/pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/2_simulate_f_prausnitzii</sub></sup>*

Download F. prausnitzii

```
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/312/465/GCF_003312465.1_ASM331246v1/GCF_003312465.1_ASM331246v1_genomic.fna.gz"
```

Run readsynth.py (commit 4df8ce1f90394987637b4a83443ec53e808c1af2)

```
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/2_simulate_f_prausnitzii
python3 /home/rkuster/readsynth/readsynth.py \
   -genome GCF_003312465.1_ASM331246v1_genomic.fna \
   -o ./output \
   -m1 CATG/ \
   -m2 A/CGT \
   -n 1 \
   -l 76 \
   -complete 1 \
   -mean 550 \
   -sd 50 \
   -a1 /pickett_flora/projects/read_simulation/raw_data/adapters/liu_ddRADseq/liu_ddRADseq_adapters_R1.txt \                                               
   -a2 /pickett_flora/projects/read_simulation/raw_data/adapters/liu_ddRADseq/liu_ddRADseq_adapters_R2.txt \                                               
   -a1s 72 \
   -a2s 71
```

## 2021.07.12
---
## **map reads from kraken assignments and reads from above simulation of F. prausnitzii to compare sites and depth**

*<sub><sup>/pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping</sub></sup>*

Grab the reads for each of the three Liu datasets that Kraken identified as F. prausnitzii.

```
python3 /pickett_flora/projects/read_simulation/code/scripts/pull_kraken_reads.py ../../1_profile_reads/SRR5298272/cseqs_1.fq 853
python3 /pickett_flora/projects/read_simulation/code/scripts/pull_kraken_reads.py ../../1_profile_reads/SRR5298272/cseqs_2.fq 853
```

To BWA map the reads without error due to fastq naming discrepency, rename the R2 headers to be identical with R1:

```
python3 /pickett_flora/projects/read_simulation/code/scripts/fix_R2_headers.py taxid_853_cseqs_2.fq
```

BWA map reads from RMS project and reads simulated from readsynth.py to the F. prausnitzii genome. The reads are then subsetted to only those that have proper pairing.

```
spack load bwa@0.7.17
spack load samtools@1.10
cp ../2_simulate_f_prausnitzii/GCF_003312465.1_ASM331246v1_genomic.fna ./
bwa index GCF_003312465.1_ASM331246v1_genomic.fna

cd SRR5298272
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna taxid_853_cseqs_1.fq modified_taxid_853_cseqs_2.fq > SRR5298272.sam
samtools view -b -f 0x2 SRR5298272.sam > SRR5298272.bam
samtools sort -o SRR5298272.sorted.bam SRR5298272.bam
samtools index SRR5298272.sorted.bam
samtools depth -o SRR5298272.depth.txt -H SRR5298272.sorted.bam
samtools stats SRR5298272.sorted.bam > SRR5298272_stats.txt
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py SRR5298272.depth.txt
12.762109703977144
cd ..

cd SRR5298274
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna taxid_853_cseqs_1.fq modified_taxid_853_cseqs_2.fq > SRR5298274.sam
samtools view -b -f 0x2 SRR5298274.sam > SRR5298274.bam
samtools sort -o SRR5298274.sorted.bam SRR5298274.bam
samtools index SRR5298274.sorted.bam
samtools depth -o SRR5298274.depth.txt -H SRR5298274.sorted.bam
samtools stats SRR5298274.sorted.bam > SRR5298274_stats.txt
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py SRR5298274.depth.txt
4.958139806567638
cd ..

cd SRR5360684
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna taxid_853_cseqs_1.fq modified_taxid_853_cseqs_2.fq > SRR5360684.sam
samtools view -b -f 0x2 SRR5360684.sam > SRR5360684.bam
samtools sort -o SRR5360684.sorted.bam SRR5360684.bam
samtools index SRR5360684.sorted.bam
samtools depth -o SRR5360684.depth.txt -H SRR5360684.sorted.bam
samtools stats SRR5360684.sorted.bam > SRR5360684_stats.txt
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py SRR5360684.depth.txt
9.023378056269712
cd ..

mkdir simulated_F_prausnitzii
cd simulated_F_prausnitzii
cp ../../2_simulate_f_prausnitzii/output/GCF_003312465.1_ASM331246v1_genomic.fna_R\* ./
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna GCF_003312465.1_ASM331246v1_genomic.fna_R1.fastq GCF_003312465.1_ASM331246v1_genomic.fna_R2.fastq > simulated.sam
samtools view -b -f 0x2 simulated.sam > simulated.bam                                                                
samtools sort -o simulated.sorted.bam simulated.bam                                                           
samtools index simulated.sorted.bam 
samtools depth -o simulated.depth.txt -H simulated.sorted.bam
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py simulated.depth.txt
1.3841179998500637
```

||SRR5298272|SRR5298274|SRR5298274|simulated|
|:-|:-|:-|:-|:-|
|reads|184856|31108|86706|1958|
|insert (bp)|317.0|321.5|299.1|383.9|
|insert sd (bp)|119.0|115.2|125.9|73.0|

The weighted average insert for the SRA data is 312.3346833.

The weighted average insert sd for the SRA data is 120.5860872.


## 2021.07.13
---
## **repeat simulation using error profile and n = 13 from SRR5298272 results and incomplete digest**


Run readsynth.py (commit 4cf90a729ccd292c4c78bcfff27f3b546bcfcd50)

```
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/2_simulate_f_prausnitzii
python3 /home/rkuster/readsynth/readsynth.py \
  -genome GCF_003312465.1_ASM331246v1_genomic.fna \
  -o ./output2/ \
  -m1 CATG/ \
  -m2 A/CGT \
  -n 13 \
  -l 76 \
  -complete 0 \
  -mean 550 \
  -sd 50 \
  -a1 /pickett_flora/projects/read_simulation/raw_data/adapters/liu_ddRADseq/liu_ddRADseq_adapters_R1.txt \
  -a2 /pickett_flora/projects/read_simulation/raw_data/adapters/liu_ddRADseq/liu_ddRADseq_adapters_R2.txt \
  -a1s 72 \
  -a2s 71 \
  -r1 ../1_profile_reads/SRR5298272/SRR5298272_R1.fastq \
  -r2 ../1_profile_reads/SRR5298272/SRR5298272_R2.fastq \
  -p 1 \
  -q1 qscores.SRR5298272_R1.fastq.csv \
  -q2 qscores.SRR5298272_R2.fastq.csv
```

```
spack load bwa@0.7.17
spack load samtools@1.10

cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/
mkdir simulated_F_prausnitzii_2
cd simulated_F_prausnitzii_2
cp ../../2_simulate_f_prausnitzii/output2/GCF_003312465.1_ASM331246v1_genomic.fna_R* ./
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna GCF_003312465.1_ASM331246v1_genomic.fna_R1.fastq GCF_003312465.1_ASM331246v1_genomic.fna_R2.fastq > simulated2.sam
samtools view -b -f 0x2 simulated2.sam > simulated2.bam                                                                
samtools sort -o simulated2.sorted.bam simulated2.bam                                                           
samtools index simulated2.sorted.bam 
samtools depth -o simulated2.depth.txt -H simulated2.sorted.bam
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py simulated2.depth.txt
3.6272405227291276
```

## 2021.07.14
---
## **bwa map Liu datasets to reference without subsetting kraken F. prausnitzii hits**

```
spack load bwa@0.7.17
spack load samtools@1.10

cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping

mkdir SRR5298272_non_profiled
cd SRR5298272_non_profiled
ln -fs /pickett_flora/projects/read_simulation/raw_data/liu_RMS/SRR5298272_R* ./
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna SRR5298272_R1.fastq SRR5298272_R2.fastq > SRR5298272_non_profiled.sam
samtools view -b -f 0x2 SRR5298272_non_profiled.sam > SRR5298272_non_profiled.bam
samtools sort -o SRR5298272_non_profiled.sorted.bam SRR5298272_non_profiled.bam
samtools index SRR5298272_non_profiled.sorted.bam
samtools stats SRR5298272_non_profiled.sorted.bam > SRR5298272_non_profiled_stats.txt

mkdir SRR5298274_non_profiled
cd SRR5298274_non_profiled
ln -fs /pickett_flora/projects/read_simulation/raw_data/liu_RMS/SRR5298274_R* ./
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna SRR5298274_R1.fastq SRR5298274_R2.fastq > SRR5298274_non_profiled.sam
samtools view -b -f 0x2 SRR5298274_non_profiled.sam > SRR5298274_non_profiled.bam
samtools sort -o SRR5298274_non_profiled.sorted.bam SRR5298274_non_profiled.bam
samtools index SRR5298274_non_profiled.sorted.bam
samtools stats SRR5298274_non_profiled.sorted.bam > SRR5298274_non_profiled_stats.txt

mkdir SRR5360684_non_profiled
cd SRR5360684_non_profiled
ln -fs /pickett_flora/projects/read_simulation/raw_data/liu_RMS/SRR5360684_R* ./
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna SRR5360684_R1.fastq SRR5360684_R2.fastq > SRR5360684_non_profiled.sam
samtools view -b -f 0x2 SRR5360684_non_profiled.sam > SRR5360684_non_profiled.bam
samtools sort -o SRR5360684_non_profiled.sorted.bam SRR5360684_non_profiled.bam
samtools index SRR5360684_non_profiled.sorted.bam
samtools stats SRR5360684_non_profiled.sorted.bam > SRR5360684_non_profiled_stats.txt
```

## 2021.07.15
---
## *compare correlation between simulated and actual RMS reads*

```
spack load samtools@1.10
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping
samtools depth -o compare_all.depth.txt -H ./simulated_F_prausnitzii/simulated.sorted.bam ./simulated_F_prausnitzii_2/simulated2.sorted.bam ./SRR5298272/SRR5298272.sorted.bam ./SRR5298274/SRR5298274.sorted.bam ./SRR5360684/SRR5360684.sorted.bam
```

used compare_sim_to_liu.R script in Rstudio to do a simple correlation of the simulated and Liu 2017 depths

## 2021.07.16
---
## *test impact of altering genome copy number (-n)*

Increasing the genome copy number seems to have a strong effect on the number of shorter reads drawn in the incomplete digest simulation. If the upper limit of fragment size (usually mean + 6sd) is very large (can be set manually with -f), then the shorter reads will be penalized more in the current weighting scheme. Even with the mean + 6sd condition, many smaller fragments will be rare unless the copy number is increased.

The goal is to keep -f fixed and increase copy number (-n) by factors of 10 at 1, 10, 100.

Also, the mean and sd will be altered to simulate the Liu 2017 F. prausnitzii dataset (mean: 143 adapter length + 312 = 455; sd: 121)

using readsynth.py (commit eba1aedcc24555c0201c8d5ce5ac0d9ac8a37fda)

```
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_16_compare_copy_number
cp -r ../2021_07_09_simulate_top_hit/2_simulate_f_prausnitzii/ ./1_simulate_f_prausnitzii
rm -rf output*
mkdir output_n_1
mkdir output_n_10
mkdir output_n_100
```

BWA align the reads to the F. prausnitzii reference.

```
spack load bwa@0.7.17
spack load samtools@1.10

cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_16_compare_copy_number
mkdir 2_map_to_ref

cp ../1_simulate_f_prausnitzii/GCF_003312465.1_ASM331246v1_genomic.fna ./
bwa index GCF_003312465.1_ASM331246v1_genomic.fna

mkdir simulation_n_1
mkdir simulation_n_10
mkdir simulation_n_100

cd simulation_n_1
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna ../../1_simulate_f_prausnitzii/output_n_1/GCF_003312465.1_ASM331246v1_genomic.fna_R1.fastq ../../1_simulate_f_prausnitzii/output_n_1/GCF_003312465.1_ASM331246v1_genomic.fna_R2.fastq > simulation_n_1.sam
samtools view -b -f 0x2 simulation_n_1.sam > simulation_n_1.bam
samtools sort -o simulation_n_1.sorted.bam simulation_n_1.bam
samtools index simulation_n_1.sorted.bam

cd ../simulation_n_10
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna ../../1_simulate_f_prausnitzii/output_n_10/GCF_003312465.1_ASM331246v1_genomic.fna_R1.fastq ../../1_simulate_f_prausnitzii/output_n_10/GCF_003312465.1_ASM331246v1_genomic.fna_R2.fastq > simulation_n_10.sam
samtools view -b -f 0x2 simulation_n_10.sam > simulation_n_10.bam
samtools sort -o simulation_n_10.sorted.bam simulation_n_10.bam
samtools index simulation_n_10.sorted.bam

cd ../simulation_n_100
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna ../../1_simulate_f_prausnitzii/output_n_100/GCF_003312465.1_ASM331246v1_genomic.fna_R1.fastq ../../1_simulate_f_prausnitzii/output_n_100/GCF_003312465.1_ASM331246v1_genomic.fna_R2.fastq > simulation_n_100.sam
samtools view -b -f 0x2 simulation_n_100.sam > simulation_n_100.bam
samtools sort -o simulation_n_100.sorted.bam simulation_n_100.bam
samtools index simulation_n_100.sorted.bam
```

Alignments viewed with IGV. Many trends carry over between the three levels of representation. The reads don't appear to follow the empirical data from Liu 2017.

## 2021.07.20
---
## *test new copy number simulation variation*

Using the updated readsynth partial_digest script, re-digest the F. prausnitzii genome using the same restriction enzymes as before from the Liu 2017 study (NlaIII and HpyCH4IV) as well as the approximate read length, standard deviation, and copy number (13, based on depth)

Sample distributions using F. prausnitzii:

![10 X copy number with .2 digest efficiency](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/20_n10_sampled_GCF_000005845.2_ASM584v2_genomic.fna.png)

![10 X copy number with .4 digest efficiency](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/40_n10_sampled_GCF_000005845.2_ASM584v2_genomic.fna.png)

![10 X copy number with .6 digest efficiency](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/60_n10_sampled_GCF_000005845.2_ASM584v2_genomic.fna.png)

![10 X copy number with .8 digest efficiency](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/80_n10_sampled_GCF_000005845.2_ASM584v2_genomic.fna.png)

![10 X copy number with 1 digest efficiency](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/100_n10_sampled_GCF_000005845.2_ASM584v2_genomic.fna.png)

Note how the size-selected reads (green) represent a normal distribution that intersects with the fragments in the defined size range. This is similar to the size selection reported in Peterson et al. 2012.

![Peterson et al. 2012 size distribution](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/peterson_2012.png)

This also follows the Sage Science BluePippin size distribution:

![Sage Science size distribution](https://github.com/ryandkuster/read_simulation/blob/main/misc/visuals/sage_bluepippin_2016.png)

readsynth.py commit 6d7984ab9086541cbdd4bd62bba6e1aaf22a4c6a

```
python3 /home/rkuster/readsynth/readsynth.py \
  -genome GCF_003312465.1_ASM331246v1_genomic.fna \
  -o ./output/ \
  -m1 CATG/ \
  -m2 A/CGT \
  -n 13 \
  -l 76 \
  -complete 0.9 \
  -mean 455 \
  -sd 121 \
  -a1 /pickett_flora/projects/read_simulation/raw_data/adapters/liu_ddRADseq/liu_ddRADseq_adapters_R1.txt \
  -a2 /pickett_flora/projects/read_simulation/raw_data/adapters/liu_ddRADseq/liu_ddRADseq_adapters_R2.txt \
  -a1s 72 \
  -a2s 71 \
  -t 1 \ 
```
map

```
spack load bwa@0.7.17
spack load samtools@1.10
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_20_updated_copy_number/2_map_to_ref/simulation_n_13
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna ../../1_simulate_f_prausnitzii/output/GCF_003312465.1_ASM331246v1_genomic.fna_R1.fastq ../../1_simulate_f_prausnitzii/output/GCF_003312465.1_ASM331246v1_genomic.fna_R2.fastq > simulation_n_13.sam
samtools view -b -f 0x2 simulation_n_13.sam > simulation_n_13.bam
samtools sort -o simulation_n_13.sorted.bam simulation_n_13.bam
samtools index simulation_n_13.sorted.bam

cd ..
samtools depth -o compare_all.depth.txt -H ./simulation_n_13/simulation_n_13.sorted.bam \
  ../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298272_non_profiled/SRR5298272_non_profiled.sorted.bam \
  ../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298274_non_profiled/SRR5298274_non_profiled.sorted.bam \
  ../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5360684_non_profiled/SRR5360684_non_profiled.sorted.bam \
  ../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/simulated_F_prausnitzii/simulated.sorted.bam \
  ../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/simulated_F_prausnitzii_2/simulated2.sorted.bam

echo "CHROM  POS 90_n13_sim SRR5298272_non_profiled SRR5298274_non_profiled SRR5360684_non_profiled old_sim1 old_sim2" > compare_depths.txt
sed -n '2,$p' compare_all.depth.txt >> compare_depths.txt
```

These results only show on nominal improvement vs. the first simulation (complete digest) of R = .35 vs R = .33. It is possible that copy number is affecting these outcomes. Assessing the SRR5298274_non_profiled results:

```
samtools depth -o SRR5298272_non_profiled.depth.txt -H SRR5298272_non_profiled.sorted.bam
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298272_non_profiled
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py SRR5298272_non_profiled.depth.txt 
16.194673036059225
```

The most recent simulation used n=13 based on an earlier depth reading for the profiled (possibly incorrectly paired) SRR5298272 results. To confirm that average depth actually relates to the input n=13 of readsynth.py, assess read depth results of the simulated reads.

```
spack load samtools@1.10
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_20_updated_copy_number/2_map_to_ref/simulation_n_13
samtools depth -o simulation_n_13.depth.txt -H simulation_n_13.sorted.bam
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py simulation_n_13.depth.txt 
7.644556639965016
```

Takeaway: average depth doesn't seem to relate to the copy number input into readsynth.py.

install ddRAGE version 1.7.1

```
conda create -n ddrage -c bioconda ddrage
source activate ddrage
```

## 2021.07.23
---
## comparing F. prausnitzii alignments using Jaccard distance

Pull only the sequences from the alignments for comparison:

```
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298272_non_profiled
samtools view SRR5298272_non_profiled.sorted.bam | awk '{ print $10 }' > SRR5298272_non_profiled_sequences.txt

cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2021_07_20_updated_copy_number/2_map_to_ref/simulation_n_13
samtools view simulation_n_13.sorted.bam | awk '{ print $10 }' > simulation_n_13_sequences.txt
```

Run custom python script 'jaccard_distance.py' using k-mer size 5 to compare recent simulations with a previous file mapped using bwa (not Kraken) to the F. prausnitzii reference:

```
python3 ../../../../../code/scripts/jaccard_distance.py simulation_n_13_sequences.txt ../../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298274_non_profiled/SRR5298274_non_profiled_sequences.txt 5
```
*getting kmers for first file  
getting kmers for second file  
simulation_n_13_sequences.txt vs SRR5298274_non_profiled_sequences.txt 0.4176675419966591*

```
python3 ../../../../../code/scripts/jaccard_distance.py ../../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298272_non_profiled/SRR5298272_non_profiled_sequences.txt ../../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298274_non_profiled/SRR5298274_non_profiled_sequences.txt 5
```

*getting kmers for first file  
getting kmers for second file  
SRR5298272_non_profiled_sequences.txt vs SRR5298274_non_profiled_sequences.txt 0.1776127556747974*

```
python3 ../../../../../code/scripts/jaccard_distance.py ../../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298272_non_profiled/SRR5298272_non_profiled_sequences.txt ../../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5360684_non_profiled/SRR5360684_non_profiled_sequences.txt 5
```

*getting kmers for first file  
getting kmers for second file  
SRR5298272_non_profiled_sequences.txt vs SRR5360684_non_profiled_sequences.txt 0.3746411604853273*

Repeat using k=20:

```
python3 ../../../../../code/scripts/jaccard_distance.py simulation_n_13_sequences.txt ../../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298274_non_profiled/SRR5298274_non_profiled_sequences.txt 20
```

*getting kmers for first file  
getting kmers for second file  
simulation_n_13_sequences.txt vs SRR5298274_non_profiled_sequences.txt 0.05765676657390458*

```
python3 ../../../../../code/scripts/jaccard_distance.py ../../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298272_non_profiled/SRR5298272_non_profiled_sequences.txt ../../../2021_07_09_simulate_top_hit/3_compare_rms_to_sim_mapping/SRR5298274_non_profiled/SRR5298274_non_profiled_sequences.txt 20
```

*getting kmers for first file  
getting kmers for second file  
SRR5298272_non_profiled_sequences.txt vs SRR5298274_non_profiled_sequences.txt 0.0791834345843937*

Jaccard distance (modified to count repeated k-mers) shows a higher similarity to

## 2021.07.26
---
## use Snipen BEI mock sequencing data to compare effects of RE choice and digest efficiency

using Kraken2 from commit 84b2874e0ba5ffc9abaebe630433a430cd0f69f4
using Bracken from commit b014034d5c0efa7c4cb3eb9c1f0913c09cdce728

```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10199716/SRR10199716
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10199724/SRR10199724
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10199725/SRR10199725
spack load sratoolkit@2.10.7

for i in SRR*; do fastq-dump --split-files $i ; done
```

Use ngscomposer scallop.py (commit 96f847fea7c89f74cc42f52004d43c2b0c5ac88a) trim adapters from reads to leave only genomic RE motif and beyond (EcoRI lost 'G' in ligation process)

```
cd /pickett_flora/projects/read_simulation/analyses/snipen_2021_RMS/2021_07_26_test_mock_data/1_front_trim_snipen_seqs
ln -fs /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/SRR101997*fastq ./
for i in *1.fastq ; do python3 /home/rkuster/ngscomposer/tools/scallop.py -r1 $i -f 11 ; done
for i in *2.fastq ; do python3 /home/rkuster/ngscomposer/tools/scallop.py -r1 $i -f 13 ; done
```
Create a custom kraken2/Bracken database for the 20 genomes in the Bei dataset:

```
cd /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/mock_community_ref_genomes
for i in *fna ; do ../../../code/kraken2/kraken2-build --add-to-library $i --db ../bei_20_genomes_db/ ; done
../../../code/kraken2/kraken2-build --download-taxonomy --db ../bei_20_genomes_db/
../../../code/kraken2/kraken2-build --build --db ../bei_20_genomes_db/
../../../code/Bracken/bracken-build -d ../bei_20_genomes_db/ -t 1 -x ../../../code/kraken2/
```

run kraken taxonomic profiling on reads:

```
cd /pickett_flora/projects/read_simulation/analyses/snipen_2021_RMS/2021_07_26_test_mock_data/2_kraken_snipen_seqs

mkdir SRR10199716
mkdir SRR10199724
mkdir SRR10199725

cd SRR10199716
/pickett_flora/projects/read_simulation/code/kraken2/kraken2 --paired --threads 2 --report 2021_07_28.kreport --classified-out cseqs#.fq --unclassified-out useqs#.fq --db /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/bei_20_genomes_db/ ../../1_front_trim_snipen_seqs/trimmed_se.SRR10199716_1.fastq ../../1_front_trim_snipen_seqs/trimmed_se.SRR10199716_2.fastq > 2021_07_28_run.txt
```

*658194 sequences (179.03 Mbp) processed in 5.211s (7578.1 Kseq/m, 2061.26 Mbp/m).  
  620655 sequences classified (94.30%)  
  37539 sequences unclassified (5.70%)*

```
cd ../SRR10199724
/pickett_flora/projects/read_simulation/code/kraken2/kraken2 --paired --threads 2 --report 2021_07_28.kreport --classified-out cseqs#.fq --unclassified-out useqs#.fq --db /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/bei_20_genomes_db/ ../../1_front_trim_snipen_seqs/trimmed_se.SRR10199724_1.fastq ../../1_front_trim_snipen_seqs/trimmed_se.SRR10199724_2.fastq > 2021_07_28_run.txt
```

*1607175 sequences (437.15 Mbp) processed in 13.629s (7075.6 Kseq/m, 1924.55 Mbp/m).  
  1538745 sequences classified (95.74%)  
  68430 sequences unclassified (4.26%)*

```
cd ../SRR10199725
/pickett_flora/projects/read_simulation/code/kraken2/kraken2 --paired --threads 2 --report 2021_07_28.kreport --classified-out cseqs#.fq --unclassified-out useqs#.fq --db /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/bei_20_genomes_db/ ../../1_front_trim_snipen_seqs/trimmed_se.SRR10199725_1.fastq ../../1_front_trim_snipen_seqs/trimmed_se.SRR10199725_2.fastq > 2021_07_28_run.txt
```

*789892 sequences (214.85 Mbp) processed in 7.217s (6567.2 Kseq/m, 1786.27 Mbp/m).  
  758411 sequences classified (96.01%)  
  31481 sequences unclassified (3.99%)* 

```
/pickett_flora/projects/read_simulation/code/Bracken/bracken -d /pickett_flora/projects/read_simulation/raw_data/standard_db/ -i 2021_07_28.kreport -o 2021_07_28.bracken -r 75 -l S -t 10
```

For each of the 20 genomes, the 16S copy number was referenced at https://rrndb.umms.med.umich.edu as in Snipen et al. 2021.

```
6    GCA_013372085.1_ASM1337208v1_genomic.fna.gz
3    GCA_000154225.1_ASM15422v1_genomic.fna.gz
12    GCA_000008005.1_ASM800v1_genomic.fna.gz
7    GCA_000012825.1_ASM1282v1_genomic.fna.gz
14    GCA_000016965.1_ASM1696v1_genomic.fna.gz
3    GCA_000008565.1_ASM856v1_genomic.fna.gz
4    GCA_000172575.2_ASM17257v2_genomic.fna.gz
7    GCA_000005845.2_ASM584v2_genomic.fna.gz
2    GCA_000008525.1_ASM852v1_genomic.fna.gz
6    GCA_000014425.1_ASM1442v1_genomic.fna.gz
6    GCA_000196035.1_ASM19603v1_genomic.fna.gz
4    GCA_000008805.1_ASM880v1_genomic.fna.gz
3    GCA_000008345.1_ASM834v1_genomic.fna.gz
4    GCA_000006765.1_ASM676v1_genomic.fna.gz
3    GCA_000012905.2_ASM1290v2_genomic.fna.gz
5    GCA_000017085.1_ASM1708v1_genomic.fna.gz
5    GCA_000007645.1_ASM764v1_genomic.fna.gz
7    GCA_000007265.1_ASM726v1_genomic.fna.gz
5    GCA_000007465.2_ASM746v2_genomic.fna.gz
4    GCA_000006885.1_ASM688v1_genomic.fna.gz
```


file opened in pandas and manipulated as follows to match the format of the 1520 genomes EDA approach:

Python REPL:
```
import pandas as pd

df = pd.read_csv('bei_16S_copy_numbers.txt',sep='\s+',header=None)
copy_no = round(100/df.iloc[:,0])
copy_no = [int(i) for i in copy_no]
df['copy_no'] = copy_no
df = df.drop(columns=0)
df.to_csv('sampled_files_key.txt',sep='\t',header=None,index=None)
```


```
spack load bwa@0.7.17
spack load samtools@1.10

cd /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/mock_community_ref_genomes/
for i in *fna ; do bwa index $i ; done

cd /pickett_flora/projects/read_simulation/analyses/snipen_2021_RMS/2021_07_26_test_mock_data/3_bwa_map_snipen_seqs

for query in ../2_kraken_snipen_seqs/SRR101997* ; do
  echo ${query#"../2_kraken_snipen_seqs/"}
  mkdir ${query#"../2_kraken_snipen_seqs/"}
  cd ${query#"../2_kraken_snipen_seqs/"}
  for genome in /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/mock_community_ref_genomes/*fna ; do
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    mkdir $genome_name
    cd $genome_name
    bwa mem $genome ../../../1_front_trim_snipen_seqs/trimmed_se.${query#"../2_kraken_snipen_seqs/"}_1.fastq \
    ../../../1_front_trim_snipen_seqs/trimmed_se.${query#"../2_kraken_snipen_seqs/"}_2.fastq > ${genome_name}.sam
    samtools view -b -f 0x2 ${genome_name}.sam > ${genome_name}.bam
    samtools sort -o ${genome_name}.sorted.bam ${genome_name}.bam
    samtools index ${genome_name}.sorted.bam
    samtools stats ${genome_name}.sorted.bam > ${genome_name}_stats.txt
    cd ..
  done
  cd ..
done
```

import pandas as pd
df = pd.read_csv('insert_size_average.txt',header=None,sep='\s+')
df[4].mean()

df = pd.read_csv('insert_size_sds.txt',header=None,sep='\s+')
df[5].mean()


```
cd SRR10199716/
for i in ./GCA* ; do cat ${i}/*stats.txt | grep "insert size average" >> insert_size_average.txt ; done
for i in ./GCA* ; do cat ${i}/*stats.txt | grep "insert size standard deviation" >> insert_size_sds.txt ; done
```

insert size average: 145.53499999999997  
insert size sd: 94.32499999999999

```
cd ../SRR10199724/
for i in ./GCA* ; do cat ${i}/*stats.txt | grep "insert size average" >> insert_size_average.txt ; done
for i in ./GCA* ; do cat ${i}/*stats.txt | grep "insert size standard deviation" >> insert_size_sds.txt ; done
```

insert size average: 145.07999999999998
insert size sd: 93.88500000000002

```
for i in ./GCA* ; do cat ${i}/*stats.txt | grep "insert size average" >> insert_size_average.txt ; done
for i in ./GCA* ; do cat ${i}/*stats.txt | grep "insert size standard deviation" >> insert_size_sds.txt ; done
```

insert size average: 143.64000000000004
insert size sd: 94.82999999999998

The average insert size across the three Snipen replicates was 145bp with 94bp standard deviation, so these values will be used for simulation.
I'll add 11bp (156bp total) for the adapter lengths that extend beyond the SBS start site as these were ligated before size selection and size selection occurred before adding the p5 and p7 indices. 148bp was used as the read length as the empirical data has R1 145bp and R2 151bp lengths.

Run bash script 'run_95_percent_digest_sim.sh' feeding in 'sampled_files_key.txt' as stdin1:

```
while read genome copies; do
  echo $genome $copies
  python3 /home/rkuster/readsynth/readsynth.py \
     -genome $genome \
     -m1 G/AATTC \
     -m2 T/TAA \
     -complete 0.95 \
     -o ./95_percent/ \
     -n $copies \
     -mean 156 \
     -sd 94 \
     -l 148 \
     -t 4 \
     -a1 /pickett_flora/projects/read_simulation/raw_data/adapters/snipen_ddRADseq/snipen_ddRADseq_adapters_R1.txt \
     -a2 /pickett_flora/projects/read_simulation/raw_data/adapters/snipen_ddRADseq/snipen_ddRADseq_adapters_R2.txt \
     -a1s 7 \
     -a2s 4 ;
done < $1
```