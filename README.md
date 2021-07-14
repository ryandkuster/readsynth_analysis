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

There.
*<sub><sup>/pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/1_profile_reads</sub></sup>*

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

*<sub><sup>/pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2_simulate_f_prausnitzii</sub></sup>*

Download F. prausnitzii

```
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/312/465/GCF_003312465.1_ASM331246v1/GCF_003312465.1_ASM331246v1_genomic.fna.gz"
```

Run readsynth.py (commit 4df8ce1f90394987637b4a83443ec53e808c1af2)

```
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2_simulate_f_prausnitzii
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

*<sub><sup>/pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/3_compare_rms_to_sim_mapping</sub></sup>*

Grab the reads for each of the three Liu datasets that Kraken identified as F. prausnitzii.

```
python3 /pickett_flora/projects/read_simulation/code/scripts/pull_kraken_reads.py ../../1_profile_reads/SRR5298272/cseqs_1.fq 853
python3 /pickett_flora/projects/read_simulation/code/scripts/pull_kraken_reads.py ../../1_profile_reads/SRR5298272/cseqs_2.fq 853
```

To BWA map the reads without error due to fastq naming discrepency, rename the R2 headers to be identical with R1:

```
python3 /pickett_flora/projects/read_simulation/code/scripts/fix_R2_headers.py taxid_853_cseqs_2.fq
```

BWA map reads from RMS project and reads simulated from readsynth.py to the F. prausnitzii genome.

```
spack load bwa@0.7.17
spack load samtools@1.10
cp ../2_simulate_f_prausnitzii/GCF_003312465.1_ASM331246v1_genomic.fna ./
bwa index GCF_003312465.1_ASM331246v1_genomic.fna

cd SRR5298272
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna taxid_853_cseqs_1.fq modified_taxid_853_cseqs_2.fq > SRR5298272.sam
samtools view -b SRR5298272.sam > SRR5298272.bam
samtools sort -o SRR5298272.sorted.bam SRR5298272.bam
samtools index SRR5298272.sorted.bam
samtools depth -o SRR5298272.depth.txt -H SRR5298272.sorted.bam
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py SRR5298272.depth.txt
12.762109703977144
cd ..

cd SRR5298274
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna taxid_853_cseqs_1.fq modified_taxid_853_cseqs_2.fq > SRR5298274.sam
samtools view -b SRR5298274.sam > SRR5298274.bam
samtools sort -o SRR5298274.sorted.bam SRR5298274.bam
samtools index SRR5298274.sorted.bam
samtools depth -o SRR5298274.depth.txt -H SRR5298274.sorted.bam
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py SRR5298274.depth.txt
4.958139806567638
cd ..

cd SRR5360684
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna taxid_853_cseqs_1.fq modified_taxid_853_cseqs_2.fq > SRR5360684.sam
samtools view -b SRR5360684.sam > SRR5360684.bam
samtools sort -o SRR5360684.sorted.bam SRR5360684.bam
samtools index SRR5360684.sorted.bam
samtools depth -o SRR5360684.depth.txt -H SRR5360684.sorted.bam
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py SRR5360684.depth.txt
9.023378056269712
cd ..

mkdir simulated_F_prausnitzii
cd simulated_F_prausnitzii
cp ../../2_simulate_f_prausnitzii/output/GCF_003312465.1_ASM331246v1_genomic.fna_R\* ./
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna GCF_003312465.1_ASM331246v1_genomic.fna_R1.fastq GCF_003312465.1_ASM331246v1_genomic.fna_R2.fastq > simulated.sam
samtools view -b simulated.sam > simulated.bam                                                                
samtools sort -o simulated.sorted.bam simulated.bam                                                           
samtools index simulated.sorted.bam 
samtools depth -o simulated.depth.txt -H simulated.sorted.bam
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py simulated.depth.txt
1.3841179998500637
```

## 2021.07.13
---
## **repeat simulation using error profile and n = 13 from SRR5298272 results and incomplete digest**


Run readsynth.py (commit 4cf90a729ccd292c4c78bcfff27f3b546bcfcd50)

```
cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/2_simulate_f_prausnitzii
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

cd /pickett_flora/projects/read_simulation/analyses/liu_2017_RMS/3_compare_rms_to_sim_mapping/
mkdir simulated_F_prausnitzii_2
cd simulated_F_prausnitzii_2
cp ../../2_simulate_f_prausnitzii/output2/GCF_003312465.1_ASM331246v1_genomic.fna_R\* ./
bwa mem ../GCF_003312465.1_ASM331246v1_genomic.fna GCF_003312465.1_ASM331246v1_genomic.fna_R1.fastq GCF_003312465.1_ASM331246v1_genomic.fna_R2.fastq > simulated2.sam
samtools view -b simulated2.sam > simulated2.bam                                                                
samtools sort -o simulated2.sorted.bam simulated2.bam                                                           
samtools index simulated2.sorted.bam 
samtools depth -o simulated2.depth.txt -H simulated2.sorted.bam
python3 /pickett_flora/projects/read_simulation/code/scripts/get_samtools_avg_depth.py simulated2.depth.txt
3.6272405227291276
```




