*note: this file corresponds to the statonlabprivate wiki titled 'read_simulation'*

## directory structure
- README.md **(you're reading me)**
- analyses/ **(holds all specific objectives for processing data)**
- code/ **(holds specific software used for analyses)**
- raw_data/ **(unprocessed data)**
  - 1520_genomes/
- lab_notebook.md **(daily log of activities)**

## data

### 1520_genomes

https://doi.org/10.1038/s41587-018-0008-8

This dataset is relevant for simulating reads of similar sequencing preparation (human gut microbiota).

## 2021.06.14
### automate download of all accessions

Installed the 1,520 genome accessions produced in https://doi.org/10.1038/s41587-018-0008-8; wget called on the list of accession numbers. 


## 2021.06.16
### summarize 1520 accessions

Automatically download all accessions and provide summary statistics on the assembly.

Notes: The necessity to produce summary statistics for these assemblies for later analyses gave way to creating custom python scripts to pull the data. The custom scripts are stored in the /raw_data/genomes/1520_genomes directory with documentation.

## 2021.06.28

### analyses/1520_genomes_RAD_characteristics/2021_06_28_test_1520_genomes/1_RE_digest_genomes
#### in silico digest each of the 1,520 genomes using MseI and EcoRI motifs using 'readsynth' scripts.

readsynth.py from commit 1e72bf20d58ef8b7af03edc661cca810a3598b85

For the purposes of exploratory analysis, only RE digestion (complete) will be performed using complete digestion. The MseI and EcoRI have no issue of overlapping recognition sites, but this needs to be addressed for future simulations, as duplicate loci will be represented in these events.

The output of each simulation will consist of 4 files:
- raw_digest...csv (the absolute, raw fragments produced by the simulated RE digest) 
- hist_raw_digest...pdf (histogram of the raw fragments)
- sampled...csv (the absolute, raw fragments produced by the simulated RE digest)
- hist_sampled (histogram of the sampled fragments)

*note* The raw_digest...csv files will contain all possible cut sites considering both template and reverse-complemented template orientations. The originating strand is represented in the 'strand' column using a '+' (template) or '-' (reverse). Therefore, if 10,000 fragments are produced, we can assume that there are truly only 5,0000 fragments, as long as digests are complete and there is no strict orientation of cut sites on the 5' or 3' end. If there is an orientation requirement (e.g. adapter ligation requires that MseI be only on the 5' end and EcoRI can only be on the 3' end) then it makes sense to store reads in this strand-specific manner.

By default the maximum fragment length is set to mean + 6sd, so for these test, using mean of 400 and standard deviation of 100, we can expect the maximum fragment length to be included to max out at 1000bp.


## 2021.06.29

### analyses/1520_genomes_RAD_characteristics/2021_06_28_test_1520_genomes/2_pull_digest_stats
#### for every template strand ('+') from the raw_digest...csv file, pull the total number of fragments as well as ranges of fragment sizes (e.g. 1 to 200bp, 201 to 400bp, 401 to 600bp).

## 2021.06.30
### analyses/1520_genomes_RAD_characteristics/2021_06_30_quantify_30_genomes/1_create_genome_keys
#### select 30 genomes and randomly assign copy number to each, then simulate fastq files for each abundance profile provided in the key

### 1_create_genome_keys
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

### analyses/1520_genomes_RAD_characteristics/2021_06_30_quantify_30_genomes/2_simulate_fastq_files
#### 2_simulate_fastq_files
copy or link the 30 chosen files to the local 'genomes' directory

use 'run_quantify_sim.sh' script to automate the fastq simulation of each of the 30 genomes (all files written to 'output' directory)

using readsynth.py commit 4df8ce1f90394987637b4a83443ec53e808c1af2

concatenate the resulting 30 fastq files into a single file:

```
cat ./output/*1.fastq > 2021_06_30_combined_R1.fastq
cat ./output/*2.fastq > 2021_06_30_combined_R2.fastq
```
