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

![correlation plot of all variables](https://github.com/ryandkuster/read_simulation/misc/visuals/correlation.png)

After normalizing the fragment counts per genome (count/genome(bp)), the most notable relationship is the negative correlation of gc_ratio and framents in the 1 to 200 bp range:

![scatter plot of gc and 1 - 200 bp fragments](https://github.com/ryandkuster/read_simulation/misc/visuals/gc_ratio_and_fragments_plot.png)


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

![sampling from complete digest at 1 X](https://github.com/ryandkuster/read_simulation/misc/visuals/hist_sampled_half_chr_1_complete.png)

Similarly, using an incomplete digest will produce an overabundance of fragments from some parts of the genome, and a weighting system must account for this. The figure below details an incomplete digest distribution (blue) with the weight-adjusted sampling distribution approximating 1x representation of the original genome (note the count similarity to the complete digest above).

![sampling from incomplete digest at 1 X](https://github.com/ryandkuster/read_simulation/misc/visuals/hist_sampled_half_chr_1_incomplete.png)

Importantly, simulating complete EcoR1/MseI digests where R1/R2 adapters ligate in a 5' to 3' specific manner requires weight consideration to accurately approximate genome copy number. Complete digests with orientation won't encounter overlap at any point on either template or non-template strands.

An example from today's simulation from the GCA_003468235 accession using 4 as the original genome copy number:

![4 X coverage example](https://github.com/ryandkuster/read_simulation/misc/visuals/GCA_003468235_4X.png)

Finally, a genome copy number of 20 for accession GCA_003433755:

![20 X coverage example](https://github.com/ryandkuster/read_simulation/misc/visuals/GCA_003433755_20X.png)

## 2021.07.06
---
## **community profiling using Kraken2**