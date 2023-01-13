# **readsynth_analysis**

# **directory structure**
- analysis/ **(all specific objectives for processing data)**
  - helius
  - liu_2017
  - liu_2017_benchmark
  - snipen_2021
  - snipen_2021_benchmark
  - sun_2022_benchmark

- code/ **(specific scripts used for analyses)**
  - scripts

- misc/
  - visuals/ **(visuals presented in this markdown file)**

- notebooks
  - compare_depths
  - compare_kraken
  - cut_probability
  - error
  - extract_read_features
  - fragment_summaries
  - performance
  - plot_helius_results
  - seq_distances

- raw_data/ **(unprocessed data)**
  - adapters
  - helius
  - kraken_dbs
  - liu_RMS
  - snipen_RMS
  - sun_2bRADM

- README.md **(you're reading me)**

# **raw_data**

## download Snipen 2021 data

```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10199716/SRR10199716
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10199724/SRR10199724
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10199725/SRR10199725
spack load sratoolkit@2.10.7

for i in SRR*; do fastq-dump --split-files $i ; done
```

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

The header naming won't work with bwa mem: 

```
for i in *R2.fastq ; do python3 ../../code/scripts/fix_R2_headers2.py $i && mv modified_${i} $i ; done
```

# analysis

## benchmarking analyses (directories ending in 'benchmark')

For liu_2017_benchmark, sun_2022_benchmark, and snipen_2021_benchmark directories a similar processing pipeline was performed to simulate existing sequence libraries and assess how closely these simulated datasets resembled their corresponding real sequence data.

The general order of analysis is as follows:
0.  kraken_profile_real - if no mock community standard, use kraken2/bracken to pull putative taxa for purposes of simulation
1.  extract_features_real - get the fragment lengths and digest rates from real sequence data
2.  simulate_x - based on the extracted sequence features, simulate library using x parameters
3.  extract_features_sim_x - get the fragment lenghts and read depths of the simulated sequence data
4.  depth_comparison - pull the samtools depth information from all real and simulated datasets and perform correlation analysis 


## profiling and fixed length taxonomic ratio analyses

These analyses follow a nearly identical process as the benchmarking procedures, but BWA-MEM map to concatenated reference genomes to avoid all instances of multimapping. Instead of calculating read depth per position as in the benchmarking above, fragment counts are used to estimate the relative abundances of the taxa profiled.
