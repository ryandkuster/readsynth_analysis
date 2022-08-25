# **read simulation**

# **directory structure**
- analyses/ **(all specific objectives for processing data)**
  - snipen_2021
  - snipen_2021_update
  - sun_2022
  - liu_2017
- code/ **(specific software used for analyses)**
  - archive
  - scripts
- misc/
  - visuals/ **(visuals presented in this markdown file)**
- notebooks
  - compare_depths
  - compare_kraken
  - cut_probability
  - error
  - fragment_summaries
  - performance
- raw_data/ **(unprocessed data)**
  - adapters/
  - archive/
  - barcelona_16S/
  - kraken_dbs/
  - large_metagenome/
  - liu_RMS/
  - ravi_16S/
  - snipen_RMS/
  - sun_2bRADM/
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
