## 2021.06.30
### select 30 genomes and randomly assign copy number to each, then simulate fastq files for each abundance profile provided in the key

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

### 2_simulate_fastq_files
copy or link the 30 chosen files to the local 'genomes' directory

use 'run_quantify_sim.sh' script to automate the fastq simulation of each of the 30 genomes (all files written to 'output' directory)

concatenate the resulting 30 fastq files into a single file:

```
cat ./output/*1.fastq > 2021_06_30_combined_R1.fastq
cat ./output/*2.fastq > 2021_06_30_combined_R2.fastq
```
