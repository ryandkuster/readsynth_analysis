eval "$(conda shell.bash hook)"
conda activate bracken
conda list > environment.log

for query in ../../simulate_no_human/SRR5298272_no_human_10M_output/sim_metagenome_*1.fastq ; do
  query_dir=$(basename ${query%%_R1.fastq})
  echo $query_dir
  query_name=${query%%_R1.fastq}
  R1=../${query_name}_R1.fastq
  R2=../${query_name}_R2.fastq
  echo $R1
  echo $R2
  mkdir $query_dir
  cd $query_dir
  kraken2 --paired --threads 6 --report 2022_10_04.kreport --classified-out cseqs#.fq --unclassified-out useqs#.fq --db ../../../../raw_data/kraken_dbs/k2_pluspfp_20220607 $R1 $R2 > 2022_10_04_run.txt 2> kraken2.log 
  bracken -d ../../../../raw_data/kraken_dbs/k2_pluspfp_20220607 -i 2022_10_04.kreport -o 2022_10_04.bracken -r 75 -l S -t 6 2> bracken.log
  cd ..
done
