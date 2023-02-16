eval "$(conda shell.bash hook)"
conda activate bracken
conda list > environment.log

for query in ../../extract_features_sim_mean_sd/1_cutadapt_trim/trimmed*1.fastq ; do
  query_dir=$(basename ${query%%_1.fastq})
  query_dir=$(basename ${query_dir#trimmed_se.})
  echo $query_dir
  query_name=${query%%_1.fastq}
  R1=../${query_name}_1.fastq
  R2=../${query_name}_2.fastq
  echo $R1
  echo $R2
  mkdir $query_dir
  cd $query_dir
  kraken2 --paired --threads 2 --report 2022_06_20.kreport --classified-out cseqs#.fq --unclassified-out useqs#.fq --db ../../../../raw_data/kraken_dbs/snipen_bei_db $R1 $R2 > 2022_06_20_run.txt 2> kraken2.log 
  bracken -d ../../../../raw_data/kraken_dbs/snipen_bei_db -i 2022_06_20.kreport -o 2022_06_20.bracken -r 150 -l S -t 10 2> bracken.log
  cd ..
done
