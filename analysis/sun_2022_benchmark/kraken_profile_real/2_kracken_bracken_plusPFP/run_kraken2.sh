eval "$(conda shell.bash hook)"
conda activate bracken
conda list > environment.log

for query in ../../extract_features_real/1_cutadapt_trim/mock*fq ; do
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
  kraken2 --paired --threads 2 --report 2022_06_17.kreport --classified-out cseqs#.fq --unclassified-out useqs#.fq --db /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/k2_pluspfp_20220607 $R1 $R2 > 2022_06_17_run.txt 2> kraken2.log 
  bracken -d /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/k2_pluspfp_20220607 -i 2022_06_17.kreport -o 2022_06_17.bracken -r 150 -l S -t 10 2> bracken.log
  cd ..
done