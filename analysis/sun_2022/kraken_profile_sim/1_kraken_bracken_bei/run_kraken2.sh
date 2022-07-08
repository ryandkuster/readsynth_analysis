eval "$(conda shell.bash hook)"
conda activate bracken
conda list > environment.log

kraken2 --help

for query in ../../extract_features_sim/1_cutadapt_trim/combined_error_sim_metagenome.fastq ; do
  query_dir=$(basename ${query%%.fastq})
  echo $query_dir
  R1=../${query}
  echo $R1
  mkdir $query_dir
  cd $query_dir
  kraken2 --report 2022_07_01.kreport --classified-out cseqs1.fq --unclassified-out useqs1.fq --db /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/sun_atcc_db $R1 > 2022_07_01_run.txt 2> kraken2.log 
  bracken -d /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/sun_atcc_db -i 2022_07_01.kreport -o 2022_07_01.bracken -r 32 -l S -t 1 2> bracken.log
  cd ..
done
