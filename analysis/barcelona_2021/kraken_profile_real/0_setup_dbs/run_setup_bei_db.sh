eval "$(conda shell.bash hook)"
conda activate kraken2

conda list > environment.log

cd /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/mock_community_ref_genomes

mkdir ../../kraken_dbs/snipen_bei_db

for i in *fna ; do
  kraken2-build --add-to-library $i --db ../../kraken_dbs/snipen_bei_db
done 

kraken2-build --download-taxonomy --db ../../kraken_dbs/snipen_bei_db
kraken2-build --build --db ../../kraken_dbs/snipen_bei_db

conda activate bracken
conda list > bracken_environment.log
bracken-build -d /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/snipen_bei_db -k 35 -l 150 -t 1
