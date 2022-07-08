eval "$(conda shell.bash hook)"
conda activate kraken2

conda list > environment.log

cd /pickett_flora/projects/read_simulation/raw_data/sun_2bRADM/mock_community_ref_genomes/manually_renamed_atcc_msa_1002

mkdir /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/sun_atcc_db

echo "adding to library"
for i in *fasta ; do
  kraken2-build --add-to-library $i --db /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/sun_atcc_db
done 

echo "downloading tax"
kraken2-build --download-taxonomy --db /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/sun_atcc_db

echo "building db"
kraken2-build --build --db /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/sun_atcc_db --kmer-len 31 --minimizer-len 30

conda activate bracken
conda list > bracken_environment.log
bracken-build -d /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/sun_atcc_db -k 31 -l 32 -t 1
