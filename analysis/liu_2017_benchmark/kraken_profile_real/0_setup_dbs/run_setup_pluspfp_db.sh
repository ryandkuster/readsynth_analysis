eval "$(conda shell.bash hook)"
conda activate bracken
conda list > bracken_environment.log
bracken-build -d /pickett_flora/projects/read_simulation/raw_data/kraken_dbs/k2_pluspfp_20220607 -k 35 -l 74 -t 1