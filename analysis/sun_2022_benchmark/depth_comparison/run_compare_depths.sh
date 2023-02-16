eval "$(conda shell.bash hook)"
conda list > compare_environment.log

python3 ../../../code/scripts/compare_samtools_depths.py -d \
  real_1_combined_depths.tsv \
  sim_1_combined_depths.tsv -n \
  real_1 \
  sim_1 \
