eval "$(conda shell.bash hook)"
conda list > compare_environment.log

python3 /pickett_flora/projects/read_simulation/code/scripts/compare_samtools_depths.py -d \
  mo_real_1_combined_depths.tsv \
  mo_real_2_combined_depths.tsv \
  mo_real_3_combined_depths.tsv \
  mo_sim_1_combined_depths.tsv \
  mo_sim_2_combined_depths.tsv \
  mo_sim_3_combined_depths.tsv \
  mo_sim_mean_sd_1_combined_depths.tsv \
  mo_sim_mean_sd_2_combined_depths.tsv \
  mo_sim_mean_sd_3_combined_depths.tsv \
  real_1_combined_depths.tsv \
  real_2_combined_depths.tsv \
  real_3_combined_depths.tsv \
  sim_1_combined_depths.tsv \
  sim_2_combined_depths.tsv \
  sim_3_combined_depths.tsv \
  sim_mean_sd_1_combined_depths.tsv \
  sim_mean_sd_2_combined_depths.tsv \
  sim_mean_sd_3_combined_depths.tsv -n \
  mo_real_1 \
  mo_real_2 \
  mo_real_3 \
  mo_sim_1 \
  mo_sim_2 \
  mo_sim_3 \
  mo_sim_mean_sd_1 \
  mo_sim_mean_sd_2 \
  mo_sim_mean_sd_3 \
  real_1 \
  real_2 \
  real_3 \
  sim_1 \
  sim_2 \
  sim_3 \
  sim_mean_sd_1 \
  sim_mean_sd_2 \
  sim_mean_sd_3
