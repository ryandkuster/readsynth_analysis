eval "$(conda shell.bash hook)"
conda list > environment.log


combine_depths(){
  alignment_dir=$1
  combined_out=$2
  cd $alignment_dir
  echo "#CHROM  POS combined" > $combined_out
  for i in [A-Z]*[0-9] ; do
    sed -n "2,$"p ${i}/depth_${i}.txt >> $combined_out
  done
}

# bwa mem all mappings
r1=../extract_features_real/3_samtools/mock-1ng.BcgI/
r1_out=../depth_comparison/real_1_combined_depths.tsv
combine_depths $r1 $r1_out

s1=../extract_features_sim/3_samtools/combined_error_sim_metagenome/
s1_out=../depth_comparison/sim_1_combined_depths.tsv
combine_depths $s1 $s1_out

