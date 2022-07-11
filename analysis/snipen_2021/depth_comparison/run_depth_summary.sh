eval "$(conda shell.bash hook)"
conda list > environment.log


combine_depths(){
  alignment_dir=$1
  combined_out=$2
  cd $alignment_dir
  echo "#CHROM  POS combined" > $combined_out
  for i in GCA* ; do
    sed -n "2,$"p ${i}/depth_${i}.txt >> $combined_out
  done
}

motif_combine_depths(){
  alignment_dir=$1
  combined_out=$2
  cd $alignment_dir
  echo "#CHROM  POS combined" > $combined_out
  for i in GCA* ; do
    sed -n "2,$"p ${i}/depth_${i}_motif_matching.txt >> $combined_out
  done
}

## bwa mem all mappings
#r1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_real/3_samtools/SRR10199716
#r1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/real_1_combined_depths.tsv
#combine_depths $r1 $r1_out
#
#r2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_real/3_samtools/SRR10199724
#r2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/real_2_combined_depths.tsv
#combine_depths $r2 $r2_out
#
#r3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_real/3_samtools/SRR10199725
#r3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/real_3_combined_depths.tsv
#combine_depths $r3 $r3_out
#
#
#
#s1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim/3_samtools/SRR10199716_sim
#s1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/sim_1_combined_depths.tsv
#combine_depths $s1 $s1_out
#
#s2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim/3_samtools/SRR10199724_sim
#s2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/sim_2_combined_depths.tsv
#combine_depths $s2 $s2_out
#
#s3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim/3_samtools/SRR10199725_sim
#s3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/sim_3_combined_depths.tsv
#combine_depths $s3 $s3_out
#
#
#
#smsd1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_mean_sd/3_samtools/SRR10199716_sim
#smsd1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/sim_mean_sd_1_combined_depths.tsv
#combine_depths $smsd1 $smsd1_out
#
#smsd2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_mean_sd/3_samtools/SRR10199724_sim
#smsd2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/sim_mean_sd_2_combined_depths.tsv
#combine_depths $smsd2 $smsd2_out
#
#smsd3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_mean_sd/3_samtools/SRR10199725_sim
#smsd3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/sim_mean_sd_3_combined_depths.tsv
#combine_depths $smsd3 $smsd3_out
#
#
## bwa mem matchings with motif sites only
#r1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_real/5_samtools_motif_only/SRR10199716
#r1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_real_1_combined_depths.tsv
#motif_combine_depths $r1 $r1_out
#
#r2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_real/5_samtools_motif_only/SRR10199724
#r2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_real_2_combined_depths.tsv
#motif_combine_depths $r2 $r2_out
#
#r3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_real/5_samtools_motif_only/SRR10199725
#r3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_real_3_combined_depths.tsv
#motif_combine_depths $r3 $r3_out
#
#
#
#s1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim/5_samtools_motif_only/SRR10199716_sim
#s1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_1_combined_depths.tsv
#motif_combine_depths $s1 $s1_out
#
#s2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim/5_samtools_motif_only/SRR10199724_sim
#s2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_2_combined_depths.tsv
#motif_combine_depths $s2 $s2_out
#
#s3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim/5_samtools_motif_only/SRR10199725_sim
#s3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_3_combined_depths.tsv
#motif_combine_depths $s3 $s3_out
#
#
#
#smsd1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_mean_sd/5_samtools_motif_only/SRR10199716_sim
#smsd1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_mean_sd_1_combined_depths.tsv
#motif_combine_depths $smsd1 $smsd1_out
#
#smsd2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_mean_sd/5_samtools_motif_only/SRR10199724_sim
#smsd2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_mean_sd_2_combined_depths.tsv
#motif_combine_depths $smsd2 $smsd2_out
#
#smsd3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_mean_sd/5_samtools_motif_only/SRR10199725_sim
#smsd3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_mean_sd_3_combined_depths.tsv
#motif_combine_depths $smsd3 $smsd3_out
#
#
#
#s100u1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_100_u/5_samtools_motif_only/SRR10199716_sim
#s100u1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_100u_1_combined_depths.tsv
#motif_combine_depths $s100u1 $s100u1_out
#
#s100u2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_100_u/5_samtools_motif_only/SRR10199724_sim
#s100u2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_100u_2_combined_depths.tsv
#motif_combine_depths $s100u2 $s100u2_out
#
#s100u3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_100_u/5_samtools_motif_only/SRR10199725_sim
#s100u3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_100u_3_combined_depths.tsv
#motif_combine_depths $s100u3 $s100u3_out
#
#
#
#s80c1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_80_c/5_samtools_motif_only/SRR10199716_sim
#s80c1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_80c_1_combined_depths.tsv
#motif_combine_depths $s80c1 $s80c1_out
#
#s80c2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_80_c/5_samtools_motif_only/SRR10199724_sim
#s80c2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_80c_2_combined_depths.tsv
#motif_combine_depths $s80c2 $s80c2_out
#
#s80c3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_80_c/5_samtools_motif_only/SRR10199725_sim
#s80c3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_80c_3_combined_depths.tsv
#motif_combine_depths $s80c3 $s80c3_out


s100c1=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_100_c/5_samtools_motif_only/SRR10199716_sim
s100c1_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_100c_1_combined_depths.tsv
motif_combine_depths $s100c1 $s100c1_out

s100c2=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_100_c/5_samtools_motif_only/SRR10199724_sim
s100c2_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_100c_2_combined_depths.tsv
motif_combine_depths $s100c2 $s100c2_out

s100c3=/pickett_flora/projects/read_simulation/analysis/snipen_2021/extract_features_sim_100_c/5_samtools_motif_only/SRR10199725_sim
s100c3_out=/pickett_flora/projects/read_simulation/analysis/snipen_2021/depth_comparison/mo_sim_100c_3_combined_depths.tsv
motif_combine_depths $s100c3 $s100c3_out


