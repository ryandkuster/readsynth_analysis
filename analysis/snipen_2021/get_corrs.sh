for i in extract*real/3_samtools/SRR* ; do
  fastq=$(basename $i)
  for j in ${i}/GCA* ; do
    genome=$(basename $j)
    real=extract_features_real/3_samtools/${fastq}/${genome}/depth_${genome}.txt
    sim=extract_features_sim_mean_sd/3_samtools/${fastq}_sim/${genome}/depth_${genome}.txt
    echo $real
    echo $sim
    python3 /pickett_flora/projects/read_simulation/code/scripts/compare_samtools_depths.py -d $real $sim -n real sim
  done
done
