eval "$(conda shell.bash hook)"
conda activate samtools

for query in ../3_samtools/SRR* ; do
  query_dir=$(basename $query)
  echo $query_dir
  mkdir $query_dir
  cd $query_dir
  for genome in /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/mock_community_ref_genomes/*fna ; do
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    mkdir $genome_name
    cd $genome_name
    # get the per-position depth using all postitions (-a) using only the first read for overlaps (-s)
    samtools sort -o final_sort_fixmate_${genome_name}_motif_matching.sam ../../${query}/${genome_name}/motif_only_sort_fixmate_${genome_name}.sam
    samtools depth -a -s -o depth_${genome_name}_motif_matching.txt -H final_sort_fixmate_${genome_name}_motif_matching.sam
    samtools view -b final_sort_fixmate_${genome_name}_motif_matching.sam > final_sort_fixmate_${genome_name}_motif_matching.bam
    samtools index final_sort_fixmate_${genome_name}_motif_matching.bam

    cd ..
  done
  cd ..
done
