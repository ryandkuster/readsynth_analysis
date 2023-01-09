eval "$(conda shell.bash hook)"
conda activate samtools
conda list > environment.log

for query in ../2_bwa_mem/SRR* ; do
  query_dir=$(basename $query)
  echo $query_dir
  mkdir $query_dir
  cd $query_dir
  for genome in /pickett_flora/projects/read_simulation/raw_data/snipen_RMS/mock_community_ref_genomes/*fna ; do
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    mkdir $genome_name
    cd $genome_name

    # sort to run stats before fixmate
    samtools sort -o sort_${genome_name}.bam ../../${query}/${genome_name}/${genome_name}.bam
    samtools stats sort_${genome_name}.bam > stats_${genome_name}.txt
    
    # collate is necessary to use fixmate
    samtools collate -o collate_${genome_name}.bam sort_${genome_name}.bam
    
    # run fixmate to correct for unexpected read pairings
    samtools fixmate -m collate_${genome_name}.bam fixmate_${genome_name}.bam
    samtools sort -o sort_fixmate_${genome_name}.bam fixmate_${genome_name}.bam
    
    # index for mapping in IGV
    samtools index sort_fixmate_${genome_name}.bam
    samtools stats sort_fixmate_${genome_name}.bam > stats_fixmate_${genome_name}.txt
    samtools view -h -o sort_fixmate_${genome_name}.sam sort_fixmate_${genome_name}.bam
    
    # get the per-position depth
    samtools depth -a -o depth_${genome_name}.txt -H sort_fixmate_${genome_name}.bam
    
    cd ..
  done
  cd ..
done

