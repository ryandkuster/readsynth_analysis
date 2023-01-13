eval "$(conda shell.bash hook)"
conda activate samtools
conda list > environment.log

for query in ../2_bwa_mem/combined* ; do
  query_dir=$(basename $query)
  echo $query_dir
  mkdir $query_dir
  cd $query_dir
  for genome in /pickett_flora/projects/read_simulation/raw_data/sun_2bRADM/mock_community_ref_genomes/atcc_msa_1002/*fasta ; do
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fasta}
    mkdir $genome_name
    cd $genome_name

    # drop unmapped reads, convert to bam
    samtools view -b -F 0x4 ../../${query}/${genome_name}/${genome_name}.sam > ${genome_name}.bam
    
    # sort to run stats
    samtools sort -o sort_${genome_name}.bam ${genome_name}.bam
    samtools stats sort_${genome_name}.bam > stats_${genome_name}.txt
    
    # index for mapping in IGV
    samtools index sort_${genome_name}.bam
    
    # get the per-position depth
    samtools depth -a -o depth_${genome_name}.txt -H sort_${genome_name}.bam

    cd ..
  done
  cd ..
done

