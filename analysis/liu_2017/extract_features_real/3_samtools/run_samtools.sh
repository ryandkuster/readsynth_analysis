eval "$(conda shell.bash hook)"
conda activate samtools
conda list > environment.log

for query in ../2_bwa_mem/trimmed_se.* ; do
  query_dir=$(basename $query)
  echo $query_dir
  mkdir $query_dir
  cd $query_dir
  for genome in ../../../../raw_data/liu_RMS/mock_community_estimate/SRR5298272_genomes/*fna ; do
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    mkdir $genome_name
    cd $genome_name

    # keep only those reads with correct read orientation in pairs
    #samtools view -b -f 0x2 ../../${query}/${genome_name}/${genome_name}.sam > ${genome_name}.bam
    
    # sort to run stats before fixmate
    samtools sort -o sort_${genome_name}.bam ../../${query}/${genome_name}/${genome_name}.bam
    
    # collate is necessary to use fixmate
    samtools collate -o collate_${genome_name}.bam sort_${genome_name}.bam
    
    # run fixmate to correct for unexpected read pairings
    samtools fixmate -m collate_${genome_name}.bam fixmate_${genome_name}.bam
    samtools sort -o sort_fixmate_${genome_name}.bam fixmate_${genome_name}.bam
    
    # convert to sam for custom scripts to motif filter
    samtools view -h -o sort_fixmate_${genome_name}.sam sort_fixmate_${genome_name}.bam
    
    cd ..
  done
  cd ..
done

