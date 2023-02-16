eval "$(conda shell.bash hook)"
conda activate samtools
conda list > environment.log

for query in ../2_bwa_mem/SRR* ; do
  query_dir=$(basename $query)
  echo $query_dir
  mkdir $query_dir
  cd $query_dir
  for genome in ../../../../raw_data/snipen_RMS/mock_community_ref_genomes/*fna ; do
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    mkdir $genome_name
    cd $genome_name

    # keep only those reads with correct read orientation in pairs
    samtools view -b -f 0x2 ../../${query}/${genome_name}/${genome_name}.sam > ${genome_name}.bam
    
    # sort to run stats before fixmate
    samtools sort -o sort_${genome_name}.bam ${genome_name}.bam
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
    
    echo -e "query\tflag\tref\tstart\tend\tedit_dist\tread_length\tlength" > summary_${genome_name}.tsv 
    samtools view sort_fixmate_${genome_name}.sam | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $4+length($10)-1 "\t" $12 "\t" length($10) "\t" $9}' >> summary_${genome_name}.tsv

    cd ..
  done
  cd ..
done

