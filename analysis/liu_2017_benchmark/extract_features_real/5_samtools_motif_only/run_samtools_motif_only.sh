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

    # keep only those reads with correct read orientation in pairs
    samtools view -b -f 0x2 ../../${query}/${genome_name}/${genome_name}.sam > ${genome_name}.bam
    
    # sort to run stats before fixmate
    samtools sort -o sort_${genome_name}.bam ${genome_name}.bam
    
    # collate is necessary to use fixmate
    samtools collate -o collate_${genome_name}.bam sort_${genome_name}.bam
    
    # run fixmate to correct for unexpected read pairings
    samtools fixmate -m collate_${genome_name}.bam fixmate_${genome_name}.bam
    samtools sort -o sort_fixmate_${genome_name}.bam fixmate_${genome_name}.bam
    
    # index for mapping in IGV
    samtools index sort_fixmate_${genome_name}.bam
    samtools stats sort_fixmate_${genome_name}.bam > stats_fixmate_${genome_name}.txt
    samtools view -h -o sort_fixmate_${genome_name}.sam sort_fixmate_${genome_name}.bam
    
    conda deactivate
    # use custom python script to narrow sam file to only those reads aligning perfectly to motif sites from simulation
    python3 /pickett_flora/projects/read_simulation/code/scripts/filter_sam_by_motif_sites.py \
        sort_fixmate_${genome_name}.sam \
        ../../../4_fragment_recreation/cut_sites/cut_sites_${genome_name}.fna.csv

    conda activate samtools
    # get the per-position depth
    samtools sort -o final_sort_fixmate_${genome_name}_motif_matching.sam sort_fixmate_${genome_name}_motif_matching.sam
    samtools depth -a -o depth_${genome_name}_motif_matching.txt -H final_sort_fixmate_${genome_name}_motif_matching.sam
    
    cd ..
  done
  cd ..
done

