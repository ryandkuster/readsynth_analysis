eval "$(conda shell.bash hook)"
conda activate bwa_and_samtools
conda list > environment.log

for query in ../1_cutadapt_trim/trimmed_se.*1.fastq ; do
  query_dir=$(basename ${query%%_R1.fastq})
  query_dir=$(basename ${query_dir#trimmed_se.})
  echo $query_dir
  query_name=${query%%R1.fastq}
  R1=../../${query_name}R1.fastq
  R2=../../${query_name}R2.fastq
  echo $R1
  echo $R2
  mkdir $query_dir
  cd $query_dir
  for genome in /pickett_flora/projects/metagenome_abundances/raw_data/liu_RMS/mock_community_estimate/SRR5298272_genomes_no_human_combined/*fna ; do
    echo $genome
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    mkdir $genome_name
    cd $genome_name

    # map reads and only keep correctly oriented pairs
    #bwa mem -t 4 $genome $R1 $R2 | samtools view -b -f 0x2 > ${genome_name}.bam
    bwa mem $genome $R1 $R2 | samtools view -b -f 0x2 > ${genome_name}.bam

    # create a sam version for custom python scripts
    samtools view ${genome_name}.bam > ${genome_name}.sam
    cd ..
  done
  cd ..
done
