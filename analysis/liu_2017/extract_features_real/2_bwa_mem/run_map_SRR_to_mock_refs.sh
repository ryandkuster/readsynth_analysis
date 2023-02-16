eval "$(conda shell.bash hook)"
conda activate bwa_and_samtools
conda list > environment.log

for query in ../1_cutadapt_trim/adapted*1.fastq ; do
  query_dir=$(basename ${query%%_R1.fastq})
  query_dir=$(basename ${query_dir#adapted_})
  echo $query_dir
  query_name=${query%%R1.fastq}
  R1=../../${query_name}R1.fastq
  R2=../../${query_name}R2.fastq
  echo $R1
  echo $R2
  mkdir $query_dir
  cd $query_dir
  for genome in ../../../../raw_data/liu_RMS/mock_community_estimate/SRR5298272_genomes/*fna ; do
    echo $genome
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    mkdir $genome_name
    cd $genome_name
    bwa mem $genome $R1 $R2 | samtools view -b -f 0x2 > ${genome_name}.bam
    cd ..
  done
  cd ..
done
