eval "$(conda shell.bash hook)"
conda activate bwa
conda list > environment.log

for query in ../1_cutadapt_trim/trimmed*1.fastq ; do
  query_dir=$(basename ${query%%_1.fastq})
  query_dir=$(basename ${query_dir#trimmed_se.})
  echo $query_dir
  query_name=${query%%_1.fastq}
  R1=../../${query_name}_1.fastq
  R2=../../${query_name}_2.fastq
  echo $R1
  echo $R2
  mkdir $query_dir
  cd $query_dir
  for genome in ../../../../raw_data/snipen_RMS/mock_community_ref_genomes/*fna ; do
    echo $genome
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    mkdir $genome_name
    cd $genome_name
    bwa mem $genome $R1 $R2 > ${genome_name}.sam
    cd ..
  done
  cd ..
done
