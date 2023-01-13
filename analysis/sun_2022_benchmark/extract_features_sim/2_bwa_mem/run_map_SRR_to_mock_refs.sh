eval "$(conda shell.bash hook)"
conda activate bwa
conda list > environment.log

for query in ../1_cutadapt_trim/combined*.fastq ; do
  query_dir=$(basename ${query%%.fastq})
  echo $query_dir
  mkdir $query_dir
  cd $query_dir
  for genome in /pickett_flora/projects/read_simulation/raw_data/sun_2bRADM/mock_community_ref_genomes/atcc_msa_1002/*fasta ; do
    echo $genome
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fasta}
    mkdir $genome_name
    cd $genome_name
    bwa mem $genome ../../${query} > ${genome_name}.sam
    cd ..
  done
  cd ..
done
