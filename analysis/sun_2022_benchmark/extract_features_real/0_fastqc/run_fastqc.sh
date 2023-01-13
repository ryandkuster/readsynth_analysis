eval "$(conda shell.bash hook)"
conda activate fastqc

for fastq in /pickett_flora/projects/read_simulation/raw_data/sun_2bRADM/2brad_msa/*fq ; do
  fastqc $fastq -o ./
done
