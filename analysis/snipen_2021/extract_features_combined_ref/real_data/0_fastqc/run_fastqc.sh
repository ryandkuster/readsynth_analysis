eval "$(conda shell.bash hook)"
conda activate fastqc

for fastq in ../../../../../raw_data/snipen_RMS/SRR*fastq ; do
  fastqc $fastq -o ./
done
