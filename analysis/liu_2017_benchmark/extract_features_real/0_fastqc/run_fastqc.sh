eval "$(conda shell.bash hook)"
conda activate fastqc
conda list > environment.log

for fastq in ../../../../raw_data/liu_RMS/*[12]*fastq ; do
  fastqc $fastq -o ./
done
