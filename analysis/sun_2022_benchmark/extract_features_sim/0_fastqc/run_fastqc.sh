eval "$(conda shell.bash hook)"
conda activate fastqc

for fastq in ../../simulate/output_mock-1ng.BcgI/combined*.fastq ; do
  fastqc $fastq -o ./
done
