eval "$(conda shell.bash hook)"
conda activate fastqc
conda list > environment.log

for fastq in /pickett_flora/projects/read_simulation/raw_data/liu_RMS/SRR5298272*[12]*fastq ; do
  fastqc $fastq -o ./
done
