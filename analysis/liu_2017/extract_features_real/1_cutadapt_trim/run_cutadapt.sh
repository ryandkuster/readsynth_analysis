eval "$(conda shell.bash hook)"
conda activate cutadapt
conda list &> environment.log

cutadapt -u 2 -o trimmed_se.SRR5298272_R1.fastq SRR5298272_R1.fastq 
cp SRR5298272_R2.fastq trimmed_se.SRR5298272_R2.fastq

for fastq in trimmed*1.fastq
do
  cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTCTCNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGCAGCGAGAGTGTGTAGATCTCGGTGGTCGCCGTATCATT \
  -o adapted_${fastq} \
  -p adapted_${fastq%%1.fastq}2.fastq \
  $fastq \
  ${fastq%%1.fastq}2.fastq \
  &> cutadapt_${fastq%%_R1.fastq}.log
done 
