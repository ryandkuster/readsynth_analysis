# conda activate cutadapt
conda list &> environment.log

for fastq in *_1.fastq
do
  cutadapt \
  -a CTCAGGACTCATCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG \
  -A GGTACGCAGTCCGTACGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  -o adapted_${fastq} \
  -p adapted_${fastq%%_1.fastq}_2.fastq \
  $fastq \
  ${fastq%%_1.fastq}_2.fastq \
  &> cutadapt_${fastq%%_1.fastq}.log
done 

for fastq in adapted*1.fastq ; do cutadapt -u 11 -o trimmed_se.${fastq#adapted_} $fastq ; done 
for fastq in adapted*2.fastq ; do cutadapt -u 13 -o trimmed_se.${fastq#adapted_} $fastq ; done 