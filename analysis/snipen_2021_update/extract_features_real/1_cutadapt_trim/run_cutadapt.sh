eval "$(conda shell.bash hook)"
conda activate cutadapt
conda list &> environment.log

ln -s ../../../../raw_data/snipen_RMS/*[12].fastq ./

for fastq in *16_1.fastq
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

for fastq in *24_1.fastq
do
  cutadapt \
  -a TTACTCAGGACTCATCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG \
  -A GAATTGGTACGCAGTCATCAGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  -o adapted_${fastq} \
  -p adapted_${fastq%%_1.fastq}_2.fastq \
  $fastq \
  ${fastq%%_1.fastq}_2.fastq \
  &> cutadapt_${fastq%%_1.fastq}.log
done 

for fastq in *25_1.fastq
do
  cutadapt \
  -a TTACTCAGGACTCATCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG \
  -A GAATTGGTACGCAGTCGCTACCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  -o adapted_${fastq} \
  -p adapted_${fastq%%_1.fastq}_2.fastq \
  $fastq \
  ${fastq%%_1.fastq}_2.fastq \
  &> cutadapt_${fastq%%_1.fastq}.log
done 

for fastq in adapted*1.fastq ; do cutadapt -u 11 -o trimmed_se.${fastq#adapted_} $fastq ; done 
for fastq in adapted*2.fastq ; do cutadapt -u 13 -o trimmed_se.${fastq#adapted_} $fastq ; done 
