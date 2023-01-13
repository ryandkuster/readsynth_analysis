sim_dir=$1
echo $sim_dir

eval "$(conda shell.bash hook)"

#ln -s ${sim_dir}sim*fastq ./
ln -s ${sim_dir}error*fastq ./

conda activate cutadapt
conda list &> environment.log

for fastq in *1.fastq
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

cutadapt -u 4 -o trimmed_se.${fastq} adapted_${fastq}
cp adapted_${fastq%%1.fastq}2.fastq trimmed_se.${fastq%%1.fastq}2.fastq
