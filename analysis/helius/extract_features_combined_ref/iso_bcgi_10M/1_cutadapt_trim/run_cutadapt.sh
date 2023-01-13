sim_dir=$1
echo $sim_dir

eval "$(conda shell.bash hook)"

ln -s ${sim_dir}sim*fastq ./

conda activate cutadapt
conda list &> environment.log

for fastq in *1.fastq
do
  cutadapt -u 2 -o trimmed_se.${fastq} ${fastq}
  cutadapt -u 2 -o trimmed_se.${fastq%%1.fastq}2.fastq ${fastq%%1.fastq}2.fastq
done 
