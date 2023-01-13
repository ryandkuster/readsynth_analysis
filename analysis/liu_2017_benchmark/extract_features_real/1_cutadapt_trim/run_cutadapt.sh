eval "$(conda shell.bash hook)"
conda activate cutadapt
conda list &> environment.log

cutadapt -u 2 -o trimmed_se.SRR5298272_R1.fastq SRR5298272_R1.fastq 
cp SRR5298272_R2.fastq trimmed_se.SRR5298272_R2.fastq

cp SRR5298274_R1.fastq trimmed_se.SRR5298274_R1.fastq
cp SRR5298274_R2.fastq trimmed_se.SRR5298274_R2.fastq

cutadapt -u 3 -o trimmed_se.SRR5360684_R1.fastq SRR5360684_R1.fastq
cutadapt -u 1 -o trimmed_se.SRR5360684_R2.fastq SRR5360684_R2.fastq
