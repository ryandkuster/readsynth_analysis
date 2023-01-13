eval "$(conda shell.bash hook)"
conda activate cutadapt
conda list &> environment.log

ln -s ../../simulate_mean_sd/SRR10199716_output/error_sim_metagenome_R1.fastq ./SRR10199716_sim_1.fastq
ln -s ../../simulate_mean_sd/SRR10199716_output/error_sim_metagenome_R2.fastq ./SRR10199716_sim_2.fastq
ln -s ../../simulate_mean_sd/SRR10199724_output/error_sim_metagenome_R1.fastq ./SRR10199724_sim_1.fastq
ln -s ../../simulate_mean_sd/SRR10199724_output/error_sim_metagenome_R2.fastq ./SRR10199724_sim_2.fastq
ln -s ../../simulate_mean_sd/SRR10199725_output/error_sim_metagenome_R1.fastq ./SRR10199725_sim_1.fastq
ln -s ../../simulate_mean_sd/SRR10199725_output/error_sim_metagenome_R2.fastq ./SRR10199725_sim_2.fastq

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
