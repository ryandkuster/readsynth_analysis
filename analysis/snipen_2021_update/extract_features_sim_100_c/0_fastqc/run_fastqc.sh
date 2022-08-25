eval "$(conda shell.bash hook)"
conda activate fastqc
conda list > environment.log

ln -s ../../simulate_100_c/SRR10199716_output/error_sim_metagenome_R1.fastq ./SRR10199716_sim_1.fastq
ln -s ../../simulate_100_c/SRR10199716_output/error_sim_metagenome_R2.fastq ./SRR10199716_sim_2.fastq
ln -s ../../simulate_100_c/SRR10199724_output/error_sim_metagenome_R1.fastq ./SRR10199724_sim_1.fastq
ln -s ../../simulate_100_c/SRR10199724_output/error_sim_metagenome_R2.fastq ./SRR10199724_sim_2.fastq
ln -s ../../simulate_100_c/SRR10199725_output/error_sim_metagenome_R1.fastq ./SRR10199725_sim_1.fastq
ln -s ../../simulate_100_c/SRR10199725_output/error_sim_metagenome_R2.fastq ./SRR10199725_sim_2.fastq

for fastq in ./*fastq ; do
  fastqc $fastq -o ./
done
