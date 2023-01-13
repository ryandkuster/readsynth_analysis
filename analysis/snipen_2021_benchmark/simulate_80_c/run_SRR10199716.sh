# using readsynth commit 88d8bb13dd4e349c326745e43c3cd9ed3085eb3b

python3 /home/rkuster/readsynth/scripts/sample_fastq.py -o ./ -r1 ../../../raw_data/snipen_RMS/SRR10199716_1.fastq -r2 ../../../raw_data/snipen_RMS/SRR10199716_2.fastq -p 10

/home/rkuster/readsynth/readsynth.py \
  -g abundances.csv \
  -o SRR10199716_output/ \
  -u 143 \
  -sd 90 \
  -n 658_194 \
  -l1 145 \
  -l2 151 \
  -m1 ecori \
  -m2 msei \
  -c 0.80 \
  -a1 for_ecori_ravi_et_al_2018.tsv \
  -a2 rev_msei_ravi_et_al_2018.tsv \
  -a1s 64 \
  -a2s 64 \
  -e A \
  -q1 SRR10199716_1.fastq_sampled_scores.csv \
  -q2 SRR10199716_2.fastq_sampled_scores.csv
