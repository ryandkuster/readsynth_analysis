# using readsynth commit 8336a061d82cc3ee460adabd72bb6b36e7e63a8a

python3 /home/rkuster/readsynth/scripts/sample_fastq.py -o ./ -r1 ../../../raw_data/snipen_RMS/SRR10199725_1.fastq -r2 ../../../raw_data/snipen_RMS/SRR10199725_2.fastq -p 10

/home/rkuster/readsynth/readsynth.py \
  -g abundances.csv \
  -o SRR10199725_output/ \
  -u 241 \
  -sd 91 \
  -n 789_892 \
  -l1 145 \
  -l2 151 \
  -m1 ecori \
  -m2 msei \
  -c 0.9759 \
  -a1 for_ecori_ravi_et_al_2018.tsv \
  -a2 rev_msei_ravi_et_al_2018.tsv \
  -a1s 64 \
  -a2s 64 \
  -e A \
  -q1 SRR10199725_1.fastq_sampled_scores.csv \
  -q2 SRR10199725_2.fastq_sampled_scores.csv
