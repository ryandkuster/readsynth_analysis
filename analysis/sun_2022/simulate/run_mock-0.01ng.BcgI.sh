# using readsynth commit 8336a061d82cc3ee460adabd72bb6b36e7e63a8a

python3 /home/rkuster/readsynth/scripts/sample_fastq.py -o ./ -r1 ../../../raw_data/REPLACE -r2 ../../../raw_data/REPLACE -p 10

/home/rkuster/readsynth/readsynth.py \
  -g abundances.csv \
  -o REPLACE_output/ \
  -d REPLACE_dist.json \
  -n 56_197 \
  -l1 32 \
  -l2 32 \
  -iso bcgi \
  -a1 REPLACE.tsv \
  -a2 REPLACE.tsv \
  -a1s 2 \
  -a2s 2 \
  -q1 REPLACE.fastq_sampled_scores.csv \
  -q2 REPLACE.fastq_sampled_scores.csv
