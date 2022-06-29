# using readsynth commit 8336a061d82cc3ee460adabd72bb6b36e7e63a8a

#python3 /home/rkuster/readsynth/scripts/sample_fastq.py -o ./ -r1 ../../../raw_data/sun_2bRADM/2brad_msa/mock-1ng.BcgI.fq -r2 ../../../raw_data/sun_2bRADM/2brad_msa/mock-1ng.BcgI.fq -p 10

/home/rkuster/readsynth/readsynth.py \
  -g abundances.csv \
  -o output_mock-1ng.BcgI/ \
  -u 32 \
  -sd 10 \
  -n 2_407_174 \
  -l1 32 \
  -l2 32 \
  -iso bcgi \
  -a1 for_2brad.tsv \
  -a2 rev_2brad.tsv \
  -a1s 4 \
  -a2s 4 \
  -q1 mock-1ng.BcgI.fq_sampled_scores.csv \
  -q2 mock-1ng.BcgI.fq_sampled_scores.csv \
