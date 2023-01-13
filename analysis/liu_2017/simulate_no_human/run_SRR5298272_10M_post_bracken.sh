# using readsynth commit 88d8bb13dd4e349c326745e43c3cd9ed3085eb3b

#python3 /home/rkuster/readsynth/scripts/sample_fastq.py -o ./ -r1 ../../../raw_data/liu_RMS/SRR5298272_R1.fastq -r2 ../../../raw_data/liu_RMS/SRR5298272_R2.fastq -p 10

python3 /home/rkuster/readsynth/readsynth.py \
  -g abundances_post_bracken.csv \
  -o SRR5298272_post_bracken_10M_output \
  -n 10_000_000 \
  -u 406 \
  -sd 100 \
  -l 75 \
  -m1 nlaiii \
  -m2 hpych4iv \
  -a1 for_nlaiii_liu_et_al_2017.tsv \
  -a2 rev_hpych4IV_liu_et_al_2017.tsv \
  -test
