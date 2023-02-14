# using readsynth commit 8336a061d82cc3ee460adabd72bb6b36e7e63a8a

#python3 /home/rkuster/readsynth/scripts/sample_fastq.py -o ./ -r1 ../../../raw_data/snipen_RMS/SRR10199716_1.fastq -r2 ../../../raw_data/snipen_RMS/SRR10199716_2.fastq -p 10

python3 /home/rkuster/readsynth/readsynth.py \
  -g abundances_no_human.csv \
  -o SRR5298272_no_human_output/ \
  -n 50_000_000 \
  -u 300 \
  -sd 50 \
  -c 0.80 \
  -l1 75 \
  -l2 75 \
  -m1 nlaiii \
  -m2 hpych4IV \
  -a1 for_nlaiii_liu_et_al_2017.tsv \
  -a2 rev_hpych4IV_liu_et_al_2017.tsv \
  -a1s 73 \
  -a2s 71
