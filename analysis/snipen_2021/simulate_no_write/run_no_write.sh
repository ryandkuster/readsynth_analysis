# using readsynth commit 88d8bb13dd4e349c326745e43c3cd9ed3085eb3b

/home/rkuster/readsynth/readsynth.py \
  -g abundances.csv \
  -o no_write_output/ \
  -lp 1200 \
  -n 658_194 \
  -l1 145 \
  -l2 151 \
  -m1 ecori \
  -m2 msei \
  -c 0.95 \
  -a1 for_ecori_ravi_et_al_2018.tsv \
  -a2 rev_msei_ravi_et_al_2018.tsv \
  -a1s 64 \
  -a2s 64 \
  -test
