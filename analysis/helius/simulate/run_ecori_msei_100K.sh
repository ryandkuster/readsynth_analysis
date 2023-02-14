# using readsynth commit 88d8bb13dd4e349c326745e43c3cd9ed3085eb3b

python3 /home/rkuster/readsynth/readsynth.py \
  -g abundances.csv \
  -o ecori_msei_100K_output/ \
  -n 100_000 \
  -u 400 \
  -sd 100 \
  -l 150 \
  -m1 ecori \
  -m2 msei