# using readsynth commit 88d8bb13dd4e349c326745e43c3cd9ed3085eb3b

python3 /home/rkuster/readsynth/readsynth.py \
  -g abundances.csv \
  -o hhai_msei_output/ \
  -n 1_000_000 \
  -u 400 \
  -sd 100 \
  -l 150 \
  -m1 hhai \
  -m2 msei
