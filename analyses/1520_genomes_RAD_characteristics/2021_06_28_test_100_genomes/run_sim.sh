for i in ../genomes/*.fna; do
  python3 /home/rkuster/readsynth/readsynth.py \
    -genome $i \
    -m1 T/TAA G/AATTC \
    -m2 T/TAA G/AATTC \
    -t 1 \
    -complete 1 \
    -o ./output \
    -n 1000 \
    -mean 400 \
    -sd 100 ;
done

