# stdin argument 1 is the sample key for the 30 files

while read genome copies; do
    echo $genome $copies
  python3 /home/rkuster/readsynth/readsynth.py \
    -genome ../genomes/${genome} \
     -m1 G/AATTC \
     -m2 T/TAA \
     -complete 1 \
     -o ./output/ \
     -n $copies \
     -mean 300 \
     -sd 50 \
     -a1 ddRADseq_adapters_R1.txt \
     -a2 ddRADseq_adapters_R2.txt \
     -a1s 33 \
     -a2s 34 ;
done < $1
