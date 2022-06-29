echo "#CHROM  POS combined" > combined_depths.tsv

for i in GCA* ; do
  sed -n "2,$"p ${i}/depth_${i}.txt >> combined_depths.tsv
done
