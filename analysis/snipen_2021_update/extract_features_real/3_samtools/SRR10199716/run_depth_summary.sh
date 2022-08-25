echo "#CHROM  POS combined" > motif_combined_depths.tsv

for i in GCA* ; do
  sed -n "2,$"p ${i}/motif_depth_${i}.txt >> motif_combined_depths.tsv
done
