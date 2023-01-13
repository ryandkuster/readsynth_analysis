# refs is the list of all taxa hits from the naive search using kraken2 ('post_bracken_references.txt')
refs=$1

while read p; do
  #var1=$(< ../simulate_no_human/SRR5298272_post_bracken_10M_output/counts_${p}.gz.csv wc -l)
  var1=$(grep -c ",[CA][AC][TG][GT],0,[01]" ../simulate_no_human/SRR5298272_post_bracken_10M_output/raw_digest_${p}.fna.gz.csv)
  var2=$(grep -c ${p%%.fna} ../extract_features_combined/sim_10M_post_bracken/3_fragment_recreation/liu_sim_10M_post_bracken.csv)
  #echo $p
  #echo $var1
  #echo $var2
  echo ${p},${var1},${var2} >> inflated_taxa.csv
done < $refs
