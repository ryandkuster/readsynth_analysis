rm insert_mean.csv
rm insert_sd.csv
rm insert_counts.csv
touch insert_mean.csv
touch insert_sd.csv
touch insert_counts.csv

for i in GCA* ; do
  sed -n "12"p ${i}/stats_fixmate_${i}.txt >> insert_counts.csv
  sed -n "40"p ${i}/stats_fixmate_${i}.txt >> insert_mean.csv
  sed -n "41"p ${i}/stats_fixmate_${i}.txt >> insert_sd.csv
done
