conda list > environment.log

mkdir cut_sites

for genome in ../../../../raw_data/snipen_RMS/mock_community_ref_genomes/*fna ; do
  python3 ../../../../code/scripts/cut_site_summaries.py -o ./cut_sites -m1 ecori -m2 msei -genome $genome
done

for query in ../3_samtools/SRR* ; do
  query_dir=$(basename $query)
  echo $query_dir
  mkdir $query_dir
  cd $query_dir
  for genome in ../../../../raw_data/snipen_RMS/mock_community_ref_genomes/*fna ; do
    genome_name=$(basename $genome)
    genome_name=${genome_name%%.fna}
    python3 ../../../../code/scripts/recreate_fragments.py ../../3_samtools/${query_dir}/${genome_name}/summary_${genome_name}.tsv ../cut_sites/cut_sites_${genome_name}.fna.csv
  done
  cd ..
done

