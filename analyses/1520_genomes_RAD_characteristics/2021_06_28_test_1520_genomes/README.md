## 2021.06.28

### 1_RE_digest_genomes
#### in silico digest each of the 1,520 genomes using MseI and EcoRI motifs using 'readsynth' scripts.

readsynth.py from commit 1e72bf20d58ef8b7af03edc661cca810a3598b85

For the purposes of exploratory analysis, only RE digestion (complete) will be performed using complete digestion. The MseI and EcoRI have no issue of overlapping recognition sites, but this needs to be addressed for future simulations, as duplicate loci will be represented in these events.

The output of each simulation will consist of 4 files:
- raw_digest...csv (the absolute, raw fragments produced by the simulated RE digest)
- hist_raw_digest...pdf (histogram of the raw fragments)
- sampled...csv (the absolute, raw fragments produced by the simulated RE digest)
- hist_sampled (histogram of the sampled fragments)

*note* The raw_digest...csv files will contain all possible cut sites considering both template and reverse-complemented template orientations. The originating strand is represented in the 'strand' column using a '+' (template) or '-' (reverse). Therefore, if 10,000 fragments are produced, we can assume that there are truly only 5,0000 fragments, as long as digests are complete and there is no strict orientation of cut sites on the 5' or 3' end. If there is an orientation requirement (e.g. adapter ligation requires that MseI be only on the 5' end and EcoRI can only be on the 3' end) then it makes sense to store reads in this strand-specific manner.

By default the maximum fragment length is set to mean + 6sd, so for these test, using mean of 400 and standard deviation of 100, we can expect the maximum fragment length to be included to max out at 1000bp.


## 2021.06.29

### 2_pull_digest_stats
#### for every template strand ('+') from the raw_digest...csv file, pull the total number of fragments as well as ranges of fragment sizes (e.g. 1 to 200bp, 201 to 400bp, 401 to 600bp).
