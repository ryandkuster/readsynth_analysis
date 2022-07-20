import pandas as pd
import sys

'''
END GOAL: get the depth resulting from fragments with no internal cut sites
          adjust depth based on how frequent those fragments are constrained

in order:
 *  assign reads to taxonomic ids (using kraken2/bracken, primarily)
    pull out the reads corresponding to each taxid
    download the reference genome for each taxid
    for mock/simulated data, training data will use abundance as label
    map those reads to their reference genome (bwa, then samtools cleanup)
 *  get the fragment start/end postions at expected motif sites
 *      get the fragment lengths
 *      get the fragment internal cut sites
 *      get the fragment sequences (for kmer breakdown)
 *  samtools depth (no overlapping R1/R2 included) to collapse fragments
    consider feature matrix including a "percent at this length" feature
    adjust read depth using machine learning/linear regression
    assign an adjusted abundance to each fragment

other things to consider:
    how do these estimated abundances compare with raw reads?
    how do these estimated abundances compare with just using median depth?
'''
