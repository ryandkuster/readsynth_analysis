### data collected from:
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA482748

#### automated downloading of genome assembly statistics
PRJNA482748 accession sheet provided as sys.argv[1] in download_assembly_stats.py.

This script downloads (wget) all GenBank assembly statistics for each of the 1,520 genomes listed into the *assembly_stats/* directory and produces a summary csv of all the accessions.

#### automated downloading of genome assembly statistics
PRJNA482748 accession sheet provided as sys.argv[1] in download_genomes.py.

This script downloads (wget) all assembled genomes as 'fna.gz' into the *genomes/* directory and creates a collection of md5checksums for each.
