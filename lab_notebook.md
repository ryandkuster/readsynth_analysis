### Read Simulation

#### 2021.06.14
#### automate download of sratoolkit accessions

Objective: Install the 1,520 genome accessions produced in https://doi.org/10.1038/s41587-018-0008-8, sratoolkit will be called on a list of accession numbers. 

Test 1: A smaller sample set of accessions will be provided to confirm the script works.
Test 2: A selected number of the 1,520 accessions will be downloaded to estimate overall size and to confirm proper folder structure.

Results: sratoolkit unnecessary, but will need to use the Entrez utility or automate wget downloads for each genome.

#### 2021.06.16
#### automate download of 1520 accessions

Objective: automatically download all accessions and provide summary statistics on the assembly.

Notes: The necessity to produce summary statistics for these assemblies for later analyses gave way to creating custom python scripts to pull the data. The custom scripts are stored in the /raw_data/genomes/1520_genomes directory with documentation.

Results:
