*note: this file corresponds to the statonlabprivate wiki titled 'read_simulation'*

## directory structure
- README.md **(you're reading me)**
- analyses/ **(holds all specific objectives for processing data)**
- code/ **(holds specific software used for analyses)**
- raw_data/ **(unprocessed data)**
  - 1520_genomes/
- lab_notebook.md **(daily log of activities)**

## data

### 1520_genomes

https://doi.org/10.1038/s41587-018-0008-8

This dataset is relevant for simulating reads of similar sequencing preparation (human gut microbiota).

#### 2021.06.14
#### automate download of all accessions

Installed the 1,520 genome accessions produced in https://doi.org/10.1038/s41587-018-0008-8; wget called on the list of accession numbers. 


#### 2021.06.16
#### summarize 1520 accessions

Automatically download all accessions and provide summary statistics on the assembly.

Notes: The necessity to produce summary statistics for these assemblies for later analyses gave way to creating custom python scripts to pull the data. The custom scripts are stored in the /raw_data/genomes/1520_genomes directory with documentation.

