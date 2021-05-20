A snakemake workflow for calling variants using GATK on NGS data

sequences_metadata.tsv needs to be a tsv file and have within it the following columns named: unique_id, Sample, fq1, fq2, rg

- unique_id: character string from set [\w]
- Sample: character string from set [\w]
- fq1/fq2: file path of the fastq file
- rg: character string of the read group information eg @RG\tID:flowcell_lane\tSM:sample\tPL:illumina\tLB:library\tPU:flowcell_lane

A good method for creating unique_id is to use "{flowcell}_{lane}"

see [this article from the GATK about read groups](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) for more information

## Setting up the pipeline

```.bash
# On the Biochem servers:

# Load conda
$ module load miniconda/Miniconda3_4.8.3
# Create the environment
$ conda env create -f environment.yaml -p env

# create links to the needed reference files (GATK resource bundles)

```

## Running the pipeline:

Make sure that conda is installed and available and activate the environment


```.bash
# On biochem servers:

# Load conda
$ module load miniconda/Miniconda3_4.8.3

# activate the environment
$ conda activate ./env
```


Once the environment is activated make sure the sequences_metadata.tsv file is present and you can try a dry run using:

```.bash
$ Snakemake -nr
```

If that succeeds then the pipeline can be started using
```.bash
$ Snakemake
```
