# How I got bowtie2 to work

## Conda environment creation
Start by creating a conda environment using `conda env create -f control_environment.yml`.

Activate using `conda activate control_exp`.

## Set up
Make sure `exp_1.fa` and `random_genome_1.fa.gz` is in pwd.

Use `mkdir index` to create directory to store bowtie2 index.

Run `bowtie2-build random_genome_1.fa.gz index/random_genome_1` to create an index inside of `index` directory.

## Running bowtie2
Run command `bowtie2 -x index/random_genome_1 -f -U exp_1.fa -S exp_1_output.sam`.

`-x` flag indicates which index to use.

`-f` tells bowtie2 that the query is in fasta format and `-U` indicates to the query.

`-S` indicates name for output file (in SAM format).