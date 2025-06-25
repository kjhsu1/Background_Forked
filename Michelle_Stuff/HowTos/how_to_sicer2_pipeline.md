How to Use the Sicer2 Pipeline: sicer2_pipeline.py

There are 4 required arguments.

In this order,

path to genome (in .fa or .fa.gz)
path to directory containing all experimental reads
path to directory containing all control reads
name of genome

This pipline creates 3 files.
  index = needed for bowtie2
  reults = this contains your results from sicer2
  sam_bam_bed = this contains all the sam, bam and bed files created from samtools and sicer2
