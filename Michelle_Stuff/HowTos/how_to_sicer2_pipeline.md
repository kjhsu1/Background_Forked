How to Use the Sicer2 Pipeline: sicer2_pipeline.py

There are 4 required arguments.

In this order,

path to genome (in .fa or .fa.gz) <br>
path to directory containing all experimental reads<br>
path to directory containing all control reads<br>
name of genome<br>

This pipline creates 3 files.<br>
  index = needed for bowtie2<br>
  reults = this contains your results from sicer2<br>
  sam_bam_bed = this contains all the sam, bam and bed files created from samtools and sicer2
