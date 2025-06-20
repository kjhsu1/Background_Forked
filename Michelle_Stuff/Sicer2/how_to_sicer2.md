How to use Sicer2

You do need your python in your yml file to run 3.10 as Sicer2 does not support python 3.13

1. Convert BAM file to BED file
make sure both samtools and bedtools are in yml dependencies
samtools view -h <name.bam> | bedtools bamtobed -i stdin <name.bed>

2. Install Sicer2
conda install -y numpy scipy
pip install SICER2
