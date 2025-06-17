import os
import sys
import argparse

parser = argparse.ArgumentParser(description='runs multiple coverages at once')
parser.add_argument('genome', type=str, help='genome')
parser.add_argument('coverages', type=str, help='dir of diff coverages')
arg = parser.parse_args()

path = '../../../Kenta_Stuff/Generated_Data/Reads_FASTA/Different_Coverages_Unbiased/'

for file in os.listdir(arg.coverages):
    file_path = ''.join([path,file])
    out = file.strip('.fa')

    os.system('rm -rf index')
    os.system('mkdir index')
    # print('mkdir done')
    os.system('bowtie2-build ../random_genome_1.fa.gz index/random_genome_1')
    # print('index done')
    os.system(f'bowtie2 -x index/random_genome_1 -f -U {file_path} -S {out}_output.sam')
    # print('bowtie2 done')


    os.system(f'samtools view -bS {out}_output.sam | samtools sort -o {out}_sorted.bam')
    os.system(f'makeTagDirectory {out}_tag_dir {out}_sorted.bam')
    os.system(f'findPeaks {out}_tag_dir > {out}_homer.txt')


sys.exit()
"""
Bowtie2 section
"""

os.system('mkdir index')
os.system('bowtie2-build random_genome_1.fa.gz index/random_genome_1')
os.system('bowtie2 -x index/random_genome_1 -f -U exp_1.fa -S exp_1_output.sam')

"""
Homer section
"""

os.system('samtools view -bS exp_1_output.sam | samtools sort -o exp_1_sorted.bam')
os.system('makeTagDirectory Exp_1_tag_dir exp_1_sorted.bam')
os.system('findPeaks Exp_1_tag_dir > exp_1_homer.txt')
os.system('ls -F')