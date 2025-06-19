import os
import sys
import argparse

parser = argparse.ArgumentParser(description='runs multiple coverages at once')
parser.add_argument('genome', help='genome')
parser.add_argument('exp', help='dir of diff exp coverages')
parser.add_argument('ctrl', help='dir of control fasta files')
parser.add_argument('--output', default='.', help='output destination')
arg = parser.parse_args()

exp_path = '../../Shared/Reads_FASTA/Experimental_unbiased/'
ctrl_path = '../../Shared/Reads_FASTA/unbiased_control_reads/'
wd = arg.output

os.system(f'rm -rf {wd}/index')
os.system(f'mkdir {wd}/index')
os.system(f'bowtie2-build rg_1.fa.gz {wd}/index/rg_1')
# print('initialization done')
print()
# sys.exit()

'''
## temp ##
ctrl = ''.join([ctrl_path, 'unbiased_control_cov_1_random_genome_1_reads.fa'])
ctrl_out = 'unbiased_control_cov_1_random_genome_1_reads'
## temp ##
'''
counter = 0
i = 0

for ctrl_file in os.listdir(arg.ctrl):
    ctrl = ''.join([ctrl_path, ctrl_file])
    ctrl_out = ctrl_file.strip('.fa')

    os.system(f'bowtie2 -x {wd}/index/rg_1 -f -U {ctrl} -S {wd}/{ctrl_out}_output.sam')
    os.system(f'samtools view -bS {wd}/{ctrl_out}_output.sam | samtools sort -o {wd}/{ctrl_out}_sorted.bam')
    os.system(f'makeTagDirectory {wd}/{ctrl_out}_tag_dir {wd}/{ctrl_out}_sorted.bam')

    for exp_file in os.listdir(arg.exp):
        exp = ''.join([exp_path, exp_file])
        exp_out = exp_file.strip('.fa')
        
        # ensures these files are only made once
        if i == 0:
            os.system(f'bowtie2 -x {wd}/index/rg_1 -f -U {exp} -S {wd}/{exp_out}_output.sam')
            os.system(f'samtools view -bS {wd}/{exp_out}_output.sam | samtools sort -o {wd}/{exp_out}_sorted.bam')
            os.system(f'makeTagDirectory {wd}/{exp_out}_tag_dir {wd}/{exp_out}_sorted.bam')

        os.system(f'findPeaks {wd}/{exp_out}_tag_dir -i {wd}/{ctrl_out}_tag_dir > {wd}/homer_{exp_out}_{ctrl_out}.txt')
    i += 1
    # if counter==1: break
    # else: counter+=1

os.system(f'rm -rf {wd}/results')
os.system(f'mkdir {wd}/results')
os.system(f'mv {wd}/homer* {wd}/results')

'''
Notes


###############
Bowtie2 section
###############

os.system('mkdir index')
os.system('bowtie2-build random_genome_1.fa.gz index/rg_1')
os.system('bowtie2 -x index/rg_1 -f -U exp_1.fa -S exp_1_output.sam')
os.system('bowtie2 -x index/rg_1 -f -U ctrl_1.fa -S ctrl_1_output.sam')

#############
Homer section
#############

os.system('samtools view -bS exp_1_output.sam | samtools sort -o exp_1_sorted.bam')
os.system('samtools view -bS ctrl_1_output.sam | samtools sort -o ctrl_1_sorted.bam')
os.system('makeTagDirectory exp_1_tag_dir exp_1_sorted.bam')
os.system('makeTagDirectory ctrl_1_tag_dir ctrl_1_sorted.bam')
os.system('findPeaks exp_1_tag_dir -i ctrl_1_tag_dir > homer_exp1_ctrl1.txt')
os.system('ls -F')

'''