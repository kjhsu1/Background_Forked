import os
import sys # used for debugging
import argparse

parser = argparse.ArgumentParser(description='runs multiple coverages at once')
parser.add_argument('genome', help='genome')
parser.add_argument('exp', help='dir of diff exp coverages')
parser.add_argument('ctrl', help='dir of control fasta files')
parser.add_argument('genome_size', help='what type and size the genome is')
parser.add_argument('--output', default='.', help='output destination')
arg = parser.parse_args()

# working directory (aka output directory)
wd = arg.output

################################
## initialization for bowtie2 ##
################################

os.system(f'rm -rf {wd}/index')
os.system(f'mkdir {wd}/index')
os.system(f'bowtie2-build {arg.genome} {wd}/index/genome > {wd}/index/log.txt')
os.system(f'mkdir {wd}/sam_bam_bed')

#################################
### running bowtie2 and sicer2 ###
#################################

i = 0 # used to track first run

for ctrl_file in os.listdir(arg.ctrl):
    ctrl_out = ctrl_file.strip('.fa')

    os.system(f'bowtie2 -x {wd}/index/genome -f -U {arg.ctrl}/{ctrl_file} -S {wd}/{ctrl_out}_output.sam')
    os.system(f'samtools view -bS {wd}/{ctrl_out}_output.sam | samtools sort -o {wd}/{ctrl_out}_sorted.bam')
    os.system(f'mv {wd}/{ctrl_out}_sorted.bam {wd}/sam_bam_bed')
    os.system(f'mv {wd}/{ctrl_out}_output.sam {wd}/sam_bam_bed')

    for exp_file in os.listdir(arg.exp):
        exp_out = exp_file.strip('.fa')
        
        # ensures these files are only made once
        if i == 0:
            os.system(f'bowtie2 -x {wd}/index/genome -f -U {arg.exp}/{exp_file} -S {wd}/{exp_out}_output.sam')
            os.system(f'samtools view -bS {wd}/{exp_out}_output.sam | samtools sort -o {wd}/{exp_out}_sorted.bam')
            os.system(f'mv {wd}/{exp_out}_sorted.bam {wd}/sam_bam_bed')
            os.system(f'mv {wd}/{exp_out}_output.sam {wd}/sam_bam_bed')
            
        os.system(f'sicer -t {wd}/sam_bam_bed/{exp_out}_sorted.bam -c {wd}/sam_bam_bed/{ctrl_out}_sorted.bam --species {arg.genome_size} --o sicer2_{exp_out}_{ctrl_out}')
    i += 1

os.system(f'rm -rf {wd}/results')
os.system(f'mkdir {wd}/results')
os.system(f'mv {wd}/sicer2_{exp_out}_{ctrl_out} {wd}/results')
os.system(f'mv {wd}/{exp_out}_sorted.bed {wd}/sam_bam_bed')