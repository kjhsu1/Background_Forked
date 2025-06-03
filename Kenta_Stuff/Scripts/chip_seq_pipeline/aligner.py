"""
This program takes in FASTA of CHIP-seq reads and aligns them to a reference genome using Bowtie2
Output: SAM file

NOTE: 
    - You need to have created a index file for your genome before running this program
    - You only need to create one for a genome once (and can use that for any proceeding bowtie2 alignment)
    - You can create one using command:
        bowtie2-build [path/to/gzipped/genome/fasta/file] [path/to/new/index/file]
        - ex. bowtie2-build random_genome_1.fa.gz index/random_genome_1
    - New index file name should be genome file name without suffix (ie. random_genome_1)
"""

"""
Imports
_______
"""

import subprocess
import argparse

"""
Arguments
________
"""

parser = argparse.ArgumentParser(description="bowtie2 arguments")
parser.add_argument('--genome_index_path', type=str, required=True,
    help='Path to the Genome Index File')
parser.add_argument('--fasta_path', type=str, required=True, 
    help='Path to the FASTA for CHIP-seq reads')
parser.add_argument('--output_path', type=str, required=True,
    help='Path for the output; Make Sure the Path Base is the File Name You Want For the Output SAM')

arg = parser.parse_args()

genome_index_path = arg.genome_index_path
fasta_path = arg.fasta_path
output = arg.output_path


"""
Functions
_________
"""

def run_bowtie2(genome_index_path, reads_fasta_path, output_path):
    """
    Runs bowtie2 and aligns a FASTA of a CHIP-seq experiment reads 
    """
    cmd = [
        "bowtie2",                      # the tool
        "-x", genome_index_path,        # path to index (no file extension)
        "-f",                           # input is FASTA (not FASTQ)
        "-U", reads_fasta_path,         # input reads
        "-S", output_path               # output SAM file
    ]

    try:
        subprocess.run(cmd, check=True)
        print("Alignment Success")
    except subprocess.CalledProcessError as e:
        print(f"Alignment Failed using Bowtie2: {e}")

"""
Main
_____
"""

# run
run_bowtie2(genome_index_path, fasta_path, output)