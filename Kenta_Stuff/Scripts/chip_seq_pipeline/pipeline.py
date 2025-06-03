"""
Can Take Two CHIP-seq reads as FASTA (one for experiment, one for control) and Call Peaks With MACS2
"""

"""
Imports
______
"""

import os
import argparse
import subprocess

"""
User Input Arguments
___________
"""

parser = argparse.ArgumentParser( 
    description='Can Take Two CHIP-seq reads as FASTA (one for experiment, one for control) and Call Peaks With MACS2'
)
parser.add_argument('--exp_fasta', type=str, required=True,
                    help='Path to Experiment CHIP-seq Read FASTA; Prefix of This Filename Will Be Used in SAM and BAMs')
parser.add_argument('--control_fasta', type=str, required=True,
                    help='Path to Control/Background CHIP-seq Read FASTA; Prefix of This Filename Will Be Used in SAM and BAMs')
parser.add_argument('--genome_index', type=str, required=True,
                    help='Path to Genome Index File For Bowtie2 Alignment')
parser.add_argument('--control_output_path', type=str, required=True,
                    help='Path to Output Path For Control/Background; Just the Directory Name is Fine; \
                    SAM, All BAM Variations Will Be Stored Here')
parser.add_argument('--exp_output_path', type=str, required=True,
                    help='Path to Output Path For Treatment/Experiment; Just the Directory Name is Fine; \
                    SAM, All BAM Variations Will Be Stored Here')
parser.add_argument('--peaks_output_path', type=str, required=True,
                    help='Path to Peak Caller Output Directory')
parser.add_argument('--genome_size', type=str, required=True,
                    help='Sum Up # of Bases of All Chromosomes in Genome')
arg = parser.parse_args()

# FASTA Paths 
exp_fasta = arg.exp_fasta
control_fasta = arg.control_fasta
# Genome Index File Path
genome_index = arg.genome_index
# Directory to Store SAM, BAM
exp_output_path = arg.exp_output_path
control_output_path = arg.control_output_path
# Peak Caller Output Directory
peaks_output_path = arg.peaks_output_path
# Genome Size
genome_size = arg.genome_size

"""
Functions
________
"""

def run_bowtie2(genome_index_path, reads_fasta_path, output_path):
    """
    Runs bowtie2 and aligns a FASTA of a CHIP-seq experiment reads

    Parameters:
        - genome_index_path: Full Path to Genome Index File
        - reads_fasta_path: Full Path to Chip-seq Reads FASTA
        - output_path: Path to Directory That Will Store SAM, BAMs (no file name)
    """
    basename = os.path.basename(reads_fasta_path)
    prefix, _ = os.path.splitext(basename)
    newname = prefix + '.sam'
    sam_path = os.path.join(output_path, newname) # derived new path for output SAM file

    cmd = [
        "bowtie2",                      # the tool
        "-x", genome_index_path,        # path to index (no file extension)
        "-f",                           # input is FASTA (not FASTQ)
        "-U", reads_fasta_path,         # input reads
        "-S", sam_path               # output SAM file
    ]

    try:
        subprocess.run(cmd, check=True)
        print("Alignment Success")
    except subprocess.CalledProcessError as e:
        print(f"Alignment Failed using Bowtie2: {e}")

def run_sam_to_bams(reads_fasta_path, output_path):
    """
    Convert SAM to BAMs Using Program sam_to_bams.py
    Input:
        - reads_fasta_path: Full Path to Chip-seq Reads FASTA
        - output_path: Path to Directory That Will Store SAMs and BAMs
    """

    basename = os.path.basename(reads_fasta_path)
    prefix, _ = os.path.splitext(basename)
    newname = prefix + '.sam'
    sam_path = os.path.join(output_path, newname) # derived path for SAM file

    # derive BAM path
    newname = prefix + '.bam'
    bam_path = os.path.join(output_path, newname)

    cmd = [
        'python3', 'sam_to_bams.py',
        '--sam_path', sam_path,
        '--bam_output', bam_path
    ]

    subprocess.run(cmd, check=True)

def run_bam_to_peaks(exp_output_path, control_output_path, exp_reads_fasta_path, control_reads_fasta_path, genome_size, peaks_output_path):
    """
    Calls Peaks With MACS2 Through Program, bam_to_peaks_called.py
    """
    # path string manipulation for exp
    basename = os.path.basename(exp_reads_fasta_path)
    prefix, _ = os.path.splitext(basename)
    exp_newname = 'sorted.' + prefix + '.bam'

    # same for control
    basename = os.path.basename(control_reads_fasta_path)
    prefix, _ = os.path.splitext(basename)
    control_newname = 'sorted.' + prefix + '.bam'

    # derive the paths for the sorted BAMs
    exp_sorted_bam_path = os.path.join(exp_output_path, exp_newname)
    control_sorted_bam_path = os.path.join(control_output_path, control_newname)

    cmd = [
        'python3', 'bam_to_peaks_called.py',
        '--exp_sorted_bam_path', exp_sorted_bam_path,
        '--control_sorted_bam_path', control_sorted_bam_path,
        '--genome_size', genome_size,
        '--output_dir', peaks_output_path
    ]

    subprocess.run(cmd, check=True)

"""
Main
________
"""

def main():
    # align treatment first
    run_bowtie2(genome_index, exp_fasta, exp_output_path)
    # then align background
    run_bowtie2(genome_index, control_fasta, control_output_path)

    # Convert Treatment SAM to BAMs
    run_sam_to_bams(exp_fasta, exp_output_path)
    # Next Control
    run_sam_to_bams(control_fasta, control_output_path)

    # Call Peaks
    run_bam_to_peaks(exp_output_path, control_output_path, exp_fasta, control_fasta, genome_size, peaks_output_path)

main()