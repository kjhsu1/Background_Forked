# make a sampler

"""
TO-DO FOR NEXT TIME 
1. CHECK IF WE NEED TO ADD ANYTHING TO THE PROGRAM
2. CHECK IF NO MISTAKES WERE MADE
3. TIDY UP THE PROGRAM
"""

import sys
import random
import matplotlib.pyplot as plt
import LIB
import json

"""
Assumptions for our Model
_________________________

- We will just make the assumption that all fragments are uniform in this model
- CHIP-seq reads are anywhere from 100-500bp approx.
"""


"""
Global Variables
________________
"""

# this particular genome has 3 chroms, chr1, chr2, chr3, with lengths 100, 200, 300 respectively
fasta = '../Genomes/random_genome_1.fa.gz'
coverage = 1 
num_bg_peaks = 1
num_fg_peaks = 1
k = 20  # fragment length


'''
Notes for this genome:
______________________

- total bp is 600
- with coverage =1, and k=20, we should see 30 samples getting pulled (CHECKED)
- also the chrom frequencies should be 1/6, 2/6, 3/6
    - meaning when we look at percentage of reads generated with respect to which chromosomes...
    - we should observe that it converges to above probabilities (CHECKED)

- to check if peaks are added properly we can just graph them
'''


"""
Functions
_________
"""


def add_peaks(pmf, num_peaks):
    """
    Take in a bins pmf, and adds num_peaks peaks to it
    """ 
    peaks = [2, 4, 6, 8, 9, 8, 6, 4, 2] # can change this however you want 
    used_peaks = [] # store used peaks
    p_index = random.randint(0, len(pmf) - len(peaks)+1) # init

    for peak in range(num_peaks):
        while p_index in used_peaks:
            p_index = random.randint(0, len(pmf) - len(peaks)) # can't create peaks on the last 6 bases
        for i in range(p_index, p_index + len(peaks)): # create the peaks in the bin
            print(i) # debug
            pmf[i] = pmf[i] + peaks[i-p_index] # should we add or multiply?
        # don't want peaks to overlap, and want a space of at least 1 between peaks
        start = p_index - 1
        end = p_index + len(peaks) + 1
        used_peaks.extend(range(start, end))
    return pmf 

def normalize_bins(pmf):
    """
    Takes in  a bins pmf, normalize so all bins sum to 1
    """
    total = sum(pmf)
    for i, bin in enumerate(pmf):
        pmf[i] = (pmf[i]/total)
    return pmf

def sample_from_bins(pmf, num_sample):
    """
    Takes in a bins pmf, and samples num_sample times from the distribution.
    """
    bins_indices = list(range(len(pmf))) # ex. range(3) = [0, 1, 2]
    samples = random.choices(bins_indices, weights=pmf, k=num_sample)
    return samples

def create_pmf(chrom_len, num_bg_peaks, num_fg_peaks, k):
    """
    Create a pmf for all possible fragments combining the background and foreground pmfs. 
    """
    
    num_bins = chrom_len - k + 1

    # initialize bins 
    pmf = [1] * num_bins

    # add background peaks
    pmf = add_peaks(pmf, num_bg_peaks)

    # add foreground peaks
    pmf = add_peaks(pmf, num_fg_peaks)

    pmf = normalize_bins(pmf)

    return pmf

def create_pmf_all_chroms(fasta):
    """
    Create pmf for each chromosome in the genome
    
    Input: Path to FASTA

    Output: Dictionary; 
        key: chromosome (ex. chr1, chr2) 
        value:  pmf for the chromosome as a list (ex. [0.001, 0.002, 0.001, ...])
    """
    genome_pmfs = {}
    for id, seq in LIB.read_fasta(fasta):
        genome_pmfs[id] = create_pmf(len(seq), num_bg_peaks, num_fg_peaks, k)
    # print(json.dumps(genome_pmfs, indent=4))
    return genome_pmfs

def chrom_bias(fasta):
    """
    Find the sampling bias for each chrom based on length relative to the sum total of genomic bps
    """
    chrom_bias = {}
    total_bp = 0
    for id, seq in LIB.read_fasta(fasta):
        total_bp += len(seq)
        chrom_bias[id] = len(seq)
    for id in chrom_bias:
        chrom_bias[id] = chrom_bias[id] / total_bp
    return chrom_bias

def sample_genome(fasta, genome_pmfs):
    """
    Sample from a genome pmf
    
    Parameters: Path to fasta, genomic P.M.F

    Output: dictionary of reads
        key: chrom (ex. 'chr1', 'chr2')
        value: list of lists of reads, (ex. [ [1,2,3,4], [4,5,6,7] ])
            - in the example, it denotes two reads; both 4 bps long spanning those coords 
    """
    # find chrom bias
    chrom_bias = {}
    total_bp = 0
    for id, seq in LIB.read_fasta(fasta):
        total_bp += len(seq)
        chrom_bias[id] = len(seq)
    for id in chrom_bias:
        chrom_bias[id] = chrom_bias[id] / total_bp
    
    # Do the experiment
    num_sample = int((total_bp * coverage) / k)
    chroms = list(chrom_bias.keys())
    biases = list(chrom_bias.values())
    reads_dict = {} # store samples for each chrom
    for chrom in chroms:
        reads_dict[chrom] = []

    for i in range(num_sample):
        picked_chrom = random.choices(chroms, weights=biases)[0]
        sample_index = sample_from_bins(genome_pmfs[picked_chrom], 1)[0]
        sample = list(range(sample_index, sample_index+k))
        reads_dict[picked_chrom].append(sample)

    # print(json.dumps(reads_dict, indent=4))
    return reads_dict

genome_pmf = create_pmf_all_chroms(fasta) # create genome pmf based on genome fasta
exp = sample_genome(fasta, genome_pmf) # generate sample from it

#print(json.dumps(exp, indent=4))

"""
Functions for Debugging
_______________________
"""

def total_num_reads(exp_sample):
    """
    Get total number of reads from a sampling experiment
    """
    read_sum = 0
    for chrom in exp_sample:
        read_sum += len(exp_sample[chrom])
    return read_sum

def chrom_distribution(exp_sample):
    """
    Get distribution of chromosomes in a sampling experiment

    WILL NOT RETURN VALUE, WILL JUST PRINT
    """
    total = total_num_reads(exp_sample)
    dict = {}
    for chrom in exp_sample:
        chrom_total = len(exp_sample[chrom])
        dict[chrom] = chrom_total / total
    print('Experimental Distribution', '\n',json.dumps(dict, indent=4))
    
    # can compare with expected distribution
    c_bias = chrom_bias(fasta)
    print('Expected Distribution', '\n', json.dumps(c_bias, indent=4))

def graph_pmf(pmf, title="Experimental P.M.F"):
    """
    Graphs the P.M.F
    """
    bin_coords = range(len(pmf))
    probs = pmf

    plt.figure()
    plt.plot(bin_coords, probs)
    plt.xlabel('Read/K-mer index')
    plt.ylabel('Probability')
    plt.title(title)
    plt.show()

def graph_all_genome_pmf(genome_pmf):
    for chrom in genome_pmf.keys():
        pmf = genome_pmf[chrom]
        graph_pmf(pmf, title=f'P.M.F for Experimental, {chrom}')

graph_all_genome_pmf(genome_pmf)

'''
bins_dict = {}
for b in bins_indices:
    bins_dict[b]: 0 # add new dict key-value pair for each bin number

for sample in samples:
    if sample in bins_dict:
        bins_dict[sample]
'''











