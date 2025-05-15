# make a sampler

import sys
import random
import matplotlib.pyplot as plt
import LIB
import json

"""
- We will just make the assumption that all fragments are uniform in this model
- CHIP-seq reads are anywhere from 100-500bp approx.
"""

'''
genome_size = sys.argv[1] # 6e10 
coverage = sys.argv[2] # 10 = 10X
'''

"""
Input Variables
"""
fasta = '../Genomes/random_genome_1.fa.gz'
coverage = 1
num_bg_peaks = 1
num_fg_peaks = 1
k = 20  # fragment length

"""
Functions
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
            p_index = random.randint(0, len(pmf) - len(peaks)+1) # can't create peaks on the last 6 bases
        for i in range(p_index, p_index + len(peaks)): # create the peaks in the bin
            pmf[i] = pmf[i] + peaks[i-p_index] # should we add or multiply?
        # don't want peaks to overlap, and want a space of at least 1 between peaks
        start = p_index - len(peaks) + 1
        end = p_index + len(peaks) + 1
        used_peaks.extend(range(start, end+1))
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

# we need to be able to create pmfs for each chr
genome_pmfs = {}

# create pmfs for all chr
for id, seq in LIB.read_fasta(fasta):
    genome_pmfs[id] = create_pmf(len(seq), num_bg_peaks, num_fg_peaks, k)
# print(json.dumps(genome_pmfs, indent=4))

# find the sampling bias for each chrom based on relative length
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
















'''
# let's graph the pmf
bin_coords = range(genome_size - k + 1)
pmf_bins = exp_pmf

plt.figure()
plt.plot(bin_coords, pmf_bins)
plt.xlabel('Bin index')
plt.ylabel('Probability')
plt.title('Experimental P.M.F')
plt.show()
'''

'''
bins_dict = {}
for b in bins_indices:
    bins_dict[b]: 0 # add new dict key-value pair for each bin number

for sample in samples:
    if sample in bins_dict:
        bins_dict[sample]
'''











