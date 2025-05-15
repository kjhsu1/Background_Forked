# make a sampler

import sys
import random
import matplotlib.pyplot as plt

'''
genome_size = sys.argv[1] # 6e10 
coverage = sys.argv[2] # 10 = 10X
'''
genome_size = 320
coverage = 5
num_peaks = 2

# CHIP-seq reads are 100-500bp approx.
# let's just say all are 300bp
k = 30
# how many 300bp frag to sample
num_sample = int((genome_size * coverage) / k)
# sliding window 
# L - k + 1 
num_bins = genome_size - k + 1
# assign prob. for each bin 
bins = [1] * num_bins

# create random peaks 
peaks = [2, 4, 6, 8, 9, 8, 6, 4, 2] # can change this however you want 
used_peaks = [] # store used peaks
p_index = random.randint(0, len(bins) - len(peaks)+1) # init

for peak in range(num_peaks):
    while p_index in used_peaks:
        p_index = random.randint(0, len(bins) - len(peaks)+1) # can't create peaks on the last 6 bases
    for i in range(p_index, p_index + len(peaks)): # create the peaks in the bin
        bins[i] = peaks[i-p_index]
    # don't want peaks to overlap, and want a space of at least 1 between peaks
    start = p_index - len(peaks) + 1
    end = p_index + len(peaks) + 1
    used_peaks.extend(range(start, end+1))

# add true signal distribtion to the pmf 


# make sure all bins sum to 1 as probability
total = sum(bins)
for i, bin in enumerate(bins):
    bins[i] = (bins[i]/total)

# should have created the bins that add up to 1, each probabilities
# in other words we have the pmf of the population

# now we can sample 
bins_indices = list(range(len(bins))) # ex. range(3) = [0, 1, 2]
samples = random.choices(bins_indices, weights=bins, k=num_sample)

print(num_sample)
print(len(bins))
print(samples)

# print(sum(bins))

# let's graph the pmf
'''
bin_coords = range(len(bins))
pmf_bins = bins

plt.figure()
plt.plot(bin_coords, pmf_bins)
plt.xlabel('Bin index')
plt.ylabel('Probability')
plt.title('PMF of bins')
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











