import random
import argparse

parser = argparse.ArgumentParser(description='creates list of positions and its bias peaks')
parser.add_argument('genome_size', type=int, help='size of genome')
parser.add_argument('n_peaks', type=int, help='number of peaks in genome (has to be >= genome size)')
parser.add_argument('max_peak_height', type=int, help='max height of peaks')
parser.add_argument('n_var_regions', type=int, help='number of regions with bias')
parser.add_argument('max_multiplier', type=int, help='max multiplier for bias')
parser.add_argument('-s', '--seed', type=int, help='sets random seed')
args = parser.parse_args()

# generate peaks
def make_peaks(genome_size, n_peaks, max_peak_height):

    # makes sure number of peaks does not exceed genome size
    if n_peaks > genome_size:
        return "Invalid number of peaks"

    # assigns a pos in the genome for each peak
    pos_peaks = list()
    for pos in range(n_peaks):
        pos = random.randint(1,genome_size)
        while pos in pos_peaks:   # checks to make sure that peaks are not in the same location
            pos = random.randint(1,genome_size)
        pos_peaks.append(pos)
    pos_peaks.sort()
    #print(pos_peaks)    # check

    # assigns a peak height for each peak from 1 to a user set max
    ls_peaks = list()
    for peak in pos_peaks:
        peak = [peak]   # makes peak a list to allow for appending a corresponding height
        peak.append(random.randint(1, max_peak_height)) # assigns random height
        #print(peak) # check to see each peak has a height
        ls_peaks.append(peak)
        #print(f"{ls_peaks}\n") # check list has correct pos and heights
    return ls_peaks # List contains peak positions and heights

#print(make_peaks(genome_size, n_peaks, max_peak_height))

# generating multiplier
def make_multiplier(genome_size, n_var_regions, max_multiplier):
    
    # makes sure number of var regions does not exceed genome size
    if n_var_regions > genome_size:
        return "Invalid number of variable regions"
    
    ls_regions = list()
    for n in range(n_var_regions):  # generates regions of variability
        pos = random.randint(1, genome_size)
        while pos in ls_regions:    # makes sure variable regions are not overlapping
            pos = random.randint(1, genome_size)
        ls_regions.append(pos)
    ls_regions.sort()
    #print(ls_regions)   # check pos is created correctly

    ls_multiplier = list()
    for pos in ls_regions: # similar to finding peaks
        pos = [pos]
        pos.append(random.randint(1, max_multiplier))
        ls_multiplier.append(pos)
    return ls_multiplier

#print(make_multiplier(genome_size, n_var_regions, max_multiplier))

def make_background(ls_peaks, ls_multiplier):

    background = list()
    for peak in ls_peaks:
        for multiplier in ls_multiplier:
            if peak[0] == multiplier[0]:
                final_peak = peak[1] * multiplier[1]
                #print(peak,multiplier)
                #print(f"{final_peak}\n")
                background.append([peak[0], final_peak])
                break
    return background

random.seed(args.seed)

# gets peaks
ls_peaks = make_peaks(args.genome_size, args.n_peaks, args.max_peak_height)
# get bias
ls_multiplier = make_multiplier(args.genome_size, args.n_var_regions, args.max_multiplier)
# gets background
background = make_background(ls_peaks, ls_multiplier)
print(ls_peaks)
print(ls_multiplier)
print(background)