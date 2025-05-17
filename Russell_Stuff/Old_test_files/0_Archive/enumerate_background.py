import argparse
import sys

parser = argparse.ArgumentParser(description='creates list of positions and its bias peaks')
parser.add_argument('genome_size', type=int, help='size of genome')
parser.add_argument('n_peaks', type=int, help='number of peaks in genome (has to be >= genome size)')
# parser.add_argument('h_peak', type=int, help='height of peaks')
# parser.add_argument('multiplier', type=int, help='multiplier for background')
args = parser.parse_args()

genome_size = args.genome_size
n_peaks = args.n_peaks
# h_peak = args.h_peak
# multiplier = args.multiplier

if genome_size < n_peaks:
    print('The number of peaks cannot be greater than the genome size')
    sys.exit()
if genome_size <= 0:
    print("Genome size has to an integer greater than 0")
    sys.exit()

'''
for h in range(h_peak):
    h += 1
    h = h * multiplier
    print(h)

h = 1
'''

def two_or_less(genome_size, n_peaks):
    ls_peaks = []
    if n_peaks == 1:
        for n in range(genome_size):
            ls_peaks.append([n+1])
    if n_peaks == 2:
        for n in range(genome_size):
            for m in range(genome_size):
                if n != m and n < m: 
                    peaks = [n+1,m+1]
                    ls_peaks.append(list(peaks))
    return ls_peaks

def greater_than_2(genome_size, n_peaks):
    ls_peaks = []   # final list of all pos of peaks in a each possiblity
    for i in range(genome_size):
        pos = i+1
        peaks = []  # temp list containing pos of peaks for current possibility
        last_num = 0   # stores the last number of first set of each peaks list
        
        # keep looping until the pos is equal to the genome size
        while pos <= genome_size:
            # stops the empty list of peaks from giving an error
            try:
                if pos == peaks[-1]:    # makes sure the number added to list is not a repeat
                    pos += 1
            except IndexError:
                pass
            if len(peaks) == n_peaks:
                
                if last_num == 0:   # assigns the last number to a variable of the first list
                    last_num = peaks[-1]
                # replaces pos with the last_num and gets rid of the last 2 indexes
                if peaks[-1] == genome_size:
                    for n in range(n_peaks-1):
                        peaks.pop()
                    pos = last_num
                    last_num = 0
                # gets rid of last index
                else:
                    peaks.pop()
                    peaks.append(pos)

            else:
                peaks.append(pos)
                pos += 1

            if len(peaks) == n_peaks:
                ls_peaks.append(list(peaks))
            
    return ls_peaks

if n_peaks > 2: print(greater_than_2(genome_size,n_peaks))
else: print(two_or_less(genome_size, n_peaks))