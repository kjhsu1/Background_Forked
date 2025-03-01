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

'''
parser.add_argument('max_peak_height', type=int, help='max height of peaks')
parser.add_argument('max_multiplier', type=int, help='max multiplier for bias')

n_peaks = args.n_peaks
'''

'''
for h in range(h_peak):
    h += 1
    h = h * multiplier
    print(h)
'''
h = 1

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
    ls_peaks = []
    for i in range(genome_size):
        pos = i+1
        peaks = []
        x = 0   # temp variable name
        while pos <= genome_size:
            try:
                if pos == peaks[-1]:
                    pos += 1
            except IndexError:
                pass
            if len(peaks) == n_peaks:
                if x == 0:
                    x = peaks[-1]
                if peaks[-1] == genome_size:
                    for n in range(n_peaks-1):
                        peaks.pop()
                    pos = x
                    x = 0
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