import argparse

parser = argparse.ArgumentParser(description='creates list of positions and its bias peaks')
parser.add_argument('n_peaks', type=int, help='number of peaks in genome (has to be >= genome size)')
parser.add_argument('genome_size', type=int, help='size of genome')
parser.add_argument('h_peak', type=int, help='height of peaks')
parser.add_argument('multiplier', type=int, help='multiplier for background')
args = parser.parse_args()

n_peaks = args.n_peaks
genome_size = args.genome_size
h_peak = args.h_peak
multiplier = args.multiplier

'''
parser.add_argument('max_peak_height', type=int, help='max height of peaks')
parser.add_argument('max_multiplier', type=int, help='max multiplier for bias')

n_peaks = args.n_peaks
'''

ls_peaks = list()
for h in range(h_peak):
    h += 1
    h = h * multiplier
    print(h)

for pos in range(genome_size):
    for i in range(n_peaks):
        pos += 1
        peak = (pos,h)
        ls_peaks.append(peak)
    break
#print(ls_peaks)