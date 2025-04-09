import argparse
import random

parser = argparse.ArgumentParser(description='Generates peaks')
parser.add_argument('gfile', type=str, help='file containing genome')
parser.add_argument('--sbp', type=int, default=1000, help='space between peaks')
parser.add_argument('--peaktype', type=str, default='narrow', help='narrow or wide peaks')
parser.add_argument('--seed', type=int, default=1, help='set random seed')
args = parser.parse_args()

random.seed(1)

def get_genome(genome_file):
    with open(genome_file) as fp:
        genome = []
        for line in fp:
            line = line.strip('\n')
            genome.append(line)
        genome = ''.join(genome)
    return genome

genome = get_genome(args.gfile)
sbp = args.sbp # space between peaks

def make_peaks(genome,sbp):
    ls_peaks = []
    peak = True
    size = 0
    count = 0
    counter = 0 # need better name
    for n in range(len(genome)):
        if count == sbp:
            peak = True
        if counter != size:
            counter += 1
            continue
        if peak == False:
            count += 1
            continue
        size = random.randint(50,500)
        ls_peaks.append(genome[n:n+size])
        print(size, n, count, counter)
        peak = False
        count = 0
        counter = 0
        if n > 3: break
    return ls_peaks

peaks = make_peaks(genome,sbp)
print(len(peaks[1]))