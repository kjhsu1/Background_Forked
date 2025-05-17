import sys
import math

def genome(genome_size):
    genome = []
    for pos in range(genome_size):
        genome.append(0)
    return genome
    

def mod_genome(genome, n_peaks):
    APC = 0
    n_APC = math.comb(genome_size,n_peaks)
    while APC != n_APC:
        peaks = 0
        new_genome = list(genome)
        pos = 0
        while peaks != n_peaks:
            new_genome[pos] = 1
            pos += 1
            new_genome[pos] = 1
            peaks += 1
            print(new_genome)
            while True:
                new_genome[pos] = 0
                pos += 1
                new_genome[pos] = 1
                print(new_genome,pos)
                if pos == n_peaks+1: break
        APC += 1
        print()
    return

genome_size = int(sys.argv[1])
n_peaks = int(sys.argv[2])
genome = genome(genome_size)
test = mod_genome(genome,n_peaks)
test

