import random
import argparse

parser = argparse.ArgumentParser(description='creates random genome')
parser.add_argument('genome_size', type=int, help='size of genome')
parser.add_argument('-s','--seed', type=int, help='set random seed')
parser.add_argument('-o','--out', type=str, default='genome.txt', help='output filename')
args = parser.parse_args()

genome_size = args.genome_size

random.seed(args.seed)

def make_genome(size):
    count = 0
    genome = []
    for pos in range(size):
        genome.append(random.choice('ACGT'))
        count += 1
        if count == 50:
            genome.append('\n')
            count = 0
    genome = ''.join(genome)
    return genome

with open(args.out, 'w') as fp:
    fp.write(make_genome(genome_size))

'''
Note:
    The genome wraps every 50 lines
    This genome also has little to no actual properties of a real genome.
'''