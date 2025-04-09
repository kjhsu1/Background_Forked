import argparse
import random
import background_maker as bm

parser = argparse.ArgumentParser(description='randomly samples genome')
parser.add_argument('size', type=int, help='size of genome')
parser.add_argument('coverage', type=int, help='depth of coverage')
parser.add_argument('--seed', type=int, help='seed for random')
args = parser.parse_args()

size = args.size
coverage = args.coverage
random.seed(args.seed)

def genome_array(size):
    genome = []
    for n in range(size):
        genome.append(0)
    return genome
genome = genome_array(size)


bm.giza(genome,10,)



def sample(genome,coverage):
    for n in range(coverage):
        for pos in range(len(genome)):
            genome[pos] += random.randint(0,1)
    return genome
genome = sample(genome,coverage)

