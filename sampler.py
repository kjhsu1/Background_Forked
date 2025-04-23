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

completion = coverage/100

def genome_array(size):
    genome = []
    for n in range(size):
        genome.append(0)
    return genome
genome = genome_array(size)
# print(genome)

background = bm.giza(genome)
# print(background)
# print(genome)

def sim_sample(genome,coverage):
    simple_sample = genome.copy()
    for n in range(coverage):
        for pos in range(len(simple_sample)):
            simple_sample[pos] += random.randint(0,1)
    return simple_sample

def sampler_1(genome,coverage,background):
    genome_copy_1 = genome.copy()
    count = 0
    percentage = 0
    for n in range(coverage):
        for pos in range(len(genome_copy_1)):
            if background[pos] >= 1:
                for i in range(background[pos]):
                    genome_copy_1[pos] += random.randint(0,1)
            elif background[pos] < 1:
                if background[pos] >= random.random():
                    genome_copy_1[pos] += random.randint(0,1)
        count += 1
        if count%completion == 0:
            percentage += 1
            print(percentage,'%')
    return genome_copy_1

def sampler_2(genome,coverage,background):
    genome_copy_2 = genome.copy()
    for n in range(coverage):
        for pos in range(len(genome_copy_2)):
            genome_copy_2[pos] += random.randint(0,1)
    for pos in range(len(background)):
        genome_copy_2[pos] = int(genome_copy_2[pos] * background[pos])
    return genome_copy_2

sim_sampled_bg = sim_sample(genome,coverage)
sampled_1_bg = sampler_1(genome,coverage,background)
sampled_2_bg = sampler_2(genome,coverage,background)


print(background, sim_sampled_bg, sampled_1_bg, sampled_2_bg, sep='\n')

# for pos in range(len(background)):
#     if background[pos] < 1:
#         if background[pos] >= random.random():
#             print(background[pos], random.random())

'''
Notes:
Is a proability of .75 considered double the probablity becuase peaks and no
peaks is a coin flip. So twice the coverage is flipping two coins and getting
at least one head?
'''