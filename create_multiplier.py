import argparse

parser = argparse.ArgumentParser(description='creates background multiplier')
parser.add_argument('--file', type=str, help='genome file')
parser.add_argument('-o','--out', type=str, default='multiplier.txt', help='output filename')
args = parser.parse_args()

def get_genome():
    with open(args.file) as fp:
        genome = []
        for line in fp:
            genome.append(line)
        ''.join(genome)
        genome = genome[0]
    return genome

genome = get_genome()
rov = 7 # number of regions of variability

# .001x .25x .5x 1x 2x 4x 1000x
multipliers = ['7','?','@','A','B','C','K'] # ASCII A = dec 0

def assigning_multiplier(genome, multipliers, rov):
    count = 0
    regions =[]
    region_multiplier = []
    for region in range(0,len(genome),rov):
        regions.append(genome[region:region+rov])
        if count == len(multipliers):
            count = 0
        region_multiplier.append(multipliers[count])
        count += 1
    return regions, region_multiplier

regions, rm = assigning_multiplier(genome,multipliers,rov)

with open(args.out,'w') as fp:
    count = 0
    for piece in regions:
        fp.write(piece)
    fp.write('\n')
    for part in regions:
        for nuc in part:
            fp.write(f"{rm[count]}")
        count += 1

# print(regions, rm)
# print(len(regions), len(rm))

'''
Notes:
    What should the values of the multiplier be?
    The size of the regions with different multipliers are mostly fixed
    Only can take in genomes in file format right now
        Not sure if it is necessary to have different genome input methods
'''