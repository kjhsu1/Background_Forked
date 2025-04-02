import argparse

parser = argparse.ArgumentParser(description='creates background multiplier')
parser.add_argument('--file', type=str, help='genome file')
parser.add_argument('-o','--out', type=str, default='multiplier.txt', help='output filename')
args = parser.parse_args()

genome_file = args.file

# extracts genome from file and returns a string
def get_genome():
    with open(genome_file) as fp:
        genome = []
        for line in fp:
            line = line.strip('\n')
            genome.append(line)
        genome = ''.join(genome)
    return genome

genome = get_genome()
rov = 7 # size of regions of variability

# .001x .25x .5x 1x 2x 4x 1000x
# ASCII A = logbase2(0)
multipliers = ['7','?','@','A','B','C','K']

# assigns a multiplier to a region
# multiplier values are in a predetermined order
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
    count = 0 # counts the characters printed
    pos = 0 # used for desinating which value is printed
    counter = 0 # need better name, counts the times a multiplier symbol is printed
    fp.write('@genome\n') # line that follows '@' is the genome
    for nuc in genome:
        fp.write(nuc)
        count += 1
        if count == 50: # 50 because wrap at that many characters
            fp.write('\n')
            if pos == 0: fp.write('+values\n') # line that follows '+' is the value
            else:        fp.write('+\n')
            for n in range(50):
                fp.write(rm[pos])
                counter += 1
                if counter == rov:
                    pos += 1
                    counter = 0
            count = 0
            fp.write('\n')
            if pos != len(rm)-1: fp.write('@\n')

'''
Notes:
    Only can take in genomes in file format right now
        Not sure if it is necessary to have different genome input methods
'''