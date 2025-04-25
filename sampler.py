import argparse
import random
import background_maker as bm

parser = argparse.ArgumentParser(description='randomly samples genome')
parser.add_argument('size', type=int, help='size of genome')
parser.add_argument('coverage', type=int, help='depth of coverage')
parser.add_argument('--seed', type=int, help='seed for random')
parser.add_argument('-ph', '--peakheight', type=int, default=16, help='height of peak')
parser.add_argument('-pw', '--peakwidth', type=int, default=3, help='width of peak')
parser.add_argument('--sbp', type=int, default=5, help='space between peaks')
args = parser.parse_args()

size = args.size
coverage = args.coverage
random.seed(args.seed)
ph = args.peakheight
pw = args.peakwidth
sbp = args.sbp


def genome_array(size):
    genome = []
    for n in range(size):
        genome.append(0)
    return genome
genome = genome_array(size)
# print(genome)

background = bm.giza(genome)
# print(background)

def sim_bgsampler(genome,coverage):
    simple_sample = genome.copy()
    for n in range(coverage):
        for pos in range(len(simple_sample)):
            simple_sample[pos] += random.randint(0,1)
    return simple_sample

def slow_bgsampler(genome,coverage,background):
    genome_copy_1 = genome.copy()
    count = 0

    percentage = 0
    completion = coverage/100


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

def fast_bgsampler(genome,coverage,background):
    genome_copy_2 = genome.copy()
    for n in range(coverage):
        for pos in range(len(genome_copy_2)):
            genome_copy_2[pos] += random.randint(0,1)
    for pos in range(len(background)):
        genome_copy_2[pos] = int(genome_copy_2[pos] * background[pos])
    return genome_copy_2

def peaks(genome,ph,pw,sbp):
    ls_peaks = genome.copy()
    count = 0


    peak_value = [ph]
    half_peak = []
    value = ph

    if pw%2 == 0: pos_in_half = int(pw / 2) - 1
    else:         pos_in_half = int(pw / 2)

    slope = ph / (pos_in_half + 1)
    for i in range(pos_in_half):
        value = int(value - slope)
        if value == 0: peak_value.append(1)
        else:          peak_value.append(value)
    peak_value.sort()
    for i in range(1, len(peak_value)):
        half_peak.append(peak_value[-(i+1)])
    if pw%2 == 0: peak_value.append(ph)
    peak_value.extend(half_peak)
    print(peak_value)


    for pos in range(len(ls_peaks)):
        if count < pw:
            ls_peaks[pos] = peak_value[count]
            count += 1
        elif count < pw+sbp:
            count += 1
        else:
            count = 0
    return ls_peaks

print(peaks(genome,ph,pw,sbp), f'ph = {ph}', f'pw = {pw}', f'sbp = {sbp}', sep='\n')





# sim_sampled_bg = sim_bgsampler(genome,coverage)
# sampled_1_bg = slow_bgsampler(genome,coverage,background)
# sampled_2_bg = fast_bgsampler(genome,coverage,background)

# print(background, sim_sampled_bg, sampled_1_bg, sampled_2_bg, sep='\n')



'''
Notes:

'''