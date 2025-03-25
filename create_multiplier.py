import argparse

parser = argparse.ArgumentParser(description='creates background multiplier')
parser.add_argument('--file', type=str, help='genome file')
parser.add_argument('--genome', type=str, help='genome string (used for testing)')
args = parser.parse_args()

def genome():
    if args.file:
        with open(args.file) as fp:
            g_file = []
            for line in fp:
                g_file.append(line)
            ''.join(g_file)
    print(g_file)

genome()