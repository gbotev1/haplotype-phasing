#!/usr/bin/env python3
# Author: Georgie Botev

import sys
from argparse import ArgumentParser


def read_in():
    # Read in matrix from command-line
    # The first k lines are the k reconstructed haplotypes
    # The next k lines are the k real haplotypes
    # The real and reconstructed haplotypes should be separated
    # by an extra linebreak.
    parser = ArgumentParser(
        description='A Python3 script to calculate the vector error.')
    parser.add_argument('input',
                        nargs='*',
                        help='The input matrix in the required form.',
                        type=str)
    args = vars(parser.parse_args())
    found_haps, _, real_haps = ' '.join(args['input']).partition('\\n\\n')
    found_haps = found_haps.split('\\n')
    found_haps = [x.split(' ') for x in found_haps]
    real_haps = real_haps.split('\\n')
    real_haps = [x.split(' ') for x in real_haps]
    return found_haps, real_haps


# Question for later:
# What to do if we have a gap?
def d(hap1, hap2):
    num_diffs = 0
    for col in range(len(hap1)):
        if hap1[col] != hap2[col]:
            num_diffs += 1
    return num_diffs


def main():
    found_haps, real_haps = read_in()
    VE = {}
    k = len(real_haps)
    n = len(real_haps[0])
    # Initialization step
    for hap1 in range(k):
        if found_haps[hap1][0] in [x[0] for x in real_haps]:
            VE[hap1] = 0
    # Extension step
    for col in range(n):
        for hap1 in range(k):
            if found_haps[hap1][col] in [x[col] for x in real_haps]:
                VE[hap1] = col * k
                for hap2 in range(k):
                    if found_haps[hap2][col - 1] in [x[col - 1] for x in real_haps]:
                        if VE[hap1] > VE[hap2] + d(found_haps[hap1], found_haps[hap2]):
                            VE[hap1] = VE[hap2] + d(found_haps[hap1], found_haps[hap2])
    VErrors = k * n
    for hap1 in range(k):
        if found_haps[hap1][-1] in [x[-1] for x in real_haps]:
            if VErrors > VE[hap1]:
                VErrors = VE[hap1]
    print(VErrors)


if __name__ == "__main__":
    main()
