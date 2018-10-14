import sys
import re

### BEGIN DEFINITIONS

def count_errors_diploid(start, hap1, hap2):
    # Initialize counters for different kinds of errors
    switch = 0
    flip = 0

    # Initialize comparison string
    xorHap = ""

    # Set size equal to length of hap1 (short encoding)
    size = len(hap1)

    # Create comparison string
    for i in xrange(size):
        if hap1[i] == hap2[start + i]:
            xorHap += "0"
        else:
            xorHap += "1"

    # Process comparison string
    i = 0
    while i < size:
        if (int(xorHap[i]) == 1):
            flag = False
            while i < size - 1 and xorHap[i] == xorHap[i + 1]:
                i += 1
                flag = True
            if flag:
                switch += 1
            else:
                flip += 1
        i += 1

    return switch + flip

### BEGIN SCRIPT

# Read in the two haplotypes to compare
with open("h1.txt", "r") as infile:
    real = infile.readlines()
with open("h2.txt", "r") as infile:
    found = infile.readlines()

# Keep track of best pairing seen so far
bestPair = [sys.maxint, ""]

for i in xrange(len(found)):
    hap1 = re.split("\t", found[i])
    for j in xrange(len(real)):
        hap2 = real[j]
        currScore = count_errors_diploid(int(hap1[0]), hap1[1], hap2)
        if currScore < bestPair[0]:
            bestPair[0] = currScore
            bestPair[1] = hap2
    print("%d : %s") % (bestPair[0], bestPair[1])
    bestPair[0] = sys.maxint
