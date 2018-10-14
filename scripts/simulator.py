#Charlotte Darby cdarby@jhu.edu
#updated 09-12-18

import argparse
import random
import numpy
import bisect
import sys

def generateHaplotypes(k, numsnp, snprate, fraglengthmax):
	snp_loc = [fraglengthmax-1] #coordinates are 0-indexed
	Z = numpy.random.geometric(p=snprate, size=numsnp-1)
	for z in Z:
		snp_loc.append(snp_loc[-1] + z)

	haplotypes = []
	if k == 2:
		h = numpy.ndarray.tolist(numpy.random.choice([0,1],numsnp))
		haplotypes.append(h) #generate a random haplotype
		haplotypes.append([abs(x-1) for x in h]) #and complement
	else:
		for _ in xrange(k):
			haplotypes.append([])
		for _ in xrange(numsnp): #for each site, there is at least one 0 and at least one 1, but others are randomly selected
			alleles = numpy.ndarray.tolist(numpy.random.choice([0,1],k-2)) + [0,1]
			random.shuffle(alleles)
			for i in xrange(len(haplotypes)):
				haplotypes[i].append(alleles[i])

	return (snp_loc, haplotypes)


def generateFragments(snp_loc, haplotypes, nreads, k, distribution, out):
	if distribution == "normal":
		lengths = numpy.random.normal(loc=mean, scale=sd, size=nreads)
	elif distribution == "lognormal":
		lengths = numpy.random.lognormal(mean=mean, sigma=sd, size=nreads)
	haps = numpy.random.choice(range(k), nreads)

	with open(out, "w") as F:
		for (L,h) in zip(lengths, haps):
			start = random.randrange(genome_len)
			end = start + int(max(min(L, len_max), len_min))
			start_snp = bisect.bisect(snp_loc, start)
			end_snp = bisect.bisect(snp_loc, end)
			if end_snp > start_snp + 1: #covers at least 2 SNP
				F.write(str(start_snp))
				F.write("\t")
				F.write("".join([str(x) for x in haplotypes[h][start_snp:end_snp]]))
				F.write("\n")


def generateFragmentsPaired(snp_loc, haplotypes, nreads, k, readlength, out):
	lengths = numpy.random.normal(loc=mean, scale=sd, size=nreads)
	haps = numpy.random.choice(range(k), nreads)

	with open(out, "w") as F:
		for (L,h) in zip(lengths, haps):
			start = random.randrange(genome_len)
			end = start + int(max(min(L, len_max), len_min))
			f = []
			start_snp = bisect.bisect(snp_loc, start)
			first_snp = start_snp
			end_snp = bisect.bisect(snp_loc, start+readlength)
			if start_snp != end_snp:
				f += haplotypes[h][start_snp:end_snp]

			start_snp = bisect.bisect(snp_loc, end-readlength)
			skip = start_snp - end_snp
			end_snp = bisect.bisect(snp_loc, end)
			if start_snp != end_snp:
				# add dashes for gaps in between mates if there are snps on both sides
				if len(f) > 0 and skip > 0: f.append("-"*skip) 
				f += haplotypes[h][start_snp:end_snp]
			if len(f) > 1: 
				F.write(str(first_snp))
				F.write("\t")
				F.write("".join([str(x) for x in f]))
				F.write("\n")


parser = argparse.ArgumentParser(description='Simulate PAIRED OR SINGLE END READ haplotype-informative fragments')
parser.add_argument('--matrixout', help='output filename for the fragment matrix',required=True)
parser.add_argument('--hapout', help='output filename for the actual haplotype(s) - depending on error rate, may not be the haplotypes corresponding to the MEC score',required=True)
parser.add_argument('--paired', help='Include this argument if the reads are paired-end',action='store_true') #default False

parser.add_argument('--cov',help="sequencing coverage - depending on SNP density, fragment coverage will be less",required=False,default=10)
parser.add_argument('--numsnp',help="number of SNP sites",required=False,default=100)
parser.add_argument('--snprate',help="probability that a site is an SNP",required=False,default=0.01)
parser.add_argument('--errorrate',help="Sequencing error rate",required=False,default=0.02)

parser.add_argument('--fraglengthmean',help="length of the complete sequencing fragment - mean",required=False,default=550)
parser.add_argument('--fraglengthsd',help="length of the complete sequencing fragment - stdev",required=False,default=30)
parser.add_argument('--fraglengthmin',help="min sequencing fragment length",required=False,default=500)
parser.add_argument('--fraglengthmax',help="max sequencing fragment length",required=False,default=600)
parser.add_argument('--readlength',help="length of each mate in pair IF PAIRED",required=False,default=150)
parser.add_argument('--distribution',help="statistical distribution for read length IF SINGLE: normal, lognormal",required=False,default="normal")
parser.add_argument('--k',help="k>1 number of haplotypes to simulate",required=False,default=2)

args = parser.parse_args()

cov = int(args.cov)
numsnp = int(args.numsnp)
snprate = float(args.snprate)
errorrate = float(args.errorrate)
mean = int(args.fraglengthmean)
readlength = int(args.readlength)

sd = int(args.fraglengthsd)
len_min = int(args.fraglengthmin)
len_max = int(args.fraglengthmax)

k = int(args.k)

# validate parameters
if cov <= 0 or len_min <= 0 or len_max <= 0 or mean <= 0:
	print("Coverage and length parameters must be positive")
	sys.exit()
if len_min > len_max:
	print("Minimum fragment length must be less than or equal to maximum fragment length")
	sys.exit()
if numsnp < 2:
	print("Must be at least two SNPs")
	sys.exit
if sd < 0: 
	print("Stdev must be nonnegative")
	sys.exit()
if k < 2:
	print "k must be greater at least 2"
	sys.exit()
if errorrate < 0.0 or errorrate > 1.0:
	print "Error rate must be between 0 and 1"
	sys.exit()
if args.paired and readlength <= 0:
	print("Read length parameter must be positive")
	sys.exit()
if args.distribution not in ["normal", "lognormal"]:
	print("Invalid distribution")
	sys.exit()

(snp_loc, haplotypes) = generateHaplotypes(k, numsnp, snprate, len_max)

with open(args.hapout,"w") as F:
	for h in haplotypes:
		F.write("".join([str(x) for x in h]))
		F.write("\n")

genome_len = snp_loc[-1] + 1
nreads = int(genome_len * cov / 2 / mean)
print(nreads)

if args.paired:
	generateFragmentsPaired(snp_loc, haplotypes, nreads, k, readlength, args.matrixout)
else: 
	generateFragments(snp_loc, haplotypes, nreads, k, args.distribution, args.matrixout)
