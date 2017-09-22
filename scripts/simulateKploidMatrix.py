
#Charlotte Darby cdarby@jhu.edu
#updated 9-21-17

#Developed on Python 2.7.13 

import argparse,sys
import random
import numpy

def kploid(fragmentlength,height, width, gaprate, errorrate,k,out):
	haplotypes = []
	for _ in xrange(k):
		haplotypes.append([])
	for _ in xrange(width): #for each site, there is at least one 0 and at least one 1, but others are randomly selected
		alleles = numpy.ndarray.tolist(numpy.random.choice([0,1],k-2)) + [0,1]
		random.shuffle(alleles)
		for i in xrange(len(haplotypes)):
			haplotypes[i].append(alleles[i])

	with open(out,"w") as F:

		if fragmentlength == None or fragmentlength == width: #No gaps at the ends of fragments
			for h in haplotypes:
				for _ in xrange(int(1.0*height/k)):
					fragment = ["-" if random.random() < gaprate else str(h[i]) if random.random() > errorrate else str(abs(1-h[i])) for i in xrange(width)]
					F.write("".join(fragment))
					F.write("\n")

			for h in xrange(height-k*int(1.0*height/k)): #Number of fragments may not be divisible by k, these are the "leftovers"
				fragment = ["-" if random.random() < gaprate else str(haplotypes[h][i]) if random.random() > errorrate else str(abs(1-haplotypes[h][i])) for i in xrange(width)]
				F.write("".join(fragment))
				F.write("\n")

		else:
			fragmentlength = width if int(fragmentlength) > width else int(fragmentlength) #There are gaps at the ends of fragments
			#If fragments cannot be distributed evenly, the coverage difference between any two sites is at most 2
			for h in haplotypes:
				persite = int(1.0*int(1.0*height/k)/(width-fragmentlength+1))
				remaining = int(1.0*height/k) - persite*(width-fragmentlength)
				
				if persite > 0:
					for i in xrange(width-fragmentlength+1): #number of gaps to put in front
						for p in xrange(persite):
							fragment = ["-" for _ in xrange(i)] + [ "-" if random.random() < gaprate else str(h[ch]) if random.random() > errorrate else str(abs(1-h[ch])) for ch in xrange(i,i+fragmentlength)] + ["-" for _ in xrange(i+fragmentlength,width)]
							F.write("".join(fragment))
							F.write("\n")
				for i in random.sample(range(width-fragmentlength+1),remaining): #number of gaps to put in front
					fragment = ["-" for _ in xrange(i)] + [ "-" if random.random() < gaprate else str(h[ch]) if random.random() > errorrate else str(abs(1-h[ch])) for ch in xrange(i,i+fragmentlength)] + ["-" for _ in xrange(i+fragmentlength,width)]
					F.write("".join(fragment))
					F.write("\n")
			
			remaining = height-k*int(1.0*height/k)
			for (i,h) in zip(random.sample(range(width-fragmentlength+1),remaining),haplotypes): #number of gaps to put in front	
				fragment = ["-" for _ in xrange(i)] + [ "-" if random.random() < gaprate else str(h[ch]) if random.random() < errorrate else str(abs(1-h[ch])) for ch in xrange(i,i+fragmentlength)] + ["-" for _ in xrange(i+fragmentlength,width)]
				F.write("".join(fragment))
				F.write("\n")
	return haplotypes

def diploid(fragmentlength,height, width, gaprate, errorrate,k,frach0,out):
	haplotype = numpy.random.choice([0,1],width) #numpy.ndarray.tolist
	with open(out,"w") as F:

		if fragmentlength == None or fragmentlength == width: #No gaps at the ends of fragments
			for _ in xrange(int(height*fractionh0)): #h0
				fragment = ["-" if random.random() < gaprate else str(haplotype[i]) if random.random() > errorrate else str(abs(1-haplotype[i])) for i in xrange(width)]
				F.write("".join(fragment))
				F.write("\n")
			for _ in xrange(height-int(height*fractionh0)): #h1 = opposite of h0; change to < errorrate
				fragment = ["-" if random.random() < gaprate else str(haplotype[i]) if random.random() < errorrate else str(abs(1-haplotype[i])) for i in xrange(width)]
				F.write("".join(fragment))
				F.write("\n")

		else:
			fragmentlength = width if int(fragmentlength) > width else int(fragmentlength) #There are gaps at the ends of fragments
			#If fragments cannot be distributed evenly, the coverage difference between any two sites is at most 2
			
			persite = int(1.0*int(height*fractionh0)/(width-fragmentlength+1))
			remaining = int(height*fractionh0) - persite * (width-fragmentlength+1)

			if persite > 0:
				for i in xrange(width-fragmentlength+1): #number of gaps to put in front
					for p in xrange(persite):
						fragment = ["-" for _ in xrange(i)] + [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() > errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)] + ["-" for _ in xrange(i+fragmentlength,width)]
						F.write("".join(fragment))
						F.write("\n")
			for i in random.sample(range(width-fragmentlength+1),remaining): #number of gaps to put in front
				fragment = ["-" for _ in xrange(i)] + [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() > errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)] + ["-" for _ in xrange(i+fragmentlength,width)]
				F.write("".join(fragment))
				F.write("\n")
			
			persite = int(1.0*(height - int(height*fractionh0))/(width-fragmentlength+1)) #h1 = opposite of h0; change to < errorrate
			remaining = height - int(height*fractionh0) - persite * (width-fragmentlength+1)

			if persite > 0:
				for i in xrange(width-fragmentlength+1): #number of gaps to put in front
					for p in xrange(persite):
						fragment = ["-" for _ in xrange(i)] + [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() < errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)] + ["-" for _ in xrange(i+fragmentlength,width)]
						F.write("".join(fragment))
						F.write("\n")
			for i in random.sample(range(width-fragmentlength+1),remaining): #number of gaps to put in front
				fragment = ["-" for _ in xrange(i)] + [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() < errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)] + ["-" for _ in xrange(i+fragmentlength,width)]
				F.write("".join(fragment))
				F.write("\n")
	return [haplotype]

parser = argparse.ArgumentParser(description='simulateKploidMatrix.py --matrixout --hapout [--height --width --gaprate --errorrate --fractionh0 --fragmentlength --k]')
parser.add_argument('--matrixout', help='output filename for the fragment matrix',required=True)
parser.add_argument('--hapout', help='output filename for the actual haplotype - depending on error rate, may not be the haplotypes corresponding to the MEC score',required=True)
parser.add_argument('--height',help="fragment matrix height (number of fragments)",required=False,default=2)
parser.add_argument('--width',help="fragment matrix width (number of sites)",required=False,default=1)
parser.add_argument('--gaprate',help="fraction of entries within fragments that are gaps (randomly distributed)",required=False,default=0.0)
parser.add_argument('--errorrate',help="fraction of matrix entries that have errors (randomly distributed) given that there is not a gap at that point",required=False,default=0.0)
parser.add_argument('--fractionh0',help="fraction of rows that come from 'hap 0' (only applies to k=2)",required=False,default=0.5)
#model unevenness of haplotype in k-ploid later?
parser.add_argument('--fragmentlength',help="length of fragments between first non-gap and last non-gap",required=False,default=None)
parser.add_argument('--k',help="k>1 number of haplotypes to simulate",required=False,default=2)


args = parser.parse_args()
(height,width) = (int(args.height),int(args.width))
gaprate = float(args.gaprate)
errorrate = float(args.errorrate)
k = int(args.k)
fractionh0 = float(args.fractionh0) if args.fractionh0 != None else None

if args.fragmentlength != None and (int(args.fragmentlength) < 1 or int(args.fragmentlength) > width):
	print "Since you specified fragment length, it must be at least 1 and at most width (default is width)"
	sys.exit()
if height < 1 or width < 1:
	print "dimensions must be at least 1"
	sys.exit()
if k < 2:
	print "k must be greater or equal to 2"
	sys.exit()
if errorrate < 0.0 or errorrate > 1.0 or gaprate < 0.0 or gaprate > 1.0:
	print "rates must be between 0 and 1"
	sys.exit()
if k == 2:
	H = diploid(args.fragmentlength,height, width, gaprate, errorrate,k,fractionh0,args.matrixout)
else:
	if fractionh0 != None and (fractionh0 < 0.0 or fractionh0 > 1.0):
		print "Since you are using k=2, fractionh0 must be between 0 and 1 (default is 0.5)"
		sys.exit()
	H = kploid(args.fragmentlength,height, width, gaprate, errorrate,k,args.matrixout)

with open(args.hapout,"w") as F:
	for h in H:
		F.write("".join([str(h[i]) for i in xrange(width)]))
		F.write("\n")
