
#Charlotte Darby cdarby@jhu.edu
#updated 11-16-17

import argparse
import random,numpy,math

def diploid(fragmentlength,height, width, gaprate, errorrate,k,frach0,out):
	haplotype = numpy.random.choice([0,1],width) #numpy.ndarray.tolist
	with open(args.matrixout,"w") as out:

		if fragmentlength == None or int(fragmentlength) == width: #No gaps at the ends of fragments
			for _ in xrange(int(height*fractionh0)): #h0
				fragment = ["-" if random.random() < gaprate else str(haplotype[i]) if random.random() > errorrate else str(abs(1-haplotype[i])) for i in xrange(width)]
				if args.shortfmt: out.write("0\t" + "".join(fragment) + "\n")
				else: out.write("".join(fragment) + "\n")
				

			for _ in xrange(height-int(height*fractionh0)): #h1 = opposite of h0; change to < errorrate
				fragment = ["-" if random.random() < gaprate else str(haplotype[i]) if random.random() < errorrate else str(abs(1-haplotype[i])) for i in xrange(width)]
				if args.shortfmt: out.write("0\t" + "".join(fragment) + "\n")
				else: out.write("".join(fragment) + "\n")

		else:
			fragmentlength = int(fragmentlength) #There are gaps at the ends of fragments
			end = width-fragmentlength+1 #index where rightmost fragment could start +1 [OR] number of fragment starting points 
			
			startpoints = []
			persite = height/end #rounds down
			for _ in xrange(persite):
				startpoints += range(0,end)

			remaining = height - persite * end
			#assert remaining < end
			if remaining > 0:
				if 1.0*end/remaining < 2: #choose places where fragments won't start instead of places where they will
					step = end/(end-remaining) #number of places to choose
					A = range(0,end,step)
					startpoints += [i for i in xrange(end) if i not in A][:remaining] #"complement" of the sequence 
				else: #choose places where fragments will start
					step = end/remaining
					startpoints += range(0,end,step)[:remaining]
			#assert len(startpoints) == height == remaining + persite * end
			random.shuffle(startpoints)
			H0 = startpoints[0::2]
			H1 = startpoints[1::2]
			for i in H0: #number of gaps to put in front
				prefix = ["-" for _ in xrange(i)] 
				fragment = [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() > errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)] 
				suffix = ["-" for _ in xrange(i+fragmentlength,width)]
				if args.shortfmt: out.write(str(i) + "\t" + "".join(fragment) + "\n")
				else: out.write("".join(prefix) + "".join(fragment) + "".join(suffix) + "\n")

			for i in H1: #number of gaps to put in front
				prefix = ["-" for _ in xrange(i)] 
				fragment = [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() < errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)]  #h1 = opposite of h0; change to < errorrate
				suffix = ["-" for _ in xrange(i+fragmentlength,width)]
				if args.shortfmt: out.write(str(i) + "\t" + "".join(fragment) + "\n")
				else: out.write("".join(prefix) + "".join(fragment) + "".join(suffix) + "\n")
		return [haplotype]

def kploid(fragmentlength,height, width, gaprate, errorrate,k,frach0,out):
	haplotypes = []
	for _ in xrange(k):
		haplotypes.append([])
	for _ in xrange(width): #for each site, there is at least one 0 and at least one 1, but others are randomly selected
		alleles = numpy.ndarray.tolist(numpy.random.choice([0,1],k-2)) + [0,1]
		random.shuffle(alleles)
		for i in xrange(len(haplotypes)):
			haplotypes[i].append(alleles[i])	
	with open(args.matrixout,"w") as out:

		if fragmentlength == None or int(fragmentlength) == width: #No gaps at the ends of fragments
			for _ in xrange(int(height*fractionh0)): #h0
				fragment = ["-" if random.random() < gaprate else str(haplotype[i]) if random.random() > errorrate else str(abs(1-haplotype[i])) for i in xrange(width)]
				if args.shortfmt: out.write("0\t" + "".join(fragment) + "\n")
				else: out.write("".join(fragment) + "\n")
				

			for _ in xrange(height-int(height*fractionh0)): #h1 = opposite of h0; change to < errorrate
				fragment = ["-" if random.random() < gaprate else str(haplotype[i]) if random.random() < errorrate else str(abs(1-haplotype[i])) for i in xrange(width)]
				if args.shortfmt: out.write("0\t" + "".join(fragment) + "\n")
				else: out.write("".join(fragment) + "\n")

		else:
			fragmentlength = int(fragmentlength) #There are gaps at the ends of fragments
			end = width-fragmentlength+1 #index where rightmost fragment could start +1
			
			startpoints = []
			persite = height/end #rounds down
			for _ in xrange(persite):
				startpoints += range(0,end)

			remaining = height - persite * end
						#assert remaining < end
			if remaining > 0:
				if 1.0*end/remaining < 2: #choose places where fragments won't start instead of places where they will
					step = end/(end-remaining) #number of places to choose
					A = range(0,end,step)
					startpoints += [i for i in xrange(end) if i not in A][:remaining] #"complement" of the sequence 
				else: #choose places where fragments will start
					step = end/remaining
					startpoints += range(0,end,step)[:remaining]
			#assert len(startpoints) == height == remaining + persite * end
			random.shuffle(startpoints)
			allstartpoints = [startpoints[0::i] for i in xrange(k)]
			for (S, haplotype) in zip(allstartpoints,haplotypes):
				for i in S: #number of gaps to put in front
					prefix = ["-" for _ in xrange(i)] 
					fragment = [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() > errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)] 
					suffix = ["-" for _ in xrange(i+fragmentlength,width)]
					if args.shortfmt: out.write(str(i) + "".join(fragment) + "\n")
					else: out.write("".join(prefix) + "".join(fragment) + "".join(suffix) + "\n")

		return haplotypes

parser = argparse.ArgumentParser(description='simulateDiploidMatrix.py --matrixout --hapout [--height --width --gaprate --errorrate --fractionh0 --fragmentlength --shortfmt]')
parser.add_argument('--matrixout', help='output filename for the fragment matrix',required=True)
parser.add_argument('--hapout', help='output filename for the actual haplotype - depending on error rate, may not be the haplotypes corresponding to the MEC score',required=True)
parser.add_argument('--height',help="fragment matrix height (number of fragments)",required=False,default=2)
parser.add_argument('--width',help="fragment matrix width (number of sites)",required=False,default=1)
parser.add_argument('--gaprate',help="fraction of entries within fragments that are gaps (randomly distributed; can be at ends of fragments)",required=False,default=0.0)
parser.add_argument('--errorrate',help="fraction of matrix entries that have errors (randomly distributed) given that there is not a gap at that point",required=False,default=0.0)
parser.add_argument('--fractionh0',help="fraction of rows that come from 'hap 0' - NOT USED",required=False,default=0.5)
parser.add_argument('--fragmentlength',help="length of fragments between first non-gap and last non-gap",required=False,default=None)
parser.add_argument('--shortfmt',help="Write the fragment matrix in abbreviated format?",required=False,default=False)
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
		print "Since you are using k=2, fractionh0 must be between 0 and 1 (default is 0.5) [this param is NOT USED]"
		sys.exit()
	H = kploid(args.fragmentlength,height, width, gaprate, errorrate,k,args.matrixout)

with open(args.hapout,"w") as F:
	for h in H:
		F.write("".join([str(h[i]) for i in xrange(width)]))
		F.write("\n")