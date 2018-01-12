
#Charlotte Darby cdarby@jhu.edu
#updated 11-16-17

import argparse
import random,numpy,math

def diploidPairedend(fragmentlength,insertlength,height, width, gaprate, errorrate,k,outLong,out):
	n_gaps = 0
	n_errors = 0
	haplotype = numpy.random.choice([0,1],width) #generate a random haplotype
	if outLong != None: outLong = open(outLong,"w")
	with open(args.matrixout,"w") as out:
		if args.longfmt: outLong.write(str(height) + "\n" + str(width) + "\n")
		if insertlength == None: insertlength = 0
		else: insertlength = int(insertlength)
		if fragmentlength == None: fragmentlength = (width-insertlength)/2
		else: fragmentlength = int(fragmentlength) 
		span=2*fragmentlength+insertlength
		#assert span <= width
		end = width-span+1 #index where rightmost fragment could start +1 [OR] number of fragment starting points 
		
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
		if len(startpoints) < remaining: 
			for sp in xrange(remaining-len(startpoints)):
				startpoints.append(sp) 
		H0 = startpoints[0::2]
		H1 = startpoints[1::2]
		for i in H0: #number of gaps to put in front
			prefix = ["-" for _ in xrange(i)] 
			suffix = ["-" for _ in xrange(i+span,width)]
			insert = ["-" for _ in xrange(insertlength)]
			fragment1 = [haplotype[q] for q in xrange(i,i+fragmentlength)] #ints
			fragment2 = [haplotype[q] for q in xrange(i+fragmentlength+insertlength,i+span)] #ints
			for f in xrange(len(fragment1)):
				if random.random() < gaprate:
					fragment1[f] = "-"
					n_gaps += 1
				elif random.random() < errorrate:
					fragment1[f] = str(abs(1-fragment1[f]))
					n_errors += 1
				else:
					fragment1[f] = str(fragment1[f])
			for f in xrange(len(fragment2)):
				if random.random() < gaprate:
					fragment2[f] = "-"
					n_gaps += 1
				elif random.random() < errorrate:
					fragment2[f] = str(abs(1-fragment2[f]))
					n_errors += 1
				else:
					fragment2[f] = str(fragment2[f])

			if args.longfmt: outLong.write("".join(prefix) + "".join(fragment1) + "".join(insert) + "".join(fragment2) + "".join(suffix) + "\n")
			out.write(str(i) + "\t" + "".join(fragment1) + "".join(insert) + "".join(fragment2) + "\n")

		for i in H1: #number of gaps to put in front
			prefix = ["-" for _ in xrange(i)] 
			suffix = ["-" for _ in xrange(i+span,width)]
			insert = ["-" for _ in xrange(insertlength)]
			fragment1 = [abs(1-haplotype[q]) for q in xrange(i,i+fragmentlength)] #ints
			fragment2 = [abs(1-haplotype[q]) for q in xrange(i+fragmentlength+insertlength,i+span)] #ints
			for f in xrange(len(fragment1)):
				if random.random() < gaprate:
					fragment1[f] = "-"
					n_gaps += 1
				elif random.random() < errorrate:
					fragment1[f] = str(abs(1-fragment1[f]))
					n_errors += 1
				else:
					fragment1[f] = str(fragment1[f])
			for f in xrange(len(fragment2)):
				if random.random() < gaprate:
					fragment2[f] = "-"
					n_gaps += 1
				elif random.random() < errorrate:
					fragment2[f] = str(abs(1-fragment2[f]))
					n_errors += 1
				else:
					fragment2[f] = str(fragment2[f])

			if args.longfmt: outLong.write("".join(prefix) + "".join(fragment1) + "".join(insert) + "".join(fragment2) + "".join(suffix) + "\n")
			out.write(str(i) + "\t" + "".join(fragment1) + "".join(insert) + "".join(fragment2) + "\n")
			
		if outLong != None: out.close()
		print str(n_gaps) + " gaps"
		print str(n_errors) + " errors"
		return [haplotype]


def kploidPairedend(fragmentlength,insertlength,height, width, gaprate, errorrate,k,outLong,out):
	n_gaps = 0
	n_errors = 0
	haplotypes = []
	for _ in xrange(k):
		haplotypes.append([])
	for _ in xrange(width): #for each site, there is at least one 0 and at least one 1, but others are randomly selected
		alleles = numpy.ndarray.tolist(numpy.random.choice([0,1],k-2)) + [0,1]
		random.shuffle(alleles)
		for i in xrange(len(haplotypes)):
			haplotypes[i].append(alleles[i])	
	if outLong != None: outLong = open(outLong,"w")
	with open(args.matrixout,"w") as out:
		if args.longfmt: outLong.write(str(height) + "\n" + str(width) + "\n")
		if insertlength == None: insertlength = 0
		else: insertlength = int(insertlength)
		if fragmentlength == None: fragmentlength = (width-insertlength)/2
		else: fragmentlength = int(fragmentlength) 
		span=2*fragmentlength+insertlength
		#assert span <= width
		end = width-span+1 #index where rightmost fragment could start +1 [OR] number of fragment starting points 
		
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
		if len(startpoints) < remaining: 
			for sp in xrange(remaining-len(startpoints)):
				startpoints.append(sp) 
		#assign start points to haplotypes
		allH = []
		for h in xrange(k):
			allH.append(startpoints[h::k])
		
		for h,H in enumerate(allH):
			for i in H: #number of gaps to put in front
				prefix = ["-" for _ in xrange(i)] 
				suffix = ["-" for _ in xrange(i+span,width)]
				insert = ["-" for _ in xrange(insertlength)]
				fragment1 = [haplotypes[h][q] for q in xrange(i,i+fragmentlength)] #ints
				fragment2 = [haplotypes[h][q] for q in xrange(i+fragmentlength+insertlength,i+span)] #ints
				for f in xrange(len(fragment1)):
					if random.random() < gaprate:
						fragment1[f] = "-"
						n_gaps += 1
					elif random.random() < errorrate:
						fragment1[f] = str(abs(1-fragment1[f]))
						n_errors += 1
					else:
						fragment1[f] = str(fragment1[f])
				for f in xrange(len(fragment2)):
					if random.random() < gaprate:
						fragment2[f] = "-"
						n_gaps += 1
					elif random.random() < errorrate:
						fragment2[f] = str(abs(1-fragment2[f]))
						n_errors += 1
					else:
						fragment2[f] = str(fragment2[f])

				if args.longfmt: outLong.write("".join(prefix) + "".join(fragment1) + "".join(insert) + "".join(fragment2) + "".join(suffix) + "\n")
				out.write(str(i) + "\t" + "".join(fragment1) + "".join(insert) + "".join(fragment2) + "\n")


		if outLong != None: out.close()
		print str(n_gaps) + " gaps"
		print str(n_errors) + " errors"
		return haplotypes

def diploid(fragmentlength,height, width, gaprate, errorrate,k,outLong,out):
	n_gaps = 0
	n_errors = 0
	haplotype = numpy.random.choice([0,1],width) #generate a random haplotype
	if outLong != None: outLong = open(outLong,"w")
	with open(args.matrixout,"w") as out:
		if args.longfmt: outLong.write(str(height) + "\n" + str(width) + "\n")
		
		if fragmentlength == None: fragmentlength = width
		else: fragmentlength = int(fragmentlength) 
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
		if len(startpoints) < remaining: 
			for sp in xrange(remaining-len(startpoints)):
				startpoints.append(sp) 
		H0 = startpoints[0::2]
		H1 = startpoints[1::2]
		for i in H0: #number of gaps to put in front
			prefix = ["-" for _ in xrange(i)] 
			#fragment = [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() > errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)] 
			suffix = ["-" for _ in xrange(i+fragmentlength,width)]

			fragment = [haplotype[q] for q in xrange(i,i+fragmentlength)] #ints
			for f in xrange(len(fragment)):
				if random.random() < gaprate:
					fragment[f] = "-"
					n_gaps += 1
				elif random.random() < errorrate:
					fragment[f] = str(abs(1-fragment[f]))
					n_errors += 1
				else:
					fragment[f] = str(fragment[f])

			if args.longfmt: outLong.write("".join(prefix) + "".join(fragment) + "".join(suffix) + "\n")
			out.write(str(i) + "\t" + "".join(fragment) + "\n")

		for i in H1: #number of gaps to put in front
			prefix = ["-" for _ in xrange(i)] 
			#fragment = [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() < errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)]  #h1 = opposite of h0; change to < errorrate
			suffix = ["-" for _ in xrange(i+fragmentlength,width)]

			fragment = [abs(1-haplotype[q]) for q in xrange(i,i+fragmentlength)] #ints
			for f in xrange(len(fragment)):
				if random.random() < gaprate:
					fragment[f] = "-"
					n_gaps += 1
				elif random.random() < errorrate:
					fragment[f] = str(abs(1-fragment[f]))
					n_errors += 1
				else:
					fragment[f] = str(fragment[f])

			if args.longfmt: outLong.write("".join(prefix) + "".join(fragment) + "".join(suffix) + "\n")
			out.write(str(i) + "\t" + "".join(fragment) + "\n")
			
		if outLong != None: out.close()
		print str(n_gaps) + " gaps"
		print str(n_errors) + " errors"
		return [haplotype]


def kploid(fragmentlength,height, width, gaprate, errorrate,k,outLong,out):
	n_gaps = 0
	n_errors = 0
	haplotypes = []
	for _ in xrange(k):
		haplotypes.append([])
	for _ in xrange(width): #for each site, there is at least one 0 and at least one 1, but others are randomly selected
		alleles = numpy.ndarray.tolist(numpy.random.choice([0,1],k-2)) + [0,1]
		random.shuffle(alleles)
		for i in xrange(len(haplotypes)):
			haplotypes[i].append(alleles[i])	
	if outLong != None: outLong = open(outLong,"w")
	with open(args.matrixout,"w") as out:
		if args.longfmt: outLong.write(str(height) + "\n" + str(width) + "\n")
		
		if fragmentlength == None: fragmentlength = width
		else: fragmentlength = int(fragmentlength) 
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
		if len(startpoints) < remaining: 
			for sp in xrange(remaining-len(startpoints)):
				startpoints.append(sp) 
		#assign start points to haplotypes
		allH = []
		for h in xrange(k):
			allH.append(startpoints[h::k])
		
		for h,H in enumerate(allH):
			for i in H: #number of gaps to put in front
				prefix = ["-" for _ in xrange(i)] 
				#fragment = [ "-" if random.random() < gaprate else str(haplotype[ch]) if random.random() > errorrate else str(abs(1-haplotype[ch])) for ch in xrange(i,i+fragmentlength)] 
				suffix = ["-" for _ in xrange(i+fragmentlength,width)]

				fragment = [haplotypes[h][q] for q in xrange(i,i+fragmentlength)] #ints
				for f in xrange(len(fragment)):
					if random.random() < gaprate:
						fragment[f] = "-"
						n_gaps += 1
					elif random.random() < errorrate:
						fragment[f] = str(abs(1-fragment[f]))
						n_errors += 1
					else:
						fragment[f] = str(fragment[f])

				if args.longfmt: outLong.write("".join(prefix) + "".join(fragment) + "".join(suffix) + "\n")
				out.write(str(i) + "\t" + "".join(fragment) + "\n")

		if outLong != None: out.close()
		print str(n_gaps) + " gaps"
		print str(n_errors) + " errors"
		return haplotypes


parser = argparse.ArgumentParser(description='simulateDiploidMatrix.py --matrixout --hapout [--height --width --gaprate --errorrate --fragmentlength --longfmt]')
parser.add_argument('--matrixout', help='output filename for the fragment matrix',required=True)
parser.add_argument('--hapout', help='output filename for the actual haplotype - depending on error rate, may not be the haplotypes corresponding to the MEC score',required=True)
parser.add_argument('--height',help="fragment matrix height (number of fragments)",required=False,default=2)
parser.add_argument('--width',help="fragment matrix width (number of sites)",required=False,default=1)
parser.add_argument('--gaprate',help="fraction of entries within fragments that are gaps (randomly distributed; can be at ends of fragments)",required=False,default=0.0)
parser.add_argument('--errorrate',help="fraction of matrix entries that have errors (randomly distributed) given that there is not a gap at that point",required=False,default=0.0)
#parser.add_argument('--fractionh0',help="fraction of rows that come from 'hap 0' - NOT USED",required=False,default=0.5)
parser.add_argument('--fragmentlength',help="length of fragments between first non-gap and last non-gap (with --paired True it is the length of one end-segment",required=False,default=None)
parser.add_argument('--insertlength',help="number of gaps between the end-segments of the pair",required=False,default=None)
parser.add_argument('--longfmt',help="Write the fragment matrix in long format?",required=False,default=False)
parser.add_argument('--k',help="k>1 number of haplotypes to simulate",required=False,default=2)
parser.add_argument('--pairedend',help="Are reads paired-end?",required=False,default=False)

args = parser.parse_args()
(height,width) = (int(args.height),int(args.width))
gaprate = float(args.gaprate)
errorrate = float(args.errorrate)
k = int(args.k)
#fractionh0 = float(args.fractionh0) if args.fractionh0 != None else None

if args.pairedend and args.fragmentlength != None and (int(args.fragmentlength) < 1 or 2*int(args.fragmentlength) > width):
	print "Since you specified fragment length in paired mode, it must be at least 1 and at most (width-insertlength)/2 (default is (width-insertlength)/2)"
	sys.exit()
if args.pairedend and args.fragmentlength != None and args.insertlength != None and (int(args.insertlength) < 0 or int(args.insertlength) > width - (int(args.fragmentlength)*2)):
	print "Since you specified insert length in paired mode, it must be at least 0 and at most width-(fragmentlength*2) (default is 0)"
	sys.exit()
if not args.pairedend and args.fragmentlength != None and (int(args.fragmentlength) < 1 or int(args.fragmentlength) > width):
	print "Since you specified fragment length (not in paired mode), it must be at least 1 and at most width (default is width)"
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
if args.longfmt == False:
	outfile = None
	outfile1 = args.matrixout
else:
	outfile = args.matrixout + "-long"
	outfile1 = args.matrixout
if k == 2:
	if args.pairedend: H = diploidPairedend(args.fragmentlength,args.insertlength,height, width, gaprate, errorrate,k,outfile,outfile1)
	else: H = diploid(args.fragmentlength,height, width, gaprate, errorrate,k,outfile,outfile1)
else:
	if args.pairedend: H = kploidPairedend(args.fragmentlength,args.insertlength,height, width, gaprate, errorrate,k,outfile,outfile1)
	else: H = kploid(args.fragmentlength,height, width, gaprate, errorrate,k,outfile,outfile1)

with open(args.hapout,"w") as F:
	for h in H:
		F.write("".join([str(h[i]) for i in xrange(width)]))
		F.write("\n")