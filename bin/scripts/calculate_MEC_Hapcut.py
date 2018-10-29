#cat haplotype-phasing/data/chr22_PLAT_shortformat_matrix.txt | python calculate_MEC_Hapcut.py chr22.hap
#cat [fragment matrix file] | python [script name] [HapCut2 output file]

import sys

hap_file = open(sys.argv[1])

blocks = []
for line in hap_file.readlines():
	L = line.strip().split()
	if "*" in L[0]:
		continue
	if "BLOCK" in L[0]:
		blocks.append([])
		for i in xrange(int(L[2])+int(L[4])-1):
			blocks[-1].append("-")
		continue
	blocks[-1][int(L[0])-1] = L[1]

MEC = 0
entries_total = 0
for line in sys.stdin:
	L = line.strip().split()
	offset = int(L[0])-1
	fragment = list(L[1])
	entries = len(fragment) - fragment.count("-")
	entries_total += entries
	for b in blocks:
		mismatch = 0
		if len(b) < offset+len(fragment): #fragment not contained in block span
			continue
		for (i,f) in enumerate(fragment):
			if f == "-": continue
			if b[offset+i] == "-": break #fragment not in block
			if f != b[offset+i]: mismatch += 1
		#if you made it here, this is the right block
		m = min(mismatch,entries-mismatch) 
		#accounts for possibility that fragment matches complement of haplotype
		MEC += m
		break

print "MEC:", str(MEC)
print "Fragment Matrix 0/1 Entries:", str(entries_total)