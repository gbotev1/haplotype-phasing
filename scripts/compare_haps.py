
#cat haplotype-phasing/data/chr18_hap.txt chr18.hap | python compare_haps.py
#cat [haplotype as single line of 0's and 1's] [HapCut2 output file] | python [script name]
import sys
H = list(sys.stdin.readline().strip())

blocks = [[]]
for line in sys.stdin:
	L = line.strip().split()
	if "*" in L[0] or "BLOCK" in L[0]:
		blocks.append([])
		continue
	blocks[-1].append(H[int(L[0])-1] == L[1])


mismatch_total = 0
switch_total = 0
for b in blocks:
	mismatch = 0
	switch = 0
	if len(b) < 2:
		continue
	if b[0] != b[1]:
		mismatch += 1
	for i in xrange(2,len(b)):
		if b[i] != b[i-1]:
			if b[i] == b[i-2]:
				mismatch += 1
			else: #b[i] != b[i-2]
				switch += 1
	print len(b),mismatch,switch
	mismatch_total += mismatch
	switch_total += switch	

print "Mismatch:", str(mismatch_total)
print "Switch:", str(switch_total)