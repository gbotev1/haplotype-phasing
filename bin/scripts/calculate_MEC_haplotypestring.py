#cat haplotype-phasing/data/chr18_hap.txt haplotype-phasing/data/chr22_PLAT_shortformat_matrix.txt  | python calculate_MEC_haplotypestring.py
#cat [haplotype as single line of 0's and 1's] [fragment matrix file] | python [script name]

import sys

H = list(sys.stdin.readline().strip())

MEC = 0
entries_total = 0
for line in sys.stdin:
	L = line.strip().split()
	offset = int(L[0])-1
	fragment = list(L[1])
	entries = len(fragment) - fragment.count("-")
	entries_total += entries

	mismatch = 0
	for (i,f) in enumerate(fragment):
		if f == "-": continue
		if f != H[offset+i]: mismatch += 1
	m = min(mismatch,entries-mismatch) 
	#accounts for possibility that fragment matches complement of haplotype
	MEC += m

print "MEC:", str(MEC)
print "Fragment Matrix 0/1 Entries:", str(entries_total)