
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


switch_total = 0
flip_total = 0

for b in blocks:
	if len(b) < 2: continue
	# Initialize comparison string
	xorHap = ""
	switch = 0
	flip = 0

	# Process comparison string
	i = 0
	while i < len(b):
		if b[i] == 1:
			flag = False
			while i < len(b) - 1 and b[i] == b[i + 1]:
				i += 1
				flag = True
			if flag:
				switch += 1
			else:
				flip += 1
		i += 1
	print flip,switch,len(b)
	switch_total += switch
	flip_total += flip

print "Switch:", str(switch_total)	
print "Flip:", str(flip_total)
