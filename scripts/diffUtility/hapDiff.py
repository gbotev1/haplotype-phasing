# Read in the two haplotypes to compare
with open("h1.txt", "r") as infile:
	real = infile.readlines()[0]
with open("h2.txt", "r") as infile:
	found = infile.readlines()[0]

# Initialize counters for different kinds of errors
switch = 0
flip = 0

# Save length of smaller haplotype (if they are different lengths?)
size = min(len(real), len(found))

# Initialize comparison string
xorHap = ""

# Create comparison string
for i in xrange(size):
	xorHap += str(int(real[i]) ^ int(found[i]))

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
	
# Display results	
print("Switch: %d" % switch)
print("Flip: %d" % flip)