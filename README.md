# Haplotype Phasing
A workspace for my current undergraduate research at Schatz Lab at Johns Hopkins University.

# Dependencies
This project uses Google Guava 23.0, which is included in the repository so that you do not have to download it separately. The JAR used for this project was downloaded from <https://github.com/google/guava/wiki/Release23>.

# Running
In order to run the parallel solver, follow the following steps:

1. Run either the simulateFragmentMatrix.py script with the "--shortfmt true" option to enable short fragments or the simulate.py script, and remember the absolute path to the resultant SNP matrix text file.

2. Configure the makefile to match your machine setup.

3. make clean

4. make

5. make run args="[absolute-path-to-mat.txt] [k] [alpha] [beta] [seedLength] [fragmentLength] [prettyPrint] [inclusiveSeeding]"
	a. A recommended argument configuration would be [absolute-path-to-mat.txt] 4 2.0 1 3 10 false false
