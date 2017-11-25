# Haplotype Phasing
A workspace for my current undergraduate research position at Schatz Lab at Johns Hopkins University.

# Dependencies
This project uses Google Guava 23.0, which is included in the repository so that you do not have to download it separately. The JAR used for this project was downloaded from <https://github.com/google/guava/wiki/Release23>.

# Limitations
The exponential solver cannot be run on more than 30 fragments because the Google Guava powerSet method is internally limited at this amount. The Fast Twister heuristic has no such limit.

# Running
In order to run the solvers, follow the following steps:

1. Run the simulateFragmentMatrix.py script with the "--shortfmt true" option to enable short fragments in order to run the Fast Twister heuristic. Otherwise, use the simulateKploidMatrix.py script to run the exponential solver.

2. Configure the tester.java file appropriately (based on step #1). Make sure to configure the path to the matrix of fragments according to your machine setup. It is set up to run the Fast Twister heuristic.

3. Configure the makefile to match your machine setup.

4. "make clean"

5. "make"

6. "make run"
