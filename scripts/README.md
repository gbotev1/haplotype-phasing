###simulateKploidMatrix: Utility for simulating fragment matrixes from diploid up to k-ploid to test phasing algorithms.


usage: simulateKploidMatrix.py [-h] --matrixout MATRIXOUT --hapout HAPOUT
                               [--height HEIGHT] [--width WIDTH]
                               [--gaprate GAPRATE] [--errorrate ERRORRATE]
                               [--fractionh0 FRACTIONH0]
                               [--fragmentlength FRAGMENTLENGTH] [--k K]

optional arguments:
  -h, --help            show this help message and exit  
  --matrixout MATRIXOUT
                        output filename for the fragment matrix  
  --hapout HAPOUT       output filename for the actual haplotype - depending
                        on error rate, may not be the haplotypes corresponding
                        to the MEC score  
  --height HEIGHT       fragment matrix height (number of fragments)  
  --width WIDTH         fragment matrix width (number of sites)  
  --gaprate GAPRATE     fraction of entries within fragments that are gaps
                        (randomly distributed)  
  --errorrate ERRORRATE
                        fraction of matrix entries that have errors (randomly
                        distributed) given that there is not a gap at that
                        point  
  --fractionh0 FRACTIONH0
                        fraction of rows that come from 'hap 0' (only applies
                        to k=2)  
  --fragmentlength FRAGMENTLENGTH
                        length of fragments between first non-gap and last
                        non-gap  
  --k K                 k>1 number of haplotypes to simulate  
  
####Notes  
In the whole matrix, fragment coverage between two sites differs by at most 2 (if no fractionh0 specified in the case of k=2)   
Fraction of fragments that come from each haplotype in k>2-ploid is not modeled  
There are no spaces separating characters in the output files; a newline separates fragments or haplotypes  


####Examples  
`python simulateKploidMatrix.py --matrixout mat.txt --hapout hap.txt --fragmentlength 3 --height 17 --errorrate 0.2 --width 5`
`cat hap.txt`  
00101  
`cat mat.txt`  
001--  
001--  
-110-  
-011-  
--111  
--101  
001--  
-010-  
100--  
100--  
100--  
-110-  
-100-  
-111-  
--000  
--110  
--110  

`python simulateKploidMatrix.py --matrixout mat.txt --hapout hap.txt --height 11 --gaprate 0.1 --width 5`  
`cat hap.txt`  
00100  
`cat mat.txt`  
00-00  
00100  
001-0  
001-0  
00100  
11011  
11011  
11011  
1-011  
11011  
11011   

`python simulateKploidMatrix.py --matrixout mat.txt --hapout hap.txt --height 11 --gaprate 0.1 --width 5 --k 3`  
`cat hap.txt`  
01100  
10010  
10001  
`cat mat.txt`  
01100  
-1100  
0110-  
1-010  
10010  
100-0  
10001  
10001  
10001  
-1100  
10010  