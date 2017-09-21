package src.main;

import java.util.Set;
import java.util.HashSet;

import com.google.common.collect.Sets;

/**
 * This class solves the haplotype reconstruction problem for diploid chromosomes using MEC via an exponential algorithm.
 * 
 * @author gbotev
 */
public class MECExponentialSolver {
	
	private char[][] SNPMatrix;
	
	public MECExponentialSolver(char[][] SNPMatrix) {
		this.SNPMatrix = SNPMatrix;
	}
	
	/**
	 * This method finds the best-guess haplotypes based on MEC.
	 * 
	 * Suggestions for improvement: Do not hold entire power set in memory; Iterate one at a time
	 * @return 
	 */
	public Haplotype solve() {
		Set<String> fragments = generateSetOfFragments();
		Set<Set<String>> powerSetFragments = Sets.powerSet(fragments);
		// Keep track of the best fragment
		Haplotype bestGuess = null;
		int minMEC = SNPMatrix[0].length * SNPMatrix.length;
		for (Set<String> fragment : powerSetFragments) {
			// Fragment and fragmentComplement are the haplotype split we are considering
			Set<String> fragmentComplement = Sets.difference(fragments, fragment);
			// Guess the best haplotypes and compute MEC score
			Haplotype currGuess = calculateHaplotype(fragment, fragmentComplement);
			// Check if current guess beats any seen so far
			if (currGuess.MEC() < minMEC) {
				bestGuess = currGuess;
				minMEC = bestGuess.MEC();
				bestGuess.setPartition(fragmentComplement);
			}
		}
		return bestGuess;
	}
	
	/**
	 * This method generates the best guess haplotypes and corresponding MEC score for a given partition of the fragments.
	 * @param s1 The first set of fragments.
	 * @param s2 The second set of fragments.
	 * @return The best guess haplotypes and corresponding MEC score.
	 */
	private Haplotype calculateHaplotype(Set<String> s1, Set<String> s2) {
		String h1 = bestGuessHaplotype(s1);
		String h2 = bestGuessHaplotype(s2);
		int MEC1 = MEC(h1, s1);
		int MEC2 = MEC(h2, s2);
		// Use combined MEC score for the MEC score of both haplotypes
		return new Haplotype(h1, h2, MEC1 + MEC2);
	}
	
	/**
	 * This method determines the MEC score for a given haplotype guess and set of fragments.
	 * @param haplotype The best guess for the haplotype from the given fragments.
	 * @param fragments The set of fragments used to generate the best guess.
	 * @return The MEC score for the given configuration.
	 */
	private int MEC(String haplotype, Set<String> fragments) {
		int MEC = 0;
		for (String fragment : fragments) {
			for (int i = 0; i < fragment.length(); i++) {
				// Find all places where a place in the fragment differs from the haplotype
				if (fragment.charAt(i) != '-' && haplotype.charAt(i) != fragment.charAt(i)) {
					MEC++;
				}
			}
		}
		return MEC;
	}
	
	/**
	 * This method determines the best guess for the haplotype based on the set of fragments given.
	 * @param fragments The set of fragments to consider.
	 * @return The best guess for the haplotype.
	 */
	private String bestGuessHaplotype(Set<String> fragments) {
		// Note that all strings have the same length because they come from SNP matrix
		SNPFrequency[] frequencies = new SNPFrequency[SNPMatrix[0].length];
		// Initialize all frequencies
		for (int i = 0; i < frequencies.length; i++) {
			frequencies[i] = new SNPFrequency();
		}
		// Determine frequencies
		for (String fragment : fragments) {
			for (int i = 0; i < fragment.length(); i++) {
				switch (fragment.charAt(i)) {
					case '1':
						frequencies[i].increment1();
						break;
					case '0':
						frequencies[i].increment0();
					default:
						// All cases covered
						break;
				}
			}
		}
		// Construct best guess
		StringBuilder bestGuess = new StringBuilder();
		for (SNPFrequency f : frequencies) {
			bestGuess.append(f.maxFrequency());
		}
		return bestGuess.toString();
	}
	
	/**
	 * This method extracts all of the fragments from the SNP matrix as a Set.
	 * @return The set of all fragments in the SNP matrix.
	 */
	private Set<String> generateSetOfFragments() {
		Set<String> fragments = new HashSet<String>();
		for (int row = 0; row < SNPMatrix.length; row++) {
			StringBuilder fragment = new StringBuilder();
			for (int column = 0; column < SNPMatrix[0].length; column++) {
				fragment.append(SNPMatrix[row][column]);
			}
			fragments.add(fragment.toString());
		}
		return fragments;
	}
	
}
