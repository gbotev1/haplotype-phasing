package src.main;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.TreeMultimap;

/**
 * This class implements the Fast Twister heuristic and should be interacted with via the solve() method.
 * 
 * @author gbotev
 */
public class FastTwister {

	private TreeMultimap<Integer, String> fragments;
	private int fragmentLength;
	private int numFragments;
	private Partition partition1;
	private Partition partition2;
	
	/**
	 * Accepts a tree multimap of fragments that are encoded in the compact notation.
	 * The key represents the starting index of the fragment, and the value represents
	 * the fragment itself.
	 * 
	 * @param fragments The set of fragments to process.
	 */
	public FastTwister(TreeMultimap<Integer, String> fragments) {
		// Number of keys
		this.numFragments = fragments.size();
		System.out.println(this.numFragments);
		this.fragments = fragments;
		int counter = 0;
		int max = 0;
		for (Map.Entry<Integer, String> entry : this.fragments.entries()) {
			int currLength = entry.getKey() + entry.getValue().length();
			if (currLength > max) {
				max = currLength;
			}
			for (int i = 0; i < entry.getValue().length(); i++) {
				if (entry.getValue().charAt(i) != '-') {
					counter++;
				}
			}
		}
		this.fragmentLength = max;
		System.out.println(counter);
	}
	
	/**
	 * This method determines the starting bipartition for the Fast Twister heuristic and then runs the heuristic.
	 * 
	 * It also determines the initial bipartition to feed into the twister method.
	 * 
	 * @return A bundle containing the best-guess haplotype for the given set of fragment reads.
	 */
	public Bundle solve() {
		System.out.println("Beginning initial split:");
		// Initialize the starting bipartition
		TreeMultimap<Integer, String> fragments1 = TreeMultimap.create();
		TreeMultimap<Integer, String> fragments2 = TreeMultimap.create();
		// Initialize the consensus haplotypes for both partitions
		StringBuilder sb1 = new StringBuilder(this.fragmentLength);
		StringBuilder sb2 = new StringBuilder(this.fragmentLength);
		for (int i = 0; i < this.fragmentLength; i++) {
			sb1.append('0');
			sb2.append('1');
		}
		String consensus1 = sb1.toString();
		String consensus2 = sb2.toString();
		// Initialize frequency arrays to avoid unnecessary computations
		int[] num01 = new int[this.fragmentLength];
		int[] num11 = new int[this.fragmentLength];
		int[] num02 = new int[this.fragmentLength];
		int[] num12 = new int[this.fragmentLength];
		int counter = 0;
		// Populate the bipartition with the rest of the fragment data
		for (Map.Entry<Integer, String> entry : this.fragments.entries()) {
			counter++;
			// Initialize counters
			int sum1 = 0;
			int sum2 = 0;
			// Save constants for efficiency
			int entryKey = entry.getKey();
			String entryValue = entry.getValue();
			int entryValueLength = entryValue.length();
			// Calculate distance score
			for (int i = 0; i < entryValueLength; i++) {
				// Save current character for efficiency
				char currChar = entryValue.charAt(i);
				// CONSENSUS FIRST PARTITION
				if (currChar == consensus1.charAt(entryKey + i)) {
					// Note that consensus never has '-'
					sum1 += 1;
				} else if (currChar != '-') {
					sum1 -= 1;
				} // Else, currChar == '-' AND currChar != consensus1.charAt(entryKey + i), so do nothing
				// CONSENSUS SECOND PARTITION
				if (currChar == consensus2.charAt(entryKey + i)) {
					// Note that consensus never has '-'
					sum2 += 1;
				} else if (currChar != '-') {
					sum2 -= 1;
				} // Else, currChar == '-' AND currChar != consensus2.charAt(entryKey + i), so do nothing
			}
			// Place fragment according to formula
			if (sum1 >= sum2) {
				fragments1.put(entryKey, entryValue);
				// Update the consensus haplotype incrementally
				consensus1 = consensusHaplotypeIncremental(num01, num11, entry).getHaplotype();
			} else {
				fragments2.put(entryKey, entryValue);
				// Update the consensus haplotype incrementally
				consensus2 = consensusHaplotypeIncremental(num02, num12, entry).getHaplotype();
			}
			if (counter % 5000 == 0) {
				System.out.println(counter);
			}
		}
		// Initialize partition objects
		this.partition1 = new Partition(fragments1, num01, num11, this.fragmentLength);
		this.partition2 = new Partition(fragments2, num02, num12, this.fragmentLength);
		// Perform the twisting
		return twister();
	}
	
	/**
	 * This method determines the best-guess haplotype using the Fast Twister heuristic given a starting bipartition.
	 * 
	 * @param fragments1 The first set of fragments in the starting bipartition.
	 * @param fragments2 The second set of fragments in the starting bipartition.
	 * @return A bundle containing the best-guess haplotype for the given starting bipartition of fragment reads.
	 */
	private Bundle twister() {
		System.out.println("Beginning twisting:");
		// Initialize placeholder for best haplotype
		Bundle bestGuess = null;
		// Trivial upper-bound for MEC score
		int minMEC = Integer.MAX_VALUE;
		// Keep track of MEC scores seen so far
		Map<Integer, Integer> scores = new HashMap<Integer, Integer>();
		// Keep switching until score criterion is met
		while (true) {
			// Remember current guess for efficiency
			Bundle currGuess = this.consensusHaplotype(this.partition1, this.partition2);
			System.out.println(currGuess.getMEC() + ",");
			// Check if current consensus is minimal
			if (currGuess.getMEC() < minMEC) {
				bestGuess = currGuess;
				minMEC = bestGuess.getMEC();
			}
			// Record current MEC score
			if (scores.containsKey(currGuess.getMEC())) {
				int currScore = scores.get(currGuess.getMEC());
				// Stopping criterion
				if (currScore > 5) {
					break;
				} else {
					// Increment score
					scores.put(currGuess.getMEC(), currScore + 1);
				}
			} else {
				// Does not contain key
				scores.put(currGuess.getMEC(), 1);
			}
			// Calculate switch fragments
			Set<Map.Entry<Integer, String>> switchFragments1 = this.partition1.getSwitchFragments();
			Set<Map.Entry<Integer, String>> switchFragments2 = this.partition2.getSwitchFragments();
			// Switch fragments
			this.partition1.addFragments(switchFragments2);
			this.partition2.addFragments(switchFragments1);
			this.partition1.removeFragments(switchFragments1);
			this.partition2.removeFragments(switchFragments2);
		}
		return bestGuess;
	}
	
	private void print(int[] array) {
		System.out.print("[");
		for (int i = 0; i < 100; i++) {
			System.out.print(array[i] + ", ");
		}
		System.out.print("]\n");
	}
	
	/**
	 * This method allows for the calculation of the consensus haplotype and
	 * the associated MEC score in an incremental manner for efficiency by using
	 * references to the frequency arrays.
	 * 
	 * @param num0 The frequency array for the number of zeros.
	 * @param num1 The frequency array for the number of ones.
	 * @param fragment The fragment being added to the partition under consideration.
	 * @return A bundle containing the MEC score and the consensus haplotype after
	 * the current fragment has been added to the partition under consideration.
	 */
	private Bundle consensusHaplotypeIncremental(int[] num0, int[] num1, Map.Entry<Integer, String> fragment) {
		// Update frequencies
		int currKey = fragment.getKey();
		String currFragment = fragment.getValue();
		int currFragmentLen = currFragment.length();
		for (int i = 0; i < currFragmentLen; i++) {
			if (currFragment.charAt(i) == '0') {
				num0[currKey + i] += 1;
			} else if (currFragment.charAt(i) == '1') {
				num1[currKey + i] += 1;
			}
			// Else, currFragment.charAt(i) == '-', so do nothing
		}
		// Initialize MEC score counter
		int mec = 0;
		// Initialize StringBuilder to store partial consensus haplotype as we compute it for efficiency
		StringBuilder consensus = new StringBuilder(this.fragmentLength);
		// Determine consensus haplotype and MEC score
		for (int i = 0; i < this.fragmentLength; i++) {
			if (num0[i] >= num1[i]) {
				consensus.append('0');
				// Only ones in this column contribute to MEC score
				mec += num1[i];
			} else {
				// num0[i] < num1[i]
				consensus.append('1');
				// Only zeros in this column contribute to MEC score
				mec += num0[i];
			}
		}
		return new Bundle(mec, consensus.toString());
	}
	
	private Bundle consensusHaplotype(Partition p1, Partition p2) {
		// Initialize frequency arrays for bipartition
		int[] num01 = p1.get0Frequency();
		int[] num11 = p1.get1Frequency();
		int[] num02 = p2.get0Frequency();
		int[] num12 = p2.get1Frequency();
		// Initialize MEC score counter
		int mec = 0;
		// Initialize StringBuilder to store partial consensus haplotype as we compute it for efficiency
		StringBuilder consensus = new StringBuilder(this.fragmentLength);
		// Determine consensus haplotype and MEC score
		for (int i = 0; i < this.fragmentLength; i++) {
			if (num01[i] > num11[i]) {
				consensus.append('0');
				mec += num11[i];
				mec += num02[i];
			} else if (num01[i] < num11[i]) {
				consensus.append('1');
				mec += num01[i];
				mec += num12[i];
			} else {
				// num01[i] == num11[i]
				// Break ties so that haplotypes are complements of each other
				if (num02[i] <= num12[i]) {
					// Break complete tie with zero
					consensus.append('0');
					mec += num11[i];
					mec += num02[i];
				} else {
					// num02[i] > num12[i]
					consensus.append('1');
					mec += num01[i];
					mec += num12[i];
				}
			}
		}
		return new Bundle(mec, consensus.toString());
	}
	
}
