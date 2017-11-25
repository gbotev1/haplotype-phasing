package src.main;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import com.google.common.collect.TreeMultimap;

/**
 * This class implements the Fast Twister heuristic and should be interacted with via the solve() method.
 * 
 * @author Georgie Botev
 */
public class FastTwister {

	private TreeMultimap<Integer, String> fragments;
	private int fragmentLength;
	private int numFragments;
	
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
		// Take the largest fragment's starting index and augment it by the size of that fragment
		this.fragmentLength = fragments.asMap().lastKey() + Collections.max(fragments.asMap().lastEntry().getValue(), Comparator.comparing(s -> s.length())).length();
		this.fragments = fragments;
	}
	
	/**
	 * This method determines the starting bipartition for the Fast Twister heuristic and then runs the heuristic.
	 * 
	 * It also determines the initial bipartition to feed into the twister method.
	 * 
	 * @return The best-guess haplotype for the given set of fragment reads.
	 */
	public Haplotype solve() {
		// Initialize the starting bipartition
		TreeMultimap<Integer, String> fragments1 = TreeMultimap.create();
		TreeMultimap<Integer, String> fragments2 = TreeMultimap.create();
		// Populate the bipartition with the rest of the fragment data
		for (Map.Entry<Integer, String> entry : this.fragments.entries()) {
			String consensus1 = consensusHaplotype(fragments1).getHaplotype();
			// The current key should always be greater than or equal to the current key of the current consensus haplotype
			if (entry.getValue().charAt(0) == consensus1.charAt(entry.getKey())) {
				fragments1.put(entry.getKey(), entry.getValue());
			} else {
				fragments2.put(entry.getKey(), entry.getValue());
			}
		}
		return twister(fragments1, fragments2);
	}
	
	/**
	 * This method determines the best-guess haplotype using the Fast Twister heuristic given a starting bipartition.
	 * 
	 * @param fragments1 The first set of fragments in the starting bipartition.
	 * @param fragments2 The second set of fragments in the starting bipartition.
	 * @return The best-guess haplotype for the given starting bipartition of fragment reads.
	 */
	private Haplotype twister(TreeMultimap<Integer, String> fragments1, TreeMultimap<Integer, String> fragments2) {
		// Initialize placeholder for best haplotype
		Haplotype bestGuess = null;
		// Trivial upper-bound for MEC score
		int minMEC = this.fragmentLength * this.numFragments;
		// Keep track of MEC scores seen so far
		Map<Integer, Integer> scores = new HashMap<Integer, Integer>();
		// Keep switching until score criterion is met
		while (true) {
			// Calculate current consensus haplotypes
			Bundle b1 = consensusHaplotype(fragments1);
			Bundle b2 = consensusHaplotype(fragments2);
			// Remember current guess for efficiency
			Haplotype currGuess = new Haplotype(b1.getHaplotype(), b2.getHaplotype(), b1.getMEC() + b2.getMEC());
			//System.out.println(currGuess.MEC() + ",");
			// Check if current consensus is minimal
			if (currGuess.MEC() < minMEC) {
				bestGuess = currGuess;
				minMEC = bestGuess.MEC();
			}
			// Record current MEC score
			if (scores.containsKey(currGuess.MEC())) {
				int currScore = scores.get(currGuess.MEC());
				// Stopping criterion
				if (currScore > 0) {
					break;
				} else {
					// Increment score
					scores.put(currGuess.MEC(), currScore + 1);
				}
			} else {
				// Does not contain key
				scores.put(currGuess.MEC(), 1);
			}
			// Extract data
			BadEntryBundle bundle = MEC(b1.getHaplotype(), b2.getHaplotype(), fragments1, fragments2);
			Map.Entry<Integer, String> badEntry1 = bundle.getBadEntry1();
			Map.Entry<Integer, String> badEntry2 = bundle.getBadEntry2();
			MoveStatBundle stats = bundle.getMoveStats();
			int changeRegular = stats.getRegularChange();
			int changeComplement = stats.getComplementChange();
			int changeBoth = stats.getBothChange();
			// Determine the minimum change
			int minChange = Math.min(changeRegular,
					Math.min(changeComplement, changeBoth));
			// Calculate next guess
			if (minChange == changeBoth) {
				if (badEntry1 != null && badEntry2 != null) {
					fragments1.put(badEntry2.getKey(), badEntry2.getValue());
					fragments2.put(badEntry1.getKey(), badEntry1.getValue());
					fragments1.remove(badEntry1.getKey(), badEntry1.getValue());
					fragments2.remove(badEntry2.getKey(), badEntry2.getValue());
				}
			} else if (minChange == changeRegular) {
				if (badEntry1 != null) {
					fragments2.put(badEntry1.getKey(), badEntry1.getValue());
					fragments1.remove(badEntry1.getKey(), badEntry1.getValue());
				}
			} else {
				// minChange == changeComplement
				if (badEntry2 != null) {
					fragments1.put(badEntry2.getKey(), badEntry2.getValue());
					fragments2.remove(badEntry2.getKey(), badEntry2.getValue());
				}
			}
		}
		return bestGuess;
	}
	
	/**
	 * This method determines the MEC score for the given bipartition and respective consensus haplotypes as well as some other relevant statistics.
	 * 
	 * @param haplotype1 The consensus haplotype corresponding to fragments1.
	 * @param haplotype2 The consensus haplotype corresponding to fragments2.
	 * @param fragments1 The first set of fragments in the given bipartition.
	 * @param fragments2 The second set of fragments in the given bipartition.
	 * @return A MECBundle containing the MEC score and relevant statistics or null if the given bipartition has zero MEC score.
	 */
	private BadEntryBundle MEC(String haplotype1, String haplotype2, TreeMultimap<Integer, String> fragments1, TreeMultimap<Integer, String> fragments2) {
		// Initialize counters and placeholders
		int maxScore1 = 0;
		Map.Entry<Integer, String> badEntry1 = null;
		int maxScore2 = 0;
		Map.Entry<Integer, String> badEntry2 = null;
		// Extract iterators for each set in given bipartition
		Iterator<Map.Entry<Integer, String>> fragments1Iter = fragments1.entries().iterator();
		Iterator<Map.Entry<Integer, String>> fragments2Iter = fragments2.entries().iterator();
		// Determine MEC scores and bad indices for fragments1
		for (int i = 0; i < fragments1.size(); i++) {
			int score = 0;
			Map.Entry<Integer, String> currEntry = fragments1Iter.next();
			int currKey = currEntry.getKey();
			String currFragment = currEntry.getValue();
			for (int j = 0; j < currFragment.length(); j++) {
				// Find all places where a site in the current fragment differs from haplotype1
				if (currFragment.charAt(j) != '-' && haplotype1.charAt(currKey +j) != currFragment.charAt(j)) {
					score++;
				}
			}
			// Check if score beats current max
			if (score > maxScore1) {
				maxScore1 = score;
				badEntry1 = currEntry;
			}
		}
		// Next, determine MEC scores and bad indices for fragments2
		for (int i = 0; i < fragments2.size(); i++) {
			int score = 0;
			Map.Entry<Integer, String> currEntry = fragments2Iter.next();
			int currKey = currEntry.getKey();
			String currFragment = currEntry.getValue();
			for (int j = 0; j < currFragment.length(); j++) {
				// Find all places where a site in the current fragment differs from haplotype2
				if (currFragment.charAt(j) != '-' && haplotype2.charAt(currKey + j) != currFragment.charAt(j)) {
					score++;
				}
			}
			// Check if score beats current max
			if (score > maxScore2) {
				maxScore2 = score;
				badEntry2 = currEntry;
			}
		}
		// Designate first entry as bad entry if no real bad entry exists
		if (badEntry1 == null && !fragments1.isEmpty()) {
			badEntry1 = fragments1.entries().iterator().next();
		}
		if (badEntry2 == null && !fragments2.isEmpty()) {
			badEntry2 = fragments2.entries().iterator().next();
		}
		// Finally, determine the approximate change in MEC for switching each of the bad indices
		MoveStatBundle stats = calculateMoveStats(fragments1, fragments2, badEntry1, badEntry2, maxScore1, maxScore2);
		return new BadEntryBundle(badEntry1, badEntry2, stats);
	}

	/**
	 * This method calculates relevant statistics about moving particular fragments within the bipartition.
	 * The bad entries can either be swapped, or either one can be moved without moving the other. All three
	 * of these transformations are considered and evaluated.
	 * 
	 * @param fragments1 The first set of fragments in the given bipartition.
	 * @param fragments2 The second set of fragments in the given bipartition.
	 * @param badEntry1 The bad entry to consider moving from fragments1.
	 * @param badEntry2 The bad entry to consider moving from fragments2.
	 * @param maxScore1 The MEC score associated with badEntry1.
	 * @param maxScore2 The MEC score associated with badEntry2.
	 * @return A MoveStatBundle containing the change in MEC score for each of the three possible transformations.
	 */
	private MoveStatBundle calculateMoveStats(TreeMultimap<Integer, String> fragments1, TreeMultimap<Integer, String> fragments2, Map.Entry<Integer, String> badEntry1, Map.Entry<Integer, String> badEntry2, int maxScore1, int maxScore2) {
		// Initialize counters
		int score1Move = 0;
		int score2Move = 0;
		// Switch bad fragments
		String fragment1 = (badEntry1 != null) ? badEntry1.getValue() : "";
		String fragment2 = (badEntry2 != null) ? badEntry2.getValue() : "";
		if (badEntry2 != null) {
			fragments1.put(badEntry2.getKey(), fragment2);
		}
		if (badEntry1 != null) {
			fragments2.put(badEntry1.getKey(), fragment1);
		}
		// Determine new consensus haplotypes
		String haplotype1New = consensusHaplotype(fragments1).getHaplotype();
		String haplotype2New = consensusHaplotype(fragments2).getHaplotype();
		// Calculate scores after switches
		for (int i = 0; i < fragment1.length(); i++) {
			if ((badEntry1 != null) && fragment1.charAt(i) != '-' && haplotype2New.charAt(badEntry1.getKey() + i) != fragment1.charAt(i)) {
				score1Move++;
			}
		}
		for (int i = 0; i < fragment2.length(); i++) {
			if ((badEntry2 != null) && fragment2.charAt(i) != '-' && haplotype1New.charAt(badEntry2.getKey() + i) != fragment2.charAt(i)) {
				score2Move++;
			}
		}
		// Undo switching of fragments
		if (badEntry2 != null) {
			fragments1.remove(badEntry2.getKey(), fragment2);
		}
		if (badEntry1 != null) {
			fragments2.remove(badEntry1.getKey(), fragment1);
		}
		// Calculate relevant statistics
		int changeRegular = score1Move - maxScore1;
		int changeComplement = score2Move - maxScore2;
		int changeBoth = score1Move + score2Move - maxScore1 - maxScore2;
		return new MoveStatBundle(changeRegular, changeComplement, changeBoth);
	}
	
	/**
	 * This method calculates the consensus haplotype and the associated MEC score.
	 * 
	 * @param fragments The set of fragments to process.
	 * @return A bundle containing the MEC score and the consensus haplotype for the given fragments.
	 */
	private Bundle consensusHaplotype(TreeMultimap<Integer, String> fragments) {
		// Check if fragments is empty
		if (fragments.isEmpty()) {
			StringBuilder sb = new StringBuilder(this.fragmentLength);
			for (int i = 0; i < this.fragmentLength; i++) {
				sb.append('0');
			}
			return new Bundle(0, sb.toString());
		}
		// Initialize frequency arrays
		int[] num0 = new int[this.fragmentLength];
		int[] num1 = new int[this.fragmentLength];
		// Determine frequencies
		for (Map.Entry<Integer, String> entry : fragments.entries()) {
			int currKey = entry.getKey();
			String currFragment = entry.getValue();
			int currFragmentLen = currFragment.length();
			for (int i = 0; i < currFragmentLen; i++) {
				if (currFragment.charAt(i) == '0') {
					num0[currKey + i] += 1;
				} else if (currFragment.charAt(i) == '1') {
					num1[currKey + i] += 1;
				}
				// Else, currFragment.charAt(i) == '-', so do nothing
			}
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
	
}
