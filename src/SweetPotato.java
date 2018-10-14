import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeMap;

import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeMultiset;

/**
 * This class implements a phasing algorithm inspired by the sweet potato paper
 * for k-ploid genomes.
 * 
 * @author gbotev
 */
public class SweetPotato {

	private int tag = 0;
	private Set<Fragment> fragments;
	private int numFragments;
	private int numSNP;
	private int k;
	private double alpha;
	private int beta;
	private int seedLength;
	private int fragmentLength;
	private boolean prettyPrint;
	private List<FrequencyArray> seedHaplotypes = new ArrayList<FrequencyArray>();
	private PriorityQueue<FrequencyArrayPair> faPairs = new PriorityQueue<FrequencyArrayPair>(new Comparator<FrequencyArrayPair>() {
		@Override
		public int compare(FrequencyArrayPair fap1, FrequencyArrayPair fap2) {
			return -Double.compare(fap1.getScore(), fap2.getScore());
		}
	});

	public SweetPotato(Set<Fragment> fragments, int k, double alpha, int beta, int seedLength, int fragmentLength, boolean prettyPrint) {
		this.fragments = fragments;
		// Number of keys
		this.numFragments = fragments.size();
		for (Fragment f : this.fragments) {
			if (f.endIndex() > this.numSNP) {
				this.numSNP = f.endIndex();
			}
		}
		// Starts at zero and is inclusive of last SNP site
		this.numSNP += 1;
		// Save k-ploid specification
		this.k = k;
		// Save hyperparameters
		this.alpha = alpha;
		this.beta = beta;
		// Save training parameters
		this.seedLength = seedLength;
		this.fragmentLength = fragmentLength;
		// Save formatting parameters
		this.prettyPrint = prettyPrint;
	}
	
	/**
	 * This method finds the best k seeds that are believed to come from different 
	 * haplotypes in terms of the average number of supporting fragments.
	 * @param numSeedHaplotypes The number of seed haplotypes should be passed for efficiency.
	 * @return The starting index of this.seedHaplotypes for which the next k seeds
	 * have the largest average number of supporting fragments.
	 */
	private int bestRangeHelper(int numSeedHaplotypes) {
		double bestAverage = 0.0;
		int bestStart = 0;
		for (int start = 0; start < numSeedHaplotypes - this.k; start += this.k) {
			int sum = 0;
			for (int offset = 0; offset < this.k; offset++) {
				sum += this.seedHaplotypes.get(start + offset).numSupportingFrags();
			}
			double currAverage = ((double) (sum)) / ((double) this.k);
			if (currAverage > bestAverage) {
				bestAverage = currAverage;
				bestStart = start;
			}
		}
		return bestStart;
	}
	
	/**
	 * This method initializes the best initial merges with which to initialize each Merger thread.
	 * @return An array of FrequencyArrayPairs of length exactly k; each element corresponds to
	 * an initial merge that is believed to be the best for that particular haplotype.
	 */
	private FrequencyArrayPair[] initializeBestMerges() {
		// Initialize the return array
		FrequencyArrayPair[] results = new FrequencyArrayPair[this.k];
		// Save current number of seeds for efficiency
		int numSeedHaplotypes = this.seedHaplotypes.size();
		// Grab best k seeds that are believed to come from different haplotypes
		FrequencyArray[] bestKSeeds = new FrequencyArray[this.k];
		int start = this.bestRangeHelper(numSeedHaplotypes);
		for (int i = 0; i < this.k; i++) {
			bestKSeeds[i] = this.seedHaplotypes.get(start + i);
		}
		for (int i = 0; i < bestKSeeds.length; i++) {
			FrequencyArray fa1 = bestKSeeds[i];
			double fa1Size = fa1.numSupportingFrags();
			Set<Fragment> fa1Frags = fa1.getFrags();
			// Calculate best merge with other FrequencyArrays
			for (int j = 0; j < numSeedHaplotypes; j++) {
				if (j >= start && j < start + this.k) {
					continue;
				}
				FrequencyArray fa2 = this.seedHaplotypes.get(j);
				double fa2Size = fa2.numSupportingFrags();
				// We do not have reason to believe that either
				// fa1 or fa2 will be smaller in size, so pass in
				// function arguments in any order
				double intersectionSize = Sets
						.intersection(fa1Frags, fa2.getFrags()).size();
				// Use Jaccard index
				double currIndex = intersectionSize
						/ (fa1Size + fa2Size - intersectionSize);
				if (intersectionSize > 0) {
					this.faPairs.add(new FrequencyArrayPair(fa1, fa2, currIndex));
				}
			}
			// Take the best merge determined by score
			results[i] = this.faPairs.poll();
			// Clear the priority queue for the next iteration
			this.faPairs.clear();
		}
		// Check if results contains null
		boolean hasNull = false;
		for (int i = 0; i < results.length; i++) {
			if (results[i] == null) {
				hasNull = true;
				break;
			}
		}
		// Return appropriately
		if (hasNull) {
			return null;
		} else{
			return results;
		}
	}
	
	public void phaseParallel() {
		// Determine seeds
		System.err.println("Seeding");
		this.seed(new int[this.seedLength], fragmentLength, 0);
		// Find the best initial merges
		System.err.println("Pairing");
		FrequencyArrayPair[] bestMerges = this.initializeBestMerges();
		// Merge best-guesses for seeds in parallel
		System.err.println("Merging");
		Thread[] threads = new Thread[this.k];
		// Start threads
		for (int i = 0; i < this.k; i++) {
			threads[i] = new Thread(new Merger(this.seedHaplotypes, bestMerges[i], this.numSNP, this.prettyPrint));
			threads[i].start();
		}
		// Wait for them to finish
		for (Thread t : threads) {
			try {
				t.join();
			} catch (InterruptedException e) {
				System.err.println("One of the Merger threads was interrupted.");
			}
		}
		// Indicate that phasing finished successfully!
		System.err.println("Finished phasing");
	}

	private static class ValueComparator implements Comparator<String> {
		private HashMap<String, Integer> map = new HashMap<String, Integer>();

		public ValueComparator(HashMap<String, Integer> map) {
			this.map.putAll(map);
		}

		@Override
		public int compare(String s1, String s2) {
			if (this.map.get(s1) > this.map.get(s2)) {
				return -1;
			} else if (this.map.get(s1) < this.map.get(s2)) {
				return 1;
			} else {
				// this.map.get(s1) == this.map.get(s2))
				return 0;
			}
		}
	}
	
	private boolean isMatch(String seed, int[] indices, Fragment f) {
		for (int i = 0; i < indices.length; i++) {
			if (!(indices[i] >= f.startIndex() && indices[i] <= f.endIndex() && f.toString().charAt(indices[i] - f.startIndex()) == seed.charAt(i))) {
				return false;
			}
		}
		return true;
	}
	
	private boolean isValid(int[] indices, Fragment f) {
		for (int i = 0; i < indices.length; i++) {
			if (!(indices[i] >= f.startIndex() && indices[i] <= f.endIndex() && f.toString().charAt(indices[i] - f.startIndex()) != '-')) {
				return false;
			}
		}
		return true;
	}

	private TreeMap<String, Integer> seedFrequency(int[] indices) {
		HashMap<String, Integer> seedFrequencies = new HashMap<String, Integer>(indices.length);
		for (Fragment f : this.fragments) {
			if (this.isValid(indices, f)) {
				// Seed SNP sites are valid for current fragment, so add if unique
				StringBuilder sb = new StringBuilder(indices.length);
				for (int i = 0; i < indices.length; i++) {
					sb.append(f.toString().charAt(indices[i] - f.startIndex()));
				}
				String seed = sb.toString();
				if (!seedFrequencies.containsKey(seed)) {
					// Does NOT contain seed, so add it with its frequency
					seedFrequencies.put(seed, f.frequency());
				} else {
					// CONTAINS seed, so increment its frequency
					seedFrequencies.put(seed, seedFrequencies.get(seed) + f.frequency());
				}
			}
		}
		TreeMap<String, Integer> sortedSeedFrequencies = new TreeMap<String, Integer>(new ValueComparator(seedFrequencies));
		sortedSeedFrequencies.putAll(seedFrequencies);
		return sortedSeedFrequencies;
	}
	
	private void doSeed(int[] indices) {
		TreeMap<String, Integer> seedFrequency = this.seedFrequency(indices);
		if (seedFrequency.size() >= this.k) {
			Set<FrequencyArray> currSeedGroup = new HashSet<FrequencyArray>(this.k);
			int size = 0;
			for (int i = 0; i < this.k; i++) {
				FrequencyArray fa = new FrequencyArray(this.numSNP, this.tag);
				// For top k (already sorted) find consensus
				Map.Entry<String, Integer> currEntry = seedFrequency.pollFirstEntry();
				String currSeed = currEntry.getKey();
				Set<Fragment> currSupport = this.supportingFragments(currSeed, indices);
				fa.addFragment(currSupport);
				size += currSupport.size();
				currSeedGroup.add(fa);
			}
			// Process currSeedGroup here to remove errors via twisting step
			Set<FrequencyArray> result = currSeedGroup;
			while (true) {
				// Twist until we cannot twist anymore (size does not change)
				Set<FrequencyArray> result2 = this.twist(result, this.tag, size);
				if (result2 == null) {
					// Twisting failed
					return;
				}
				if (result2.equals(result)) {
					// They are same
					result = result2;
					break;
				}
				result = result2;
			}
			this.seedHaplotypes.addAll(result);
			this.tag++;
		}
	}
	
	private Set<FrequencyArray> twist(Set<FrequencyArray> fas, int tag, int size) {
		// Extract all supporting fragments from current seed group 
		Set<Fragment> supportingFrags = new HashSet<Fragment>(size + size / 2);
		for (FrequencyArray fa : fas) {
			supportingFrags.addAll(fa.getFrags());
		}
		
		ArrayList<Fragment> consensuses = new ArrayList<Fragment>(this.k);
		ArrayList<Set<Fragment>> fasNewRaw = new ArrayList<Set<Fragment>>(this.k);
		// Initialization step
		for (int i = 0; i < this.k; i++) {
			fasNewRaw.add(new HashSet<Fragment>());
		}
		// Calculate current consensus for efficiency
		for (FrequencyArray fa : fas) {
			consensuses.add(fa.consensus());
		}
		TreeMultiset<Map.Entry<Integer, Integer>> scores = TreeMultiset.create(Map.Entry.comparingByKey());
		for (Fragment f : supportingFrags) {
			for (int i = 0; i < this.k; i++) {
				int score = consensuses.get(i).similarTo(f);
				scores.add(new AbstractMap.SimpleImmutableEntry<Integer, Integer>(score, i));
			}
			Multiset.Entry<Map.Entry<Integer, Integer>> highScore = scores.lastEntry();
			if (highScore.getCount() == 1 && highScore.getElement().getKey() >= f.length() - this.beta) {
				// DEFAULT: f.length() - 1
				// The membership of this fragment is unambiguous because it matches 1 of the consensus
				// haplotypes and there are either no errors or 1 error (when using magic number 8).
				fasNewRaw.get(highScore.getElement().getValue()).add(f);
			}
			scores.clear();
		}
		// Save for efficiency
		double cutoff = ((double)(supportingFrags.size()) / (this.alpha * this.k));
		for (int i = 0; i < this.k; i++) {
			if (fasNewRaw.get(i).size() <= cutoff) {
				// If one of the seeds does not have "enough" supporting fragments
				// Or... strictly less than k partitions
				return null;
			}
		}
		
		Set<FrequencyArray> fasNew = new HashSet<FrequencyArray>();
		for (int i = 0; i < this.k; i++) {
			FrequencyArray fa = new FrequencyArray(this.numSNP, tag);
			fa.addFragment(fasNewRaw.get(i));
			fasNew.add(fa);
		}
		return fasNew;
	}
	
	private void seed(int[] indices, int limit, int level) {
		if (level == indices.length) {
			this.doSeed(indices);
		} else {
			int start = (level == 0) ? 0 : indices[level - 1] + 1;
			for (indices[level] = start; indices[level] < this.numSNP; indices[level]++) {
				if (level == indices.length - 1 && indices[indices.length - 1] - indices[0] >= limit) {
					break;
				}
				seed(indices, limit, level + 1);
			} 
		}
	}
	
	private Set<Fragment> supportingFragments(String seed, int[] indices) {
		HashSet<Fragment> supportingFragments = new HashSet<Fragment>();
		for (Fragment f : this.fragments) {
			if (this.isMatch(seed, indices, f)) {
				supportingFragments.add(f);
			}
		}
		return supportingFragments;
	}

}
