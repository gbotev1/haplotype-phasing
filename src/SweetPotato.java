import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.PriorityBlockingQueue;

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

	private int numThreads = 16;
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
	private boolean inclusiveSeeding;
	private static List<FrequencyArray> seedHaplotypes = new ArrayList<FrequencyArray>();
	private static PriorityBlockingQueue<FrequencyArrayPair> faPairs = new PriorityBlockingQueue<FrequencyArrayPair>(
			11, new Comparator<FrequencyArrayPair>() {
				@Override
				public int compare(FrequencyArrayPair fap1,
						FrequencyArrayPair fap2) {
					return -Double.compare(fap1.getScore(), fap2.getScore());
				}
			}); // Use default initial size

	public SweetPotato(Set<Fragment> fragments, int k, double alpha, int beta,
			int seedLength, int fragmentLength, boolean prettyPrint,
			boolean inclusiveSeeding) {
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
		// Save seeding flag
		this.inclusiveSeeding = inclusiveSeeding;
	}

	/**
	 * This method initializes the PriorityQueue with all of the initial pairs.
	 */
	private void initializeBestMerges() {
		// Parallelize for-loops
		ExecutorService executorService = Executors.newWorkStealingPool();
		// Save current number of seeds for efficiency
		int numSeedHaplotypes = seedHaplotypes.size();
		// Calculate all seedHaplotype pairs
		for (int i = 0; i < numSeedHaplotypes - 1; i++) {
			FrequencyArray fa1 = seedHaplotypes.get(i);
			double fa1Size = fa1.numSupportingFrags();
			Set<Fragment> fa1Frags = fa1.getFrags();
			// Create all FrequencyArray pairs
			for (int j = i + 1; j < numSeedHaplotypes; j++) {
				FrequencyArray fa2 = seedHaplotypes.get(j);
				double fa2Size = fa2.numSupportingFrags();
				// We do not have reason to believe that either
				// fa1 or fa2 will be smaller in size, so pass in
				// function arguments in any order
				double intersectionSize = Sets
						.intersection(fa1Frags, fa2.getFrags()).size();
				// Use Jaccard index
				double currIndex = intersectionSize
						/ (fa1Size + fa2Size - intersectionSize);
				// Only add if the intersection size is greater than zero;
				// this guarantees that the merge step will execute as expected
				if (intersectionSize > 0) {
					faPairs.add(new FrequencyArrayPair(fa1, fa2, currIndex));
				}
			}
		}
		// Do not forget to shutdown the executor service
		executorService.shutdown();
	}

	public void phaseSerial() {
		// Determine seeds
		System.err.println("Seeding");
		this.seed(new int[this.seedLength], this.fragmentLength, 0);
		// Remove conflicting fragments
		// seedHaplotypes =
		// FrequencyArray.removeConflictingFragments(seedHaplotypes);
		// Find the best initial merges
		System.err.println("Pairing");
		this.initializeBestMerges();
		// Merge best-guesses for seeds in parallel
		System.err.println("Merging");
		Collection<FrequencyArrayPair> toRemove = new HashSet<FrequencyArrayPair>();
		// Parallelize for-loops
		ExecutorService executorService = Executors.newWorkStealingPool();
		// Print the starting number of pairs
		System.err.printf("Starting number of pairs: %d\n", faPairs.size());
		while (!faPairs.isEmpty()) {
			System.err.println(faPairs.size());
			FrequencyArrayPair bestMerge = faPairs.poll();
			FrequencyArray fa1 = bestMerge.getFirst();
			FrequencyArray fa2 = bestMerge.getSecond();
			FrequencyArray merge = FrequencyArray.merge(fa1, fa2);
			if (merge != null) {
				// Merge was successful, so update seedHaplotypes
				seedHaplotypes.remove(fa1);
				seedHaplotypes.remove(fa2);
				// Update only the pairs that could have changed
				// Remove all FrequencyArrayPairs that were involved in merge
				for (FrequencyArrayPair fap : faPairs) {
					if (fap.contains(fa1) || fap.contains(fa2)) {
						toRemove.add(fap);
					}
				}
				faPairs.removeAll(toRemove);
				// Now, add all FrequencyArrayPairs
				Set<Fragment> mergedFrags = merge.getFrags();
				double mergedFragsSize = mergedFrags.size();
				for (FrequencyArray fa : seedHaplotypes) {
					executorService.submit(new Runnable() {
						@Override
						public void run() {
							Set<Fragment> currFrags = fa.getFrags();
							double currFragsSize = currFrags.size();
							double intersectionSize;
							if (mergedFragsSize <= currFragsSize) {
								intersectionSize = Sets
										.intersection(mergedFrags, currFrags)
										.size();
							} else {
								// currFragsSize < mergedFragsSize
								intersectionSize = Sets
										.intersection(currFrags, mergedFrags)
										.size();
							}
							// If intersection size is zero, merge will always
							// fail, so do not add!
							if (intersectionSize > 0) {
								// Intersection size is zero, but set score to
								// Jaccard index
								faPairs.add(new FrequencyArrayPair(merge, fa,
										intersectionSize / (currFragsSize
												+ mergedFragsSize
												- intersectionSize)));
							}
						}
					});
				}
				// Add merge to seedHaplotypes
				seedHaplotypes.add(merge);
				// Clear collection
				toRemove.clear();
			}
		}
		// Do not forget to shutdown the executor service
		executorService.shutdown();
		// Sort by SADF to print in convenient order
		Collections.sort(seedHaplotypes, new Comparator<FrequencyArray>() {
			public int compare(FrequencyArray fa1, FrequencyArray fa2) {
				// return -Integer.compare(fa1.length(), fa2.length());
				return -FrequencyArray.compareFrequencyArrays2(fa1, fa2);
			}
		});
		// Keep track of unique consensuses
		HashSet<String> uniqueConsensuses = new HashSet<String>();
		for (FrequencyArray fa : seedHaplotypes) {
			String currString = (!this.prettyPrint)
					? fa.consensus().print()
					: fa.consensus().prettyPrint();
			if (!uniqueConsensuses.contains(currString)) {
				System.out.printf("%s\n", currString);
				uniqueConsensuses.add(currString);
			}
		}
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
			if (!(indices[i] >= f.startIndex() && indices[i] <= f.endIndex()
					&& f.toString().charAt(indices[i] - f.startIndex()) == seed
							.charAt(i))) {
				return false;
			}
		}
		return true;
	}

	private boolean isValid(int[] indices, Fragment f) {
		for (int i = 0; i < indices.length; i++) {
			if (!(indices[i] >= f.startIndex() && indices[i] <= f.endIndex()
					&& f.toString()
							.charAt(indices[i] - f.startIndex()) != '-')) {
				return false;
			}
		}
		return true;
	}

	private TreeMap<String, Integer> seedFrequency(int[] indices) {
		HashMap<String, Integer> seedFrequencies = new HashMap<String, Integer>(
				indices.length);
		for (Fragment f : this.fragments) {
			if (this.isValid(indices, f)) {
				// Seed SNP sites are valid for current fragment, so add if
				// unique
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
					seedFrequencies.put(seed,
							seedFrequencies.get(seed) + f.frequency());
				}
			}
		}
		TreeMap<String, Integer> sortedSeedFrequencies = new TreeMap<String, Integer>(
				new ValueComparator(seedFrequencies));
		sortedSeedFrequencies.putAll(seedFrequencies);
		return sortedSeedFrequencies;
	}

	private void doSeed(int[] indices) {
		TreeMap<String, Integer> seedFrequency = this.seedFrequency(indices);
		// Make sure there is enough variability at this seed site
		int variabilityThreshold = (this.inclusiveSeeding) ? this.k - 1: this.k;
		if (seedFrequency.size() > variabilityThreshold) {
			Set<FrequencyArray> currSeedGroup = new HashSet<FrequencyArray>(
					this.k);
			int size = 0;
			for (int i = 0; i < this.k; i++) {
				FrequencyArray fa = new FrequencyArray(this.numSNP, this.tag);
				// For top k (already sorted) find consensus
				Map.Entry<String, Integer> currEntry = seedFrequency
						.pollFirstEntry();
				String currSeed = currEntry.getKey();
				Set<Fragment> currSupport = this.supportingFragments(currSeed,
						indices);
				fa.addFragment(currSupport);
				size += currSupport.size();
				currSeedGroup.add(fa);
			}
			// Process currSeedGroup here to remove errors via twisting step
			Set<FrequencyArray> result = currSeedGroup;
			while (true) {
				// Twist until we cannot twist anymore (size does not change)
				Set<FrequencyArray> result2 = this.twist(result, this.tag,
						size);
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
			seedHaplotypes.addAll(result);
			this.tag++;
		}
	}

	private Set<FrequencyArray> twist(Set<FrequencyArray> fas, int tag,
			int size) {
		// Extract all supporting fragments from current seed group
		Set<Fragment> supportingFrags = new HashSet<Fragment>(size + size / 2);
		for (FrequencyArray fa : fas) {
			supportingFrags.addAll(fa.getFrags());
		}

		ArrayList<Fragment> consensuses = new ArrayList<Fragment>(this.k);
		ArrayList<Set<Fragment>> fasNewRaw = new ArrayList<Set<Fragment>>(
				this.k);
		// Initialization step
		for (int i = 0; i < this.k; i++) {
			fasNewRaw.add(new HashSet<Fragment>());
		}
		// Calculate current consensus for efficiency
		for (FrequencyArray fa : fas) {
			consensuses.add(fa.consensus());
		}
		TreeMultiset<Map.Entry<Integer, Integer>> scores = TreeMultiset
				.create(Map.Entry.comparingByKey());
		for (Fragment f : supportingFrags) {
			for (int i = 0; i < this.k; i++) {
				int score = consensuses.get(i).similarTo(f);
				scores.add(
						new AbstractMap.SimpleImmutableEntry<Integer, Integer>(
								score, i));
			}
			Multiset.Entry<Map.Entry<Integer, Integer>> highScore = scores
					.lastEntry();
			if (highScore.getCount() == 1 && highScore.getElement()
					.getKey() >= f.length() - this.beta) {
				// DEFAULT: f.length() - 1
				// The membership of this fragment is unambiguous because it
				// matches 1 of the consensus
				// haplotypes and there are either no errors or 1 error (when
				// using magic number 8).
				fasNewRaw.get(highScore.getElement().getValue()).add(f);
			}
			scores.clear();
		}
		// Save for efficiency
		double cutoff = ((double) (supportingFrags.size())
				/ (this.alpha * this.k));
		for (int i = 0; i < this.k; i++) {
			if (fasNewRaw.get(i).size() <= cutoff) {
				// If one of the seeds does not have "enough" supporting
				// fragments
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
				if (level == indices.length - 1
						&& indices[indices.length - 1] - indices[0] >= limit) {
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
