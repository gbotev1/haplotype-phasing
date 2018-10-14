import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;

import com.google.common.collect.Sets;

/**
 * This class attempts a Priority-based merge using the given initial merge pair
 * and list of seed haplotypes.
 * 
 * @author gbotev
 */
public class Merger implements Runnable {
	
	private ArrayList<FrequencyArray> seedHaplotypes;
	private FrequencyArrayPair seedMerge;
	private PriorityQueue<FrequencyArrayPair> faPairs;
	private int numSNP;
	private boolean prettyPrint;
	
	public Merger(List<FrequencyArray> seedHaplotypes, FrequencyArrayPair seedMerge, int numSNP, boolean prettyPrint) {
		this.seedHaplotypes = new ArrayList<FrequencyArray>(seedHaplotypes);
		this.seedMerge = seedMerge;
		this.faPairs = new PriorityQueue<FrequencyArrayPair>(new Comparator<FrequencyArrayPair>() {
			@Override
			public int compare(FrequencyArrayPair fap1, FrequencyArrayPair fap2) {
				return -Double.compare(fap1.getScore(), fap2.getScore());
			}
		});
		this.faPairs.add(this.seedMerge);
		this.numSNP = numSNP;
		this.prettyPrint = prettyPrint;
	}
	
	@Override
	public void run() {
		// Merge now
		Collection<FrequencyArrayPair> toRemove = new HashSet<FrequencyArrayPair>();
		while (!this.faPairs.isEmpty()) {
			FrequencyArrayPair bestMerge = this.faPairs.poll();
			FrequencyArray fa1 = bestMerge.getFirst();
			FrequencyArray fa2 = bestMerge.getSecond();
			FrequencyArray merge = FrequencyArray.merge(fa1, fa2);
			if (merge != null) {
				// Merge was successful, so update seedHaplotypes
				this.seedHaplotypes.remove(fa1);
				this.seedHaplotypes.remove(fa2);
				// Update only the pairs that could have changed
				// Remove all FrequencyArrayPairs that were involved in merge
				for (FrequencyArrayPair fap : this.faPairs) {
					if (fap.contains(fa1) || fap.contains(fa2)) {
						toRemove.add(fap);
					}
				}
				this.faPairs.removeAll(toRemove);
				// Now, add all FrequencyArrayPairs
				Set<Fragment> mergedFrags = merge.getFrags();
				double mergedFragsSize = mergedFrags.size();
				for (FrequencyArray fa : this.seedHaplotypes) {
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
					// If intersection size is zero, merge will always fail, so do not add!
					if (intersectionSize > 0) {
						// Intersection size is zero, but set score to Jaccard index
						this.faPairs.add(new FrequencyArrayPair(merge, fa,
								intersectionSize / (currFragsSize + mergedFragsSize - intersectionSize)));
					}
				}
				// Add merge to seedHaplotypes
				this.seedHaplotypes.add(merge);
				// Clear collection
				toRemove.clear();
				// Check if desired haplotype length already achieved, and if so, stop merging!
				if (merge.length() == this.numSNP) {
					break;
				}
			}
		}
		FrequencyArray bestGuess = Collections.max(this.seedHaplotypes, new Comparator<FrequencyArray>() {
			public int compare(FrequencyArray fa1, FrequencyArray fa2) {
				return FrequencyArray.compareFrequencyArrays2(fa1, fa2);
            }
        });
		// Print in the format requested by the user
		if (this.prettyPrint) {
			bestGuess.consensus().prettyPrint();
		} else {
			bestGuess.consensus().print();
		}
	}
	
}
