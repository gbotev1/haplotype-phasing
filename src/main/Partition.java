package src.main;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.TreeMultimap;

public class Partition {
	
	private HashSet<Map.Entry<Integer, String>> fragments;
	private int[] num0;
	private int[] num1;
	private int fragmentLength;
	
	public Partition(TreeMultimap<Integer, String> fragments, int[] num0, int[] num1, int fragmentLength) {
		// Set size of ArrayList for efficiency
		this.fragments = new HashSet<Map.Entry<Integer, String>>(fragments.size());
		// Populate HashSet from given TreeMultimap
		for (Map.Entry<Integer, String> entry : fragments.entries()) {
			this.fragments.add(entry);
		}
		// Save the given frequency arrays
		this.num0 = num0;
		this.num1 = num1;
		// Save the fragment length
		this.fragmentLength = fragmentLength;
	}
	
	public int[] get0Frequency() {
		return this.num0;
	}
	
	public int[] get1Frequency() {
		return this.num1;
	}
	
	/**
	 * Returns the sum of the absolute deviations of the frequency arrays.
	 * @param fragment
	 * @param baseline
	 * @return
	 */
	private int computeSADF(Map.Entry<Integer, String> fragment, int baseline) {
		int fragmentKey = fragment.getKey();
		String fragmentValue = fragment.getValue();
		int fragmentLength = fragmentValue.length();
		for (int i = 0; i < fragmentLength; i++) {
			if (fragmentValue.charAt(i) == '0') {
				if (this.num0[fragmentKey + i] <= this.num1[fragmentKey + i]) {
					baseline += 1;
				} else {
					// this.num0[fragmentKey + i] > this.num1[fragmentKey + i]
					baseline -= 1;
				}
			} else if (fragmentValue.charAt(i) == '1') {
				if (this.num0[fragmentKey + i] < this.num1[fragmentKey + i]) {
					baseline -= 1;
				} else {
					// this.num0[fragmentKey + i] >= this.num1[fragmentKey + i]
					baseline += 1;
				}
			} // Else, fragmentValue.charAt(i) == '-', so do nothing
		}
		return baseline;
	}
	
	public Set<Map.Entry<Integer, String>> getSwitchFragments() {
		// Initialize HashSet to keep track of smallest variation scores seen
		HashSet<Map.Entry<Integer, String>> switchFragments = new HashSet<Map.Entry<Integer, String>>();
		// Calculate SADF baseline for efficiency
		int SADFBaseline = 0;
		for (int i = 0; i < this.fragmentLength; i++) {
			SADFBaseline += Math.abs(this.num0[i] - this.num1[i]);
		}
		// Keep track of minimum variation score seen; allow for no fragments to be moved
		int maxSADF = SADFBaseline;
		// Calculate SADF scores for all fragments
		for (Map.Entry<Integer, String> entry : this.fragments) {
			int currVS = this.computeSADF(entry, SADFBaseline);
			if (currVS > maxSADF) {
				maxSADF = currVS;
				switchFragments.clear();
				switchFragments.add(entry);
			} else if (currVS == maxSADF) {
				switchFragments.add(entry);
			}
		}
		return switchFragments;
	}
	
	public void addFragments(Set<Map.Entry<Integer, String>> fragments) {
		// Adjust frequency arrays and add each fragment
		for (Map.Entry<Integer, String> entry : fragments) {
			int currKey = entry.getKey();
			String currFragment = entry.getValue();
			int currFragmentLength = currFragment.length();
			for (int i = 0; i < currFragmentLength; i++) {
				if (currFragment.charAt(i) == '0') {
					this.num0[currKey + i] += 1;
				} else if (currFragment.charAt(i) == '1') {
					this.num1[currKey + i] += 1;
				}
				// Else, currFragment.charAt(i) == '-', so do nothing
			}
			this.fragments.add(entry);
		}
	}
	
	public void removeFragments(Set<Map.Entry<Integer, String>> fragments) {
		// Adjust frequency arrays and remove each fragment
		for (Map.Entry<Integer, String> entry : fragments) {
			int currKey = entry.getKey();
			String currFragment = entry.getValue();
			int currFragmentLength = currFragment.length();
			for (int i = 0; i < currFragmentLength; i++) {
				if (currFragment.charAt(i) == '0') {
					this.num0[currKey + i] -= 1;
				} else if (currFragment.charAt(i) == '1') {
					this.num1[currKey + i] -= 1;
				}
				// Else, currFragment.charAt(i) == '-', so do nothing
			}
			this.fragments.remove(entry);
		}
	}

}
