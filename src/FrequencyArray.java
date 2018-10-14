import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Sets;

/**
 * This class efficiently keeps track of the frequencies at each SNP site.
 * 
 * @author gbotev
 *
 */
public class FrequencyArray {

	private int numSNP;
	private int activeStart;
	private int activeEnd;
	private int[] num0;
	private int[] num1;
	private Set<Fragment> supportingFragments;
	private Set<Integer> tags;
	
	public FrequencyArray(int numSNP) {
		this.numSNP = numSNP;
		// Initialize frequency arrays; by default, initialization occurs with 0
		this.num0 = new int[this.numSNP];
		this.num1 = new int[this.numSNP];
		// this.activeEnd = 0 by default
		this.activeStart = this.numSNP;
		this.supportingFragments = new HashSet<Fragment>();
		// Initialize set of tags
		this.tags = new HashSet<Integer>();
	}
	
	public FrequencyArray(int numSNP, int tag) {
		this.numSNP = numSNP;
		// Initialize frequency arrays; by default, initialization occurs with 0
		this.num0 = new int[this.numSNP];
		this.num1 = new int[this.numSNP];
		// this.activeEnd = 0 by default
		this.activeStart = this.numSNP;
		this.supportingFragments = new HashSet<Fragment>();
		// Initialize set of tags
		this.tags = new HashSet<Integer>();
		this.tags.add(tag);
	}
	
	public int length() {
		return this.activeEnd - this.activeStart + 1;
	}
	
	public int numSupportingFrags() {
		return this.supportingFragments.size();
	}
	
	public Set<Integer> getTags() {
		return tags;
	}
	
	public int getTag() {
		return tags.iterator().next();
	}
	
	@Override
	public boolean equals(Object o) {
		// Check if being compared to itself
	    if (this == o) {
	        return true;
	    }
	    // Check if object is null
	    if (o == null) {
	        return false;
	    }
	    // Check type of object
	    if (getClass() != o.getClass()) {
	        return false;
	    }
		// Typecast is safe to do now
		FrequencyArray fa = (FrequencyArray) o;
		// Check if supporting fragments are the same
		return Objects.equals(this.supportingFragments, fa.supportingFragments);
	}
	
	@Override
	public int hashCode() {
		return Objects.hash(this.supportingFragments);
	}
	
	public static int compareFrequencyArrays1(FrequencyArray fa1, FrequencyArray fa2) {
		return Integer.compare(fa1.activeStart, fa2.activeStart);
	}
	
	public static int compareFrequencyArrays2(FrequencyArray fa1, FrequencyArray fa2) {
		return Integer.compare(fa1.sadf(), fa2.sadf());
	}
	
	/**
	 * This method gets the supporting fragments for this FrequencyArray.
	 * @return The Set of supporting fragments.
	 */
	public Set<Fragment> getFrags() {
		return this.supportingFragments;
	}
	
	/**
	 * This method calculates the SADF score of the supporting fragments efficiently.
	 * @return The SADF score of the supporting fragments.
	 */
	public int sadf() {
		int sadf = 0;
		for (int i = 0; i < this.numSNP; i++) {
			sadf += Math.abs(this.num0[i] - this.num1[i]);
		}
		return sadf; 
	}
	
	/**
	 * This method calculates the MEC score of the supporting fragments efficiently.
	 * @return The MEC score of the supporting fragments.
	 */
	public int mec() {
		int mec = 0;
		for (int i = 0; i < this.numSNP; i++) {
			mec += Math.min(this.num0[i], this.num1[i]);
		}
		return mec;
	}
	
	/**
	 * This method removes the most conflicting fragments from the given set of seeds 
	 * until no more can be removed.
	 * @param fas The given set of haplotype seeds.
	 * @return The set of haplotype seeds with the conflicting fragments removed.
	 */
	public static List<FrequencyArray> removeConflictingFragments(List<FrequencyArray> fas) {
		int counter = 0;
		for (int prevSize = Integer.MAX_VALUE; prevSize != fas.size();) {
			// Keep track of the size
			prevSize = fas.size();
			Set<FrequencyArray> toRemoveFA = new HashSet<FrequencyArray>();
			for (int j = 0; j < prevSize; j++) {
				FrequencyArray fa = fas.get(j);
				Set<Fragment> toRemoveFrag = new HashSet<Fragment>();
				for (Fragment f : fa.supportingFragments) {
					String currFrag = f.toString();
					for (int i = 0; i < f.length(); i++) {
						int num0 = fa.num0[f.startIndex() + i];
						int num1 = fa.num1[f.startIndex() + i];
						char currSNP = currFrag.charAt(i);
						if (currSNP == '0') {
							if (num0 <= num1) {
								toRemoveFrag.add(f);
								break;
							}
						} else if (currSNP == '1') {
							if (num1 <= num0) {
								toRemoveFrag.add(f);
								break;
							}
						} // Else, currSNP == '-', so do nothing!
					}
				}
				for (Fragment f : toRemoveFrag) {
					counter++;
					fa.removeFragment(f);
				}
				if (fa.supportingFragments.isEmpty()) {
					toRemoveFA.add(fa);
				}
			}
			fas.removeAll(toRemoveFA);
		}
		System.err.printf("Removed %d fragments.\n", counter);
		return fas;
	}
	
	/**
	 * This method updates the frequency arrays by adding the given fragment.
	 * @param f The fragment to add.
	 */
	public void addFragment(Fragment f) {
		// Save fragment
		if (this.supportingFragments.add(f)) {
			// Update active region bounds
			if (f.startIndex() < this.activeStart) {
				this.activeStart = f.startIndex();
			}
			if (f.endIndex() > this.activeEnd) {
				this.activeEnd = f.endIndex();
			}
			// Update frequency arrays
			for (int i = 0; i < f.length(); i++) {
				if (f.toString().charAt(i) == '0') {
					this.num0[f.startIndex() + i] += f.frequency();
				} else if (f.toString().charAt(i) == '1') {
					this.num1[f.startIndex() + i] += f.frequency();
				} // Else, fragment.toString().charAt(i) == '-', so do nothing!
			}
		}
	}
	
	/**
	 * This method updates the frequency arrays by adding the given set of fragments.
	 * @param fragments The set of fragments to add.
	 */
	public void addFragment(Set<Fragment> fragments) {
		for (Fragment f : fragments) {
			this.addFragment(f);
		}
	}
	
	/**
	 * This method updates the frequency arrays by removing the given fragment.
	 * @param f The fragment to remove.
	 */
	public void removeFragment(Fragment f) {
		// Remove fragment
		if (this.supportingFragments.remove(f)) {
			// Update frequency arrays
			for (int i = 0; i < f.length(); i++) {
				if (f.toString().charAt(i) == '0') {
					this.num0[f.startIndex() + i] -= f.frequency();
				} else if (f.toString().charAt(i) == '1') {
					this.num1[f.startIndex() + i] -= f.frequency();
				} // Else, fragment.toString().charAt(i) == '-', so do nothing!
			}
			// Update active region bounds
			for (int i = 0; i < this.numSNP; i++) {
				if (this.num0[i] != 0 && this.num1[i] != 0) {
					// First SNP site with a contribution from at least one fragment
					this.activeStart = i;
					break;
				}
			}
			for (int i = this.numSNP - 1; i >= 0; i--) {
				if (this.num0[i] != 0 && this.num1[i] != 0) {
					// Last SNP site with a contribution from at least one fragment
					this.activeEnd = i;
					break;
				}
			}
		}
	}
	
	/**
	 * This method updates the frequency arrays by removing the given set of fragments.
	 * @param fragments The set of fragments to remove.
	 */
	public void removeFragment(Set<Fragment> fragments) {
		for (Fragment f : fragments) {
			this.removeFragment(f);
		}
	}
	
	/**
	 * This method merges two given FrequencyArray's together if they do not share
	 * any tags and returns a FrequencyArray with the result. Note that it assumes
	 * at least one supporting fragment is shared.
	 * @param fa1 The first FrequencyArray to merge.
	 * @param fa2 The second FrequencyArray to merge.
	 * @return The merged FrequencyArray, or null if the result was unsuccessful.
	 */
	public static FrequencyArray merge(FrequencyArray fa1, FrequencyArray fa2) {
		// Number of SNP sites should be the same for either FrequencyArray
		FrequencyArray mergedFA = new FrequencyArray(fa1.numSNP);
		// If no commons tags are shared, then merge! Since this merge method is being called, there will be
		// at least one supporting fragment shared.
		// We have reason to suspect that fa2 will be the smaller of the two sets, so use as first argument!
		if (Collections.disjoint(fa2.tags, fa1.tags)) {
			// Add both sets of supporting fragments
			mergedFA.addFragment(fa1.supportingFragments);
			mergedFA.addFragment(fa2.supportingFragments);
			// Add both sets of tags
			mergedFA.tags.addAll(fa1.tags);
			mergedFA.tags.addAll(fa2.tags);
			return mergedFA;
		}
		return null;
	}
	
	/**
	 * This method merges the given FrequencyArray with this FrequencyArray
	 * if they share any supporting fragments.
	 * 
	 * @param fa The FrequencyArray to merge.
	 * @return A boolean indicating whether the merge was successful.
	 */
	public boolean merge(FrequencyArray fa) {
		// If no commons tags are shared, and there is at least one common supporting fragment, then merge!
		if (Collections.disjoint(this.tags, fa.tags) && !Collections.disjoint(this.supportingFragments, fa.supportingFragments)) {
			this.addFragment(fa.supportingFragments);
			this.tags.addAll(fa.tags);
			return true;
		}
		return false;
	}
	
	/**
	 * This method calculates the consensus of the FragmentSet using the
	 * majority element for each SNP site and breaking ties with zero.
	 * 
	 * @return The consensus fragment.
	 */
	public Fragment consensus() {
		if (this.supportingFragments.isEmpty()) {
			return new Fragment(this.activeStart, "");
		}
		StringBuilder sb = new StringBuilder(this.activeEnd - this.activeStart + 1);
		// Remember that active region is inclusive at end points
		for (int i = this.activeStart; i <= this.activeEnd; i++) {
			if (this.num0[i] > this.num1[i]) {
				sb.append('0');
			} else if (this.num0[i] < this.num1[i]){
				sb.append('1');
			} else {
				// this.num0[i] == this.num1[i]
				sb.append('-');
			}
		}
		return new Fragment(this.activeStart, sb.toString());
	}
	
	/**
	 * This method calculates the consensus of the given set of fragments using
	 * the majority element for each SNP site and breaking ties with zero.
	 * 
	 * @param fragments
	 *            The given set of fragments for which to calculate the
	 *            consensus.
	 * @return The consensus fragment.
	 */
	public static Fragment consensus(Set<Fragment> fragments) {
		// Determine active region
		int min = Integer.MAX_VALUE;
		int max = 0;
		for (Fragment f : fragments) {
			if (f.startIndex() < min) {
				min = f.startIndex();
			}
			if (f.endIndex() > max) {
				max = f.endIndex();
			}
		}
		// Initialize consensus builder and frequency arrays
		int activeRegionLength = max - min + 1;
		StringBuilder sb = new StringBuilder(activeRegionLength);
		int[] num0 = new int[activeRegionLength];
		int[] num1 = new int[activeRegionLength];
		// Populate frequency arrays
		for (Fragment f : fragments) {
			for (int i = 0; i < f.length(); i++) {
				if (f.toString().charAt(i) == '0') {
					num0[f.startIndex() - min + i] += f.frequency();
				} else if (f.toString().charAt(i) == '1') {
					num1[f.startIndex() - min + i] += f.frequency();
				} // Else, f.toString().charAt(i) == '-', so do nothing!
			}
		}
		// Determine consensus
		for (int i = 0; i < activeRegionLength; i++) {
			if (num0[i] > num1[i]) {
				sb.append('0');
			} else if (num0[i] < num1[i]) {
				sb.append('1');
			} else {
				// this.num0[i] == this.num1[i]
				sb.append('-');
			}
		}
		return new Fragment(min, sb.toString());
	}
	
}
