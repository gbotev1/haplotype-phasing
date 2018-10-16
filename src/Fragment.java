import java.util.Objects;

/**
 * This class contains all the properties necessary to represent a fragment
 * utilizing the short encoding.
 * 
 * @author gbotev
 *
 */
public class Fragment {

	private int startIndex;
	private int endIndex;
	private String fragment;
	private int frequency;
	private int length;

	public Fragment(int startIndex, String fragment) {
		this.startIndex = startIndex;
		this.fragment = fragment;
		this.frequency = 1;
		this.length = fragment.length();
		this.endIndex = startIndex + this.length - 1;
	}
	
	public static int compareFragments(Fragment f1, Fragment f2) {
		return Integer.compare(f1.startIndex, f2.startIndex);
	}
		
	/**
	 * This method increments the frequency of the fragment.
	 */
	public void incrementFrequency() {
		this.frequency++;
	}
	
	public int similarTo(Fragment f) {
		// Determine overlapping range (i.e. active region)
		int start = Math.max(this.startIndex, f.startIndex);
		int end = Math.min(this.endIndex, f.endIndex);
		// Initialize counter for number of matches between the fragments
		int similarityScore = 0;
		// Remember that bounds are INCLUSIVE
		for (int i = start; i <= end; i++) {
			char c1 = this.fragment.charAt(i - this.startIndex);
			char c2 = f.fragment.charAt(i - f.startIndex);
			if (c1 != '-' && c2 != '-' && c1 == c2) {
				// Match at this position
				similarityScore += 1;
			}
		}
		return similarityScore;
	}
	
	@Override
	public String toString() {
		return this.fragment;
	}
	
	public String prettyPrint() {
		return String.format("%d-(%d)->%s", this.startIndex, this.length, this.fragment);
	}
	
	public String print() {
		return String.format("%d\t%s", this.startIndex, this.fragment);
	}
	
	@Override
	public int hashCode() {
		return Objects.hash(this.startIndex, this.fragment);
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
		Fragment f = (Fragment) o;
		// Check if starting index and fragment are the same
		return Objects.equals(this.startIndex, f.startIndex)
	            && Objects.equals(this.fragment, f.fragment);
	}

	/**
	 * This method gets the length of the fragment efficiently.
	 * 
	 * @return The length of the fragment.
	 */
	public int length() {
		return this.length;
	}

	/**
	 * This method gets the number of times the fragment has been repeated.
	 * 
	 * @return The frequency of the fragment.
	 */
	public int frequency() {
		return this.frequency;
	}

	/**
	 * This method gets the starting index of the fragment.
	 * 
	 * @return The starting index of the fragment.
	 */
	public int startIndex() {
		return this.startIndex;
	}

	/**
	 * This method gets the ending index of the fragment.
	 * 
	 * @return The ending index of the fragment.
	 */
	public int endIndex() {
		return this.endIndex;
	}

}
