import java.util.Objects;

/**
 * This class allows FrequencyArrayPairs inserted into the faPairs in the
 * Solver class to maintain a consistent ordering. In particular, 
 * ties between FrequencyArrayPairs with equal values are broken via the
 * FIFO ordering.
 * 
 * Adapted from https://docs.oracle.com/javase/7/docs/api/java/util/concurrent/PriorityBlockingQueue.html.
 * 
 * @author Georgie Botev
 */
public class FIFOFrequencyArrayPair implements Comparable<FIFOFrequencyArrayPair> {
	// Current FrequencyArrayPair's number
	private final long seqNum;
	// The FrequencyArrayPair
	private final FrequencyArrayPair fap;
	
	public FIFOFrequencyArrayPair(FrequencyArrayPair fap, long seqNum) {
		this.seqNum = seqNum;
		this.fap = fap;
	}
	
	/**
	 * This method retrieves the FrequencyArrayPair stored in this object.
	 * @return The FrequencyArrayPair.
	 */
	public FrequencyArrayPair getFrequencyArrayPair() {
		return this.fap;
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
		FIFOFrequencyArrayPair FIFOfap = (FIFOFrequencyArrayPair) o;
		// Check if FrequencyArrayPairs are the same
		return Objects.equals(this.fap, FIFOfap.fap);
	}
	
	@Override
	public int hashCode() {
		return Objects.hash(this.fap);
	}
	
	/**
	 * Defines a custom comparator for these objects according to the
	 * FIFO ordering if the FrequencyArrayPair scores are equal.
	 */
	public int compareTo(FIFOFrequencyArrayPair other) {
		int result = -Double.compare(this.fap.getScore(), other.fap.getScore());
		// Check if tie
		if (result == 0) {
			// Break tie by FIFO ordering
			// Note that two seqNums will never be the same via
			// the way they are assigned!
			result = seqNum < other.seqNum ? -1 : 1;
		}
		return result;
	}
	
}
