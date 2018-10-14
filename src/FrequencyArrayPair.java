import java.util.Objects;

/**
 * This class encapsulates the notion of a pair of FrequencyArray's, including
 * useful information about the pair.
 * 
 * @author gbotev
 *
 */
public class FrequencyArrayPair {

	private FrequencyArray fa1;
	private FrequencyArray fa2;
	private double score;
	
	// Default constructor
	FrequencyArrayPair(FrequencyArray fa1, FrequencyArray fa2, double score) {
		this.fa1 = fa1;
		this.fa2 = fa2;
		this.score = score;
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
		FrequencyArrayPair fap = (FrequencyArrayPair) o;
		// Check if supporting fragments are the same
		return Objects.equals(this.fa1, fap.fa1) && Objects.equals(this.fa2, fap.fa2) && Objects.equals(this.score, fap.score);
	}
	
	@Override
	public int hashCode() {
		return Objects.hash(this.fa1, this.fa2, this.score);
	}
	
	/**
	 * Getter method for the score of this pair.
	 * @return The score of this pair.
	 */
	public double getScore() {
		return this.score;
	}
	
	/**
	 * Setter method to update the score of this pair.
	 * @param score The new score to write.
	 */
	public void setScore(int score) {
		this.score = score;
	}
	
	/**
	 * Getter method for the first FrequencyArray of this pair.
	 * @return The first FrequencyArray of this pair.
	 */
	public FrequencyArray getFirst() {
		return this.fa1;
	}
	
	/**
	 * Getter method for the second FrequencyArray of this pair.
	 * @return The second FrequencyArray of this pair.
	 */
	public FrequencyArray getSecond() {
		return this.fa2;
	}
	
	/**
	 * This method checks if this pair contains the given FrequencyArray.
	 * @param fa The FrequencyArray to check if contained.
	 * @return A boolean indicating whether this pair contains the given
	 * FrequencyArray.
	 */
	public boolean contains(FrequencyArray fa) {
		return fa.equals(this.fa1) || fa.equals(this.fa2);
	}

	/**
	 * Compares the scores of the given FrequencyArrayPair's favoring higher
	 * scores.
	 * @param fap1 The first FrequencyArrayPair.
	 * @param fap2 The second FrequencyArrayPair.
	 * @return An integer indicating whether the first score is greater than
	 *         (1), less than (-1), or equal to (0) the second score.
	 */
	/*
	public static int compareTo(FrequencyArrayPair fap) {
		return -Integer.compare(this.score, fap.score);
	}
	*/
	
}
