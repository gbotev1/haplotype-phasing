package src.main;

/**
 * This class keeps track of the frequency of each of the possible entries at a particular SNP site.
 * 
 * @author gbotev
 */
public class SNPFrequency {

	private int num0;
	private int num1;
	
	public SNPFrequency() {
		this.num0 = 0;
		this.num1 = 0;
	}
	
	public char maxFrequency() {
		// Break ties with zero
		if (num0 >= num1) {
			return '0';
		} else {
			// There are  more ones
			return '1';
		}
	}
	
	public void increment0() {
		num0++;
	}
	
	public int num0() {
		return num0;
	}
	
	public void increment1() {
		num1++;
	}
	
	public int num1() {
		return num1;
	}
	
}
