package src.main;

import java.util.Set;

/**
 * This class captures the haplotype as a pair of strings.
 * 
 * @author gbotev
 */
public class Haplotype {
	
	private String h1;
	private String h2;
	private int MEC;
	private Set<String> partition;
	
	public Haplotype(String h1, String h2, int MEC) {
		this.h1 = h1;
		this.h2 = h2;
		this.MEC = MEC;
	}
	
	public Set<String> getPartition() {
		return this.partition;
	}
	
	public void setPartition(Set<String> partition) {
		this.partition = partition;
	}
	
	public String h1() {
		return this.h1;
	}
	
	public String h2() {
		return this.h2;
	}
	
	public int MEC() {
		return this.MEC;
	}
	
}
