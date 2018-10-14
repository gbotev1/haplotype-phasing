

public class Bundle {
	
	private int MEC;
	private String haplotype;
	
	public Bundle(int MEC, String haplotype) {
		this.MEC = MEC;
		this.haplotype = haplotype;
	}
	
	public int getMEC() {
		return this.MEC;
	}
	
	public String getHaplotype() {
		return this.haplotype;
	}

}
