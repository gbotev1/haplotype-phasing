package src.main;

public class MoveStatBundle {
	
	private int changeRegular;
	private int changeComplement;
	private int changeBoth;
	
	public MoveStatBundle(int changeRegular, int changeComplement, int changeBoth) {
		this.changeRegular = changeRegular;
		this.changeComplement = changeComplement;
		this.changeBoth = changeBoth;
	}
	
	public int getRegularChange() {
		return this.changeRegular;
	}
	
	public int getComplementChange() {
		return this.changeComplement;
	}
	
	public int getBothChange() {
		return this.changeBoth;
	}

}
