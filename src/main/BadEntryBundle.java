package src.main;

import java.util.Map;

public class BadEntryBundle {
	
	private Map.Entry<Integer, String> badEntry1;
	private Map.Entry<Integer, String> badEntry2;
	private MoveStatBundle bundle;
	
	public BadEntryBundle(Map.Entry<Integer, String> badEntry1, Map.Entry<Integer, String> badEntry2, MoveStatBundle bundle) {
		this.badEntry1 = badEntry1;
		this.badEntry2 = badEntry2;
		this.bundle = bundle;
	}
	
	public Map.Entry<Integer, String> getBadEntry1() {
		return this.badEntry1;
	}
	
	public Map.Entry<Integer, String> getBadEntry2() {
		return this.badEntry2;
	}
	
	public MoveStatBundle getMoveStats() {
		return this.bundle;
	}

}
