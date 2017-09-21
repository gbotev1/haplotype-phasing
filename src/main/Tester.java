package src.main;

import java.util.Scanner;

/**
 * This class allows the user to test MEC (Minimum Error Correction) used in an exponential algorithm.
 * 
 * @author gbotev
 */
public class Tester {
	
	// The user-supplied sequences to process.
	private static char[][] SNPMatrix;

	public static void main(String[] args) {
		ReadSequences();
		ProcessSNPMatrix();
	}
	
	/**
	 * This method determines the best-guess for the haplotype of the given SNP matrix based on MEC.
	 */
	private static void ProcessSNPMatrix() {
		MECExponentialSolver exponentialSolver = new MECExponentialSolver(SNPMatrix);
		long startTime = System.nanoTime();
		Haplotype haplotype = exponentialSolver.solve();
		long endTime = System.nanoTime();
		// Calculate the duration in microseconds
		long duration = (endTime - startTime) / 1000;
		System.out.println("Exponential MEC Solver\nH1: " + haplotype.h1() + "\nH2: " + haplotype.h2() + "\nMEC: " + haplotype.MEC() + "\nTime: " + duration + " µs");
		for (String s : haplotype.getPartition()) {
			System.out.println(s.toString());
		}
	}
	
	/**
	 * This method reads the user-defined sequences and stores them as arrays of integers.
	 */
	private static void ReadSequences() {
		Scanner in = new Scanner(System.in);
		// Query user for dimensions of SNP matrix
		System.out.println("Enter the number of rows in the SNP matrix.");
		int rows = in.nextInt();
		System.out.println("Enter the number of columns in the SNP matrix.");
		int columns = in.nextInt();
		// Initialize the SNP matrix
		SNPMatrix = new char[rows][columns];
		System.out.println("Input the SNP matrix row by row (no spaces between individual entries).");
		int row = 0;
		while(row < rows && in.hasNextLine()) {
			String nextLine = in.nextLine();
			// Disregard first pass empty line
			if (nextLine.isEmpty()) {
				continue;
			}
			ProcessSNPMatrixRow(row, nextLine);
			row++;
		}
		in.close();
	}
	
	/**
	 * This method adds the given row to the SNP matrix.
	 * @param rowNum The current row index.
	 * @param row The given row to process.
	 */
	private static void ProcessSNPMatrixRow(int rowNum, String row) {
		for (int i = 0; i < row.length(); i++) {
			SNPMatrix[rowNum][i] = row.charAt(i);
		}
	}
	
}
