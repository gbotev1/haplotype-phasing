package src.main;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;

import com.google.common.collect.TreeMultimap;

/**
 * This class allows the user to test MEC (Minimum Error Correction) used in an exponential algorithm.
 * 
 * @author gbotev
 */
public class Tester {
	
	// The user-supplied sequences to process.
	private static char[][] SNPMatrix;
	
	// The user-supplied sequences in short format to process.
	private static TreeMultimap<Integer, String> fragments;

	public static void main(String[] args) {
		//ReadSequencesFromFile("/Users/gbotev/Documents/GitHub Repositories/haplotype-phasing/scripts/mat.txt");
		//ProcessSNPMatrix();
		ReadNewSequencesFromFile("/Users/gbotev/Documents/GitHub Repositories/haplotype-phasing/scripts/mat.txt");
		ProcessShortFragMatrix();
	}

	/**
	 * This method determines the best-guess haplotype of the given SNP matrix under the MEC objective using a rudimentary exponential solver.
	 */
	private static void ProcessSNPMatrix() {
		ExponentialSolver exponentialSolver = new ExponentialSolver(SNPMatrix);
		long startTime = System.nanoTime();
		Haplotype haplotype = exponentialSolver.solve();
		long endTime = System.nanoTime();
		// Calculate the duration in microseconds
		long duration = (endTime - startTime) / 1000;
		System.out.println("Exponential Solver\nH1: " + haplotype.h1() + "\nH2: " + haplotype.h2() + "\nMEC: " + haplotype.MEC() + "\nTime: " + duration + " µs");
		/*
		for (String s : haplotype.getPartition()) {
			System.out.println(s.toString());
		}
		*/
	}
	
	/**
	 * This method determines the best-guess haplotype of the given SNP matrix under the MEC objective using the Fast Twister heuristic.
	 */
	private static void ProcessShortFragMatrix() {
		FastTwister fastSolver = new FastTwister(fragments);
		long startTime = System.nanoTime();
		Haplotype haplotype = fastSolver.solve();
		long endTime = System.nanoTime();
		// Calculate the duration in microseconds
		long duration = (endTime - startTime) / 1000;
		//System.out.println("Fast Twister Solver\nH1: " + haplotype.h1() + "\nH2: " + haplotype.h2() + "\nMEC: " + haplotype.MEC() + "\nTime: " + duration + " µs");
		System.out.println("Fast Twister Solver\nMEC: " + haplotype.MEC() + "\nTime: " + duration + " µs");
	}
	
	private static void ReadNewSequencesFromFile(String fileName) {
		try {
			// Load the raw file into a buffered reader
			BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
			
			// Initialize the TreeMap
			fragments = TreeMultimap.create();
			
			// Read the fragments from the file
			String currLine = null;
			StringBuilder index = new StringBuilder();
			while ((currLine = bufferedReader.readLine()) != null) {
				for (int i = 0; i < currLine.length(); i++) {
					char currChar = currLine.charAt(i);
					if (currChar == '\t') {
						// Add index and actual fragment
						fragments.put(Integer.parseInt(index.toString()), currLine.substring(i + 1));
					} else {
						index.append(currChar);
					}
				}
				// Clear the StringBuilder
				index.setLength(0);
			}
			
			// Close the bufferedReader when finished
			bufferedReader.close();     
        }
        catch(FileNotFoundException e) {
            System.out.println(String.format("Unable to open file %s.", fileName));                
        }
        catch(IOException e) {
            System.out.println(String.format("Error reading file %s.", fileName));
        }
	}
	
	/**
	 * This method reads the user-defined sequences from a file and stores them as arrays of integers.
	 * 
	 * @param fileName The absolute path of the file containing the user-defined sequences.
	 */
	private static void ReadSequencesFromFile(String fileName) {
        try {
			// Load the raw file into a buffered reader
			BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
			
			// Read the dimensions of the SNP matrix
			int rows = Integer.parseInt(bufferedReader.readLine());
			int columns = Integer.parseInt(bufferedReader.readLine());
			
			// Initialize the SNP matrix
			SNPMatrix = new char[rows][columns];
			for (int row = 0; row < rows; row++) {
				ProcessSNPMatrixRow(row, bufferedReader.readLine());
			}
			
			// Close the bufferedReader when finished
			bufferedReader.close();     
        }
        catch(FileNotFoundException e) {
            System.out.println(String.format("Unable to open file %s.", fileName));                
        }
        catch(IOException e) {
            System.out.println(String.format("Error reading file %s.", fileName));
        }
	}
	
	/**
	 * This method reads the user-defined sequences from the console and stores them as arrays of integers.
	 */
	private static void ReadSequencesFromConsole() {
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
