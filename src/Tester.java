import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.TreeSet;

/**
 * This class allows the user to phase k-ploid haplotypes.
 * 
 * @author gbotev
 */
public class Tester {

	// The user-supplied sequences in short format to process.
	private static HashSet<Fragment> fragments;

	public static void main(String[] args) {
		try {
			// Extract commmand-line arguments
			String filename = args[0];
			int k = Integer.parseInt(args[1]);
			double alpha = Double.parseDouble(args[2]);
			int beta = Integer.parseInt(args[3]);
			int seedLength = Integer.parseInt(args[4]);
			int fragmentLength = Integer.parseInt(args[5]);
			boolean prettyPrint = Boolean.parseBoolean(args[6]);
			// Run phaser
			ReadNewSequencesFromFile(filename);
			//ProcessShortFragMatrixParallel(k, alpha, beta, seedLength, fragmentLength, prettyPrint);
			ProcessShortFragMatrixSerial(k, alpha, beta, seedLength, fragmentLength, prettyPrint);
		} catch (Exception e) {
			System.err.println("Command-line arguments were not entered properly.");
		}
	}
	
	/*
	private static void ProcessShortFragMatrixParallel(int k, double alpha, int beta, int seedLength, int numFragments, boolean prettyPrint) {
		System.err.println("Sweet Potato Solver Parallel");
		SweetPotato potatoSolver = new SweetPotato(fragments, k, alpha, beta, seedLength, numFragments, prettyPrint);
		long startTime = System.nanoTime();
		potatoSolver.phaseParallel();
		long endTime = System.nanoTime();
		// Calculate the duration in microseconds
		long duration = (endTime - startTime) / 1000;
		System.err.printf("Time: %d µs\nTime: %d ms\nTime: %d sec\n\n", duration, duration / 1000, duration / 1000000);
	}
	*/
	
	private static void ProcessShortFragMatrixSerial(int k, double alpha, int beta, int seedLength, int numFragments, boolean prettyPrint) {
		System.err.println("Sweet Potato Solver Serial");
		SweetPotato potatoSolver = new SweetPotato(fragments, k, alpha, beta, seedLength, numFragments, prettyPrint);
		long startTime = System.nanoTime();
		potatoSolver.phaseSerial();
		long endTime = System.nanoTime();
		// Calculate the duration in microseconds
		long duration = (endTime - startTime) / 1000;
		System.err.printf("Time: %d µs\nTime: %d ms\nTime: %d sec\n\n", duration, duration / 1000, duration / 1000000);
	}
	
	private static void ReadNewSequencesFromFile(String fileName) {
		try {
			// Load the raw file into a buffered reader
			BufferedReader bufferedReader = new BufferedReader(
					new FileReader(fileName));

			// Initialize the HashSet
			fragments = new HashSet<Fragment>();

			// Read the fragments from the file
			String currLine = null;
			StringBuilder index = new StringBuilder();
			while ((currLine = bufferedReader.readLine()) != null) {
				for (int i = 0; i < currLine.length(); i++) {
					char currChar = currLine.charAt(i);
					if (currChar == '\t') {
						// Save constants for efficiency
						int currIndex = Integer.parseInt(index.toString());
						String currString = currLine.substring(i + 1);
						Fragment currFragment = new Fragment(currIndex, currString);
						if (fragments.contains(currFragment)) {
							// Find the fragment and increment its frequency
							for (Fragment f : fragments) {
								if (f.equals(currFragment)) {
									f.incrementFrequency();
									break;
								}
							}
						} else {
							// fragments DOES NOT contain currFragment, so add
							// it!
							fragments.add(currFragment);
						}
					} else {
						index.append(currChar);
					}
				}
				// Clear the StringBuilder
				index.setLength(0);
			}
			// Close the bufferedReader when finished
			bufferedReader.close();
		} catch (FileNotFoundException e) {
			System.err.println(String.format("Unable to open file %s.", fileName));
		} catch (IOException e) {
			System.err.println(String.format("Error reading file %s.", fileName));
		}
	}

}
