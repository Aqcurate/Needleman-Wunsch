import java.util.Arrays;
/**
 * This class is an example runner
 * for NeedlemanWunsch.java
 *
 * @author Andrew Quach
 * @author Tamir Enkhjargal
 *
 * @version 2.0.0
 */
public class Main {

    /**
     * Method that creates a NeedlemanWunsch object and
     * prints out the Strand information.
     */
    public static void main(String[] args) {
        // Create a NeedlemanWunsch object
        // Strand1 = CGUCC 
        // Strand2 = GCCC
        // Match = 5
        // Mismatch = -3
        // Indel = -5
        // Allow mismatching
        NeedlemanWunsch n1 = new NeedlemanWunsch("CGUCC", "GCCC", 5, -3, -5, true);
        // Print out the information
        n1.printStrandInfo(); 

        NeedlemanWunsch n2 = new NeedlemanWunsch("CGUCC", "GCCC", false);
        n2.printStrandInfo();
    }
}
