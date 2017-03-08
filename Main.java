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
        // Strand1 = UUAGG
        // Strand2 = CGGCC
        // Do not allow mismatching
        NeedlemanWunsch n1 = new NeedlemanWunsch("UUAGG", "CGGCC", false);
        // Print out the information
        n1.printStrandInfo(); 
    }
}
