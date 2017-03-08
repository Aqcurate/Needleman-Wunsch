import java.lang.Math;
import java.util.Arrays;

/**
 * This class uses the Needleman-Wunsch algorithm
 * to align nucleotide sequences.
 *
 * @author Andrew Quach
 * @author Tamir Enkhjargal
 *
 * @version 1.0.0
 */
public class NeedlemanWunsch {

    // Default scoring scheme for match, mismatch, and indel
    private final int MATCH;
    private final int MISMATCH;
    private final int INDEL;

    private String strand1;
    private String strand2;
    
    private boolean allowMismatch;

    public NeedlemanWunsch(String strand1, String strand2) {
        this(strand1, strand2, 1, -1, -1, true);
    }

    public NeedlemanWunsch(String strand1, String strand2, boolean allowMismatch) {
        this(strand1, strand2, 1, allowMismatch ? -1 : -999, -1, allowMismatch);
    }

    public NeedlemanWunsch(String strand1, String strand2, int match, int mismatch, int indel) {
        this(strand1, strand2, match, mismatch, indel, true);
    }

    public NeedlemanWunsch(String strand1, String strand2, int match, int mismatch, int indel, boolean allowMismatch) {
        this.strand1 = strand1;
        this.strand2 = strand2;

        this.MATCH = match;
        this.MISMATCH = mismatch;
        this.INDEL = indel;

        this.allowMismatch = allowMismatch;
    }

    /**
     * Generates solution matrix given 2 RNA strands.
     * Uses the Needleman-Wunsch algorithm.
     *
     * @return the solution matrix
     */
    public int[][] findSolution() {
        // Generate solution matrix based on lengths of both strands
        // Let strand1 be the side strand
        // Let strand2 be the top strand
        int[][] solution = new int[strand1.length()+1][strand2.length()+1];

        // Set the starting point to value of START
        solution[0][0] = 0;

        // Fill in the top row. Moving to the right always adds the value of INDEL.
        for (int i = 1; i < strand2.length()+1; i++) {
            solution[0][i] = solution[0][i-1] + INDEL;
        }

        // Fill in the left column. Moving down always adds the value of INDEL.
        for (int i = 1; i < strand1.length()+1; i++) {
            solution[i][0] = solution[i-1][0] + INDEL;
        }

        // Fill in the rest of the matrix based on a few rules.
        for (int i = 1; i < strand1.length()+1; i++) {
            for (int j = 1; j < strand2.length()+1; j++) {

                int matchValue;

                // If the characters that correspond to the grid position are equal for both strands
                // Set the matchValue to MATCH, else set the matchValue to MISMATCH
                if (strand1.charAt(i-1) == strand2.charAt(j-1)) matchValue = MATCH;
                else matchValue = MISMATCH;

                // Set the value to the maximum of these three values
                    // Position to the left + INDEL
                    // Position above + INDEL
                    // Position top-left + matchVALUE
                solution[i][j] = max(solution[i][j-1] + INDEL, solution[i-1][j] + INDEL,  solution[i-1][j-1] + matchValue);
            }
        }

        // Return solution matrix
        return solution;
    }

    /**
     * Helper method for calculating a maximum of three numbers.
     *
     * @return the maximum of the three given integers
     */
    private int max(int a, int b, int c) {
        return Math.max(Math.max(a, b), c);
    }

     /**
      * Aligns RNA strands based off a solution matrix.
      * Finds one of the 'best' paths in the solution matrix.
      * Uses the 'best' path to generate aligned RNA strands.
      * This method does so non-recursively.
      * 
      * @return the two aligned RNA strands
      */
    public String[] findPath(int[][] solution) {
        // Let strand1 be the side strand
        // Let strand2 be the top strand
        String alignedStrand1 = "";
        String alignedStrand2 = "";

        // Start from the bottom right of the solution matrix
        int i = solution.length - 1;
        int j = solution[0].length - 1;

        int best;
        boolean matchAllowed;

        // While you are not at the top/left side of the matrix
        // This prevents an OOB exception
        while (i != 0  && j != 0) {
            matchAllowed = true;
            if (strand1.charAt(i-1) != strand2.charAt(j-1) && !allowMismatch) matchAllowed = false;
            best = max(solution[i][j-1], solution[i-1][j], solution[i-1][j-1]);
            // Calculate the best path to the current position
            // If the top-left position is the best
            if (solution[i-1][j-1] == best && matchAllowed) {
                // Add the character corresponding to that position to both strands
                // This is the case for either a match or mismatch
                alignedStrand1 = strand1.charAt(i-1) + alignedStrand1;
                alignedStrand2 = strand2.charAt(j-1) + alignedStrand2;
                // Move to the new position
                i -= 1;
                j -= 1;
            // If the left position is the best
            } else if (solution[i][j-1] == best) {
                // Add '-' to strand1
                // Add the character correponding to that position to strand2
                // This represents a gap in the side strand
                alignedStrand1 = "-" + alignedStrand1;
                alignedStrand2 = strand2.charAt(j-1) + alignedStrand2;
                // Move to the new position
                j -= 1;
            // If the above position is the best
            } else {
                // Add '-' to strand2
                // Add the character corresponding to that position to strand1
                // This represents a gap in the top strand
                alignedStrand1 = strand1.charAt(i-1) + alignedStrand1;
                alignedStrand2 = "-" + alignedStrand2;
                // Move to the new position
                i -= 1;
            // If the top-left position is the best
            }
        }

        // If you are at the top of the matrix
        if (i == 0) {
            // Append "-" for every space you are away from 0,0 to strand2
            // EX: If you are at 0,3 (j = 3), add "---" to strand2
            for (int k = 0; k < j; k++) {
                alignedStrand1 = "-" + alignedStrand1;
                alignedStrand2 = strand2.charAt(j-k) + alignedStrand2;
            }
        // If you are at the left most side of the matrix
        } else {
            // Append "-" for every space you are away from 0,0 to strand1
            // EX: If you are at 3,0 (i = 3), add "---" to strand1
            for (int k = 0; k < i; k++) {
                alignedStrand1 = strand1.charAt(i-k) + alignedStrand1;
                alignedStrand2 = "-" + alignedStrand2;
            }
        }

        return new String[] {alignedStrand1, alignedStrand2};
    }
    
    public String[] recursiveFindPath(int[][] solution, int i, int j) {
        String alignedStrand1 = "";
        String alignedStrand2 = "";

        // If you are at the top of the matrix
        if (i == 0) {
            // Append "-" for every space you are away from 0,0 to strand1
            // EX: If you are at 0,3 (j = 3), add "---" to strand1
            for (int k = 0; k < j; k++) {
                alignedStrand1 = "-" + alignedStrand1;
                alignedStrand2 = strand2.charAt(j-k) + alignedStrand2;
            }

            return new String[] {alignedStrand1, alignedStrand2};
        // If you are at the left most side of the matrix
        } else if (j == 0) {
            // Append "-" for every space you are away from 0,0 to strand2
            // EX: If you are at 3,0 (i = 3), add "---" to strand1
            for (int k = 0; k < i; k++) {
                alignedStrand1 = strand1.charAt(i-k) + alignedStrand1;
                alignedStrand2 = "-" + alignedStrand2;
            }

            return new String[] {alignedStrand1, alignedStrand2};
        }

        // Calculate the best path to the current position
        // Check position to the left, above, and top-left
        int best;
        boolean matchAllowed = true;

        if (strand1.charAt(i-1) != strand2.charAt(j-1) && !allowMismatch) matchAllowed = false;
        best = max(solution[i][j-1], solution[i-1][j],  solution[i-1][j-1]);

        // If the top-left position is the best
        if (solution[i-1][j-1] == best && matchAllowed) {
            // Add the character corresponding to that position to both strands
            // This is the case for either a match or mismatch
            alignedStrand1 = "" + strand1.charAt(i-1);
            alignedStrand2 = "" + strand2.charAt(j-1);
            // Move to the new position
            i -= 1;
            j -= 1;
        // If the left position is the best
        } else if (solution[i][j-1] == best) {
            // Add '-' to strand1
            // Add the character correponding to that position to strand2
            // This represents a gap in the side strand
            alignedStrand1 = "-";
            alignedStrand2 = "" + strand2.charAt(j-1);
            // Move to the new position
            j -= 1;
        // If the above position is the best
        } else {
            // Add '-' to strand2
            // Add the character corresponding to that position to strand1
            // This represents a gap in the top strand
            alignedStrand1 = "" + strand1.charAt(i-1);
            alignedStrand2 = "-";
            // Move to the new position
            i -= 1;
        // If the top-left position is the best
        }

        String[] alignedStrands = recursiveFindPath(solution, i, j);
        alignedStrand1 = alignedStrands[0] + alignedStrand1;
        alignedStrand2 = alignedStrands[1] + alignedStrand2;

        return new String[] {alignedStrand1, alignedStrand2};
    }

    
    /**
     * Method that abstracts away findSolution and findPath.
     * Prints out the aligned strands and alignment score.
     */
    public void alignStrands() {
        int[][] solution = findSolution();
        System.out.println(Arrays.deepToString(solution));
        int score = solution[solution.length-1][solution[0].length-1];
        String[] alignedStrands = findPath(solution);

        System.out.println(alignedStrands[0]);
        System.out.println(alignedStrands[1]);

        alignedStrands = recursiveFindPath(solution, solution.length-1, solution[0].length-1);

        System.out.println("RECURSIVE FUN");
        System.out.println(alignedStrands[0]);
        System.out.println(alignedStrands[1]);

        System.out.println("The score for this alignment is: " + score);
    }
}
