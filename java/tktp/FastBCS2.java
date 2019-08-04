package tktp;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;
import java.util.stream.IntStream;

/**
 * This class provides a native Java implementation of Joe's TopKTaupath algorithm
 *
 * @author Adriano Caloiaro
 * @version 6/30/14
 */
public class FastBCS2 {
    final static boolean DEBUG = false;
    static int lastTieReset = Integer.MIN_VALUE;
    static Comparator comparator = Collections.reverseOrder();
    static private boolean PARALLELIZE = true;

    /**
     * ================================================================================================
     * Taupath functions
     * ================================================================================================
     */
    public static ConcordanceMatrixFBCS2 getPi(double[] x, double[] y, boolean parallelize) {
        PARALLELIZE = parallelize;
        return getPi(x, y);
    }

    // Create an ordered concordance matrix based on the fast backward conditional search algorithm
    public static ConcordanceMatrixFBCS2 getPi(double[] x, double[] y) {
        if (DEBUG) System.out.println("Performing fast backward conditional search ");

        final int N = x.length;
        double matrixSum;
        int[] qValuesi = new int[0];
        int[] qValuesk;

        // Create the concordance matrix with a naturally ordered permutation, i.e. 0 .. N-1. This is the base
        // from which all permuted concordance matrices are created.
        ConcordanceMatrixFBCS2 cm = ConcordanceMatrixFBCS2.create(x, y, PARALLELIZE);
        int[] pi = cm.pi();

        Vector<Integer>[] ties = new Vector[N];

        int stage = N - 1;
        lastTieReset = stage;

        // Calculate the column sums for the entire concordance matrix
        IntStream sumStream = IntStream.range(0, N);
        matrixSum = (PARALLELIZE) ? (double) sumStream.parallel().map(cm.columnSums(stage)).sum() :
                (double) sumStream.map(cm.columnSums(stage)).sum();

        // Initialize the first Taupath score
        IntStream.rangeClosed(0, stage).forEach(k -> cm.tauPath[k] = 1.0);
        cm.tauPath[stage] = (matrixSum / (N * (N - 1)));

        // Loops until the concordance matrix achieves full concordance
        while (matrixSum != (stage * (stage + 1))) {
            if (DEBUG) System.out.format("\n===============\ni=%s\n", stage);

            // Set the permutation index on cm
            cm.permute(pi);

            // The algorithm states that in the case of ties, a tie is selected randomly. By choosing the first,
            // we eliminate the element of choice and introduce determinism.
            Vector<Integer> tieList = cm.tieList(stage);
            Vector<Integer> ti = previousTies(ties, pi[stage], stage);

            int mini = tieList.firstElement();

            // Retain this tie list for future analysis when there is more than one
            if (tieList.size() > 1) ties[stage] = tieList;

            if (DEBUG) System.out.format("Transposing: i<%s> <-> %s %n", pi[stage], mini);
            permuteWhereValuesEqual(pi[stage], mini, pi, cm);

            // Subtract the least concordant column from the column sum list and return the resulting matrixSum
            // for the next stage: (m^{i-1})
            matrixSum = cm.subtractFromColSums(stage - 1, pi[stage]);

            // tie-check: Check if the column at pi[stage] tied with any other columns in a previous stage
            boolean swap = false;

            if (ti.size() > 0) {
                if (DEBUG) System.out.format("Previous ties: %s\n", ti.toString());
                int maxTieId = ti.firstElement();
                if (maxTieId > stage) {
                    for (int k : ti) {
                        if (k > stage) {
                            if (DEBUG) System.out.println("k = " + k);
                            // Generate the cumulative sums from rows i to k - i, inclusive, for columns i,k.
                            qValuesi = cm.cumulativeSums(stage, stage, k);
                            qValuesk = cm.cumulativeSums(k, stage, k);

                            // Test observations "k" and "stage" for stagewise optimality at the current stage
                            boolean allGreaterThanOrEqual = true;
                            boolean anyGreaterThan = false;
                            for (int z = 0; z < k - stage; z++) {
                                if (qValuesi[z] > qValuesk[z]) allGreaterThanOrEqual = false;
                                if (qValuesk[z] > qValuesi[z]) anyGreaterThan = true;
                            }

                            // Forward-swap: It has been determined that some other observation (pi[k]) that tied with
                            // column pi[stage] in a previous stage should be permuted to pi[stage] to ensure a locally
                            // monotone decreasing Taupath
                            if (allGreaterThanOrEqual && anyGreaterThan) {
                                if (DEBUG) System.out.format("pi: %s\n", Arrays.toString(pi));

                                // Recalculate the column sums between previousStage and k so matrixSum values are correct
                                // when taupath scores are calculated in the post-permute portion of the algorithm
                                for (int unacccountedForStage = stage; unacccountedForStage <= k; unacccountedForStage++) {
                                    matrixSum = cm.addToColSums(unacccountedForStage, pi[unacccountedForStage]);
                                }

                                // Swap the current stage observation with the previous tie that resulted in the greater Tau(k)
                                permuteWhereValuesEqual(pi[stage], pi[k], pi, cm);
                                stage = k - 1;

                                // Reset all ties that occurred prior to the current stage
                                for (int z = stage; z < N; z++) {
                                    ties[z] = null;
                                }

                                cm.tauPath[k] = matrixSum / (k * (k + 1));
                                matrixSum = cm.subtractFromColSums(stage, cm.piVal(k));

                                // Now we can calculate Tau(k) for the current stage
                                if (stage > 0) cm.tauPath[stage] = matrixSum / (stage * (stage + 1));

                                if (DEBUG)
                                    System.out.format("Reset: transpose i<%s> to %s\npi: %s\n", pi[stage], pi[k], Arrays.toString(pi));

                                swap = true;
                                break;
                            }
                        }
                    }
                }
            }// End of: tie-check
            // Post-permute
            // Decrement i to the next stage
            if (!swap) {
                // Calculate the TauPath score for the current stage
                if (stage > 0) cm.tauPath[stage - 1] = (matrixSum / (stage * (stage - 1)));

                stage--;
            }

            if (DEBUG) System.out.printf("Matrix sum: %s %n", matrixSum);
        } // End of: Backward elimination
        return cm;
    }

    /**
     * @param conc The concordance matrix
     * @return
     */
    public static double[] getTau(ConcordanceMatrixFBCS2 conc) {
        return conc.tauPath;
    }

    /**
     * Transposes values in the permutation index and column sums. This is a by-value transposition where k and l define the two
     * values in pi that will be swapped. This is likely NOT an ideal implementation and was created purely to
     * bring this class into parity with taupath_rex.R
     *
     * @param k  The index of the first value for which to search
     * @param l  The index of the second value for which to search
     * @param pi The permutation index in which transpositions will occur
     */
    private static void permuteWhereValuesEqual(int k, int l, int[] pi, ConcordanceMatrixFBCS2 cm) {
        int indexK = 0;
        int indexL = 0;

        // Find in the indices of k and l in pi
        for (int i = 0; i < pi.length; i++) {
            if (pi[i] == k) indexK = i;
            if (pi[i] == l) indexL = i;
        }

        // Permute the values
        int temp = pi[indexK];
        pi[indexK] = pi[indexL];
        pi[indexL] = temp;
        cm.pi(pi);

        // Permute the column sums in the column sum vector that correspond with the two columns that were just permuted
        temp = cm.data[cm.columnSumOffset + indexK];
        cm.data[cm.columnSumOffset + indexK] = cm.data[cm.columnSumOffset + indexL];
        cm.data[cm.columnSumOffset + indexL] = (short) temp;
    }

    /**
     * Resolves the naturally ordered indexes of the columns in the tie list and returns the maximum column
     *
     * @param ties The previous tie list to examine for "rowId"
     * @return The largest ID
     */
    public static Vector<Integer> previousTies(Vector<Integer>[] ties, int colId, int currentStage) {
        Vector<Integer> previousTies = new Vector<>();

        for (int i = currentStage+1; i <= lastTieReset; i++) {
            if (ties[i] != null) {
                if (ties[i].contains(colId)) previousTies.add(i);
            }
        }
        Collections.sort(previousTies, comparator);
        return previousTies;
    }
}
