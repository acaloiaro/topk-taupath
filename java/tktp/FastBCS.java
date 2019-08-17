package tktp;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;

/**
 * This class provides a native Java implementation of Joe's TopKTaupath algorithm
 *
 * @author Adriano Caloiaro
 */
public class FastBCS {
  final static boolean DEBUG = false;
  static int lastTieReset = Integer.MIN_VALUE;
  static Comparator comparator = Collections.reverseOrder();

  // Create an ordered concordance matrix using FastBCS
  public static int[] getPi(double[] x, double[] y) {
    final int N = x.length;
    ConcordanceMatrix cm = ConcordanceMatrix.create(x, y);
    int[] pi = cm.pi();

    Vector<Integer>[] ties = (Vector<Integer>[]) new Vector[N];
    int i = N - 1;
    lastTieReset = i;
    boolean permute = true;

    while (permute) {
      Vector<Integer> tieList = cm.tieList(i);
      Vector<Integer> ti = previousTies(ties, pi[i], i);

      // The algorithm states that in the case of a tie, one is selected randomly. By choosing the first,
      // we eliminate the element of choice.
      Integer mini = tieList.elementAt(0);
      transposeWhereValuesEqual(pi[i], mini, pi, cm);
      if (tieList.size() > 1) ties[i] = tieList;

      boolean swap = false;

      if (ti.size() > 0) {
        int maxTieId = ti.firstElement();
        if (maxTieId > i) {
          for (int k : ti) {
            if (k > i) {
              int qValuesi[] = cm.cumulativeSums(i, i, k);
              int qValuesk[] = cm.cumulativeSums(k, i, k);

              boolean allGreaterThanOrEqual = true;
all:
              for (int z = 0; z < k - i; z++) {
                if (qValuesi[z] > qValuesk[z]) {
                  allGreaterThanOrEqual = false;
                  break all;
                }
              }

              boolean anyGreaterThan = false;
any:
              for (int z = 0; z < k - i; z++) {
                if (qValuesi[z] < qValuesk[z]) {
                  anyGreaterThan = true;
                  break any;
                }
              }

              if (allGreaterThanOrEqual && anyGreaterThan) {
                transposeWhereValuesEqual(pi[i], pi[k], pi, cm);
                i = k - 1;
                lastTieReset = i;
                swap = true;
                break;
              }
            }
          }
        }
      }

      if (!swap) {
        i--;
      }

      int matrixSum = cm.matrixSum(i);

      if (matrixSum == (i * (i + 1))) {
        permute = false;
      }
    }

    return cm.pi();
  }

  /**
   * Transposes values in the permutation index. This is a by-value transposition where k and l define the two
   * values in pi that will be swapped. This is likely NOT an ideal implementation and was created purely to
   * bring this class into parity with taupath_rex.R
   *
   * @param k  The first value for which to search
   * @param l  The second value for which to search
   * @param pi The permutation index in which transpositions will occur
   */
  private static void transposeWhereValuesEqual(int k, int l, int[] pi, ConcordanceMatrix cm) {
    int indexK = 0;
    int indexL = 0;

    // Find in the indices of k and l in pi
    for (int i = 0; i < pi.length; i++) {
      if (pi[i] == k) indexK = i;
      if (pi[i] == l) indexL = i;
    }

    // Transpose the values
    int temp = pi[indexK];
    pi[indexK] = pi[indexL];
    pi[indexL] = temp;
    cm.pi(pi);
  }

  /**
   * Resolves the naturally ordered indexes of the columns in the tie list and returns the maximum column
   *
   * todo eliminate use of Vector
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
