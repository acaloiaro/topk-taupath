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
public class YuFastBCS {
  final static boolean DEBUG = false;
  static int lastTieReset = Integer.MIN_VALUE;
  static Comparator comparator = Collections.reverseOrder();

  /**
   * ================================================================================================
   * Taupath functions
   * ================================================================================================
   */
  // Create an ordered concordance matrix based on the fast backward conditional search algorithm
  public static YuConcordanceMatrix getPi(double[] x, double[] y) {
    final int N = x.length;


    // Create the concordance matrix with a naturally ordered permutation, i.e. 0 .. N-1. This is the base
    // from which all permuted concordance matrices are created.
    YuConcordanceMatrix cm = YuConcordanceMatrix.create(x, y);
    int[] pi = cm.pi();
    Vector<Integer>[] ties = (Vector<Integer>[]) new Vector[N];
    int i = N - 1;
    lastTieReset = i;
    boolean permute = true;

    while (permute) {       // This while condition is suspect. The R implementation terminates
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

    return cm;
  }


  /**
   * todo make more efficient
   * @param conc The concordance matrix
   * @return
   */
  public static double[] getTau(YuConcordanceMatrix conc) {

    int scoreCount = 1;
    int n = conc.size();

    double[] tauPathScores = new double[n];
    for(int pathlength = 2; pathlength <= n; pathlength++){
      short[] partialMatrix = new short[((pathlength * (pathlength - 1)) / 2) + pathlength];  // n(n-1)/2  ( Lower triangular matrix, without diagonal )

      for(int i = 0; i < pathlength; i++) {
        //System.out.print(i + ": ");
        for(int j = 0; j < i; j++){
          short val = conc.pv(i, j);
          //System.out.print(val + ", ");
          partialMatrix[((i * (i + 1)) / 2) + j] = val;
        }
        //System.out.println();
      }
      //if (pathlength == 66) System.out.println(Arrays.toString(partialMatrix));
      tauPathScores[scoreCount++] = calculateKendallTauCoefficient(partialMatrix, pathlength);
    }

    tauPathScores[0] = tauPathScores[1]; // Apply penalty to tauPathScores[0] ( todo get reference from topK paper )

    //System.out.println(Arrays.toString(tauPathScores));
    return tauPathScores;
  }

  private static double calculateKendallTauCoefficient(short[] concordanceMatrix, int n){
    int sum             = 0;
    int nosOffDiagonals = 0;

    for(int i = 0; i < n; i++){
      for(int j = 0; j < i; j++){
        //if (n == 66) { System.out.print(concordanceMatrix[((i * (i + 1)) / 2) + j] + ", " + j);  }
        sum += concordanceMatrix[((i * (i + 1)) / 2) + j];

        nosOffDiagonals++;
      }
      //if (n==66) System.out.println();
    }
    return ((double)sum/nosOffDiagonals);
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
  private static void transposeWhereValuesEqual(int k, int l, int[] pi, YuConcordanceMatrix cm) {
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

  public static void main(String[] a) {
      double[] x = { 1, 2, 3, 4 };
      double[] y = { 1, 2, 3, 4 };
      YuConcordanceMatrix m = getPi(x, y);
      System.out.println(m.toTable());
  }
}
