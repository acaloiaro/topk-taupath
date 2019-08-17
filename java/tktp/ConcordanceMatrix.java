package tktp;

import java.util.Arrays;
import java.util.Vector;
import java.util.function.IntConsumer;
import java.util.stream.IntStream;

/**
 * The concordance matrix is an nXn matrix where n is the number of observations in a bivariate distribution.
 * This implementation is a square, immutable, efficient in-memory representation backed by a single dimension
 * array of bytes in row-major order.
 *
 * Optimizations:
 * 1. All data is stored in a contiguous, row-major array with the following structure:
 * [
 *      [ permutation index  c1, c2, c3, … cn ],
 *      [ column sum         c1, c2, c3, … cn ],
 *      [ concordance       (r1,c1), (r1,c2), (r1,c3), (r1,cn) … ],
 *      …
 *      [ concordance       (rn,c1), (rn,c2), (rn,c3), (rn,cn) … ]
 * ]
 * 2. Concordance matrix was switched back to row-major ordering from ConcordanceMatrixStreams2. When the permutation
 * index, columns sums, and concordance data are contained in a single, linear matrix, the design favors row-major
 * ordering.

 *
 * @author Adriano Caloiaro
 */
public final class ConcordanceMatrix {
 private int size;
 public short[] data;
 public int columnSumOffset;
 private int dataOffset;

 private boolean DEBUG = false;
 public static ConcordanceMatrix create(double[] x, double[] y) {
  final int N = x.length;

  ConcordanceMatrix m = new ConcordanceMatrix();
  m.size = N;

  // Set the offsets
  m.columnSumOffset = N;
  m.dataOffset = 2 * N;

  // calculate the concordance matrix
  m.data = new short[(2 * N) + (N * N)]; // 2 * N for the permutation index column and the column sum column
  IntConsumer lambda = gid -> {
   int i = gid / N;
   int j = gid % N;

   if (i != j) {
    double product = ((y[i] - y[j]) * (x[i] - x[j]));

    if (product > 0) {
     m.setV(i, j, (byte) 1);
    } else if (product < 0) {
     m.setV(i, j, (byte) - 1);
    }
   }
  };

  // Calculate the concordance matrix either sequentially or in parallel based on job size
  IntStream.range(0, N * N).forEach(lambda);

  // initialize the base matrix with a naturally ordered permutation index, i.e. 0 .. N-1
  m.pi(new int[N]);

  // Initialize the natural ordering of pi values
  IntStream.range(0, N).forEach(i -> m.setPiAt(i, i)); // Cannot assign pi value with a JTP device because of the even/divisible by 64 bug

  return m;
 }

 /** get the value of a cell */
 public short v(int i, int j) {
  return this.data[dataOffset + (i * size + j)];
 }

 /** get the value of a cell through its permuted index */
 public short pv(int i, int j) {
  return this.data[dataOffset + (data[i] * size) + data[j]];
 }

 /** set the value of a cell */
 private void setV(int i, int j, byte value) {
  this.data[dataOffset + (i * size) + j] = value;
 }

 /** the size of one dimension of the matrix */
 public int size() {
  return this.size;
 }

 /** permutation index */
 public int[] pi() {
  return IntStream.range(0, size).map(i -> this.data[i]).toArray();
 }

 public void pi(int[] pi) {
  for (int i = 0; i < size; i++) {
   this.data[i] = (short) pi[i];
  }
 }

 public int piVal(int i) {
  return this.data[i];
 }

 public void setPiAt(int i, int value) {
  this.data[i] = (short) value;
 }

 public int[] cumulativeSums(int j, int high, int low) {
  int[] sums = new int[low - high + 1];
  int[] tmpCol = new int[low + 1];


  for (int i = 0; i < tmpCol.length; i++) {
   tmpCol[i] = pv(i, j);
  }

  // Swap when [start] falls on the diagonal. Introduced in tKtp (2014-12-12)
  if (tmpCol[low] == 0) {
   tmpCol[low] = tmpCol[high];
   tmpCol[high] = 0;
  }

  // Accumulate all sums in the range [0,high) in sums[0]
  int sumIndex = 0;
  for (int i = 0; i <= high; i++) {
   sums[sumIndex] += tmpCol[i];
  }

  // Accumulate all sums in the range [0, i) in sums[i]
  for (int i = high + 1; i <= low; i++) {
   sums[++sumIndex] = sums[sumIndex - 1] + tmpCol[i];
  }

  return sums;
 }

 public int matrixSum(int end) {
  int sum = 0;
  columnSums(end);
  for (int j = 0; j <= end; j++) {
   sum += this.data[columnSumOffset + j];
  }
  return sum;
 }

 /** Returns the sum of all columns from 0 to index, inclusive. */
 public void columnSums(int index) {
  for (int i = 0; i <= index; i++) {
   short sum = 0;
   for (int j = 0; j <= index; j++) {
    sum += pv(i, j);
   }
   this.data[columnSumOffset + i] = sum;
  }
 }

 /** Create a new concordance matrix from the permuted index */
 public ConcordanceMatrix permute(int[] pi) {
  this.pi(pi);
  return this;
 }

 /** Find the column(s) with the minimum sum */
 public Vector < Integer > tieList(int index) {
  Vector < Integer > listOfTies = new Vector < > ();

  columnSums(index);

  // find the minimum sum
  int minColumnSum = Integer.MAX_VALUE;
  for (int j = 0; j <= index; j++) {
   if (data[columnSumOffset + j] < minColumnSum) {
    minColumnSum = data[columnSumOffset + j];
   }
  }

  // find the columns with the minimum sum
  for (int j = 0; j <= index; j++) {
   if (data[columnSumOffset + j] == minColumnSum) {
    listOfTies.addElement((int) data[j]);
   }
  }

  if (DEBUG) System.out.format("tie list: %s\n", listOfTies.toString());
  return listOfTies;
 }
}

