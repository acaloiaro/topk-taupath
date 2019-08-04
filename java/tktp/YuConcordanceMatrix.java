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

public final class YuConcordanceMatrix {
    private int    size;
	public short[] data;
    public int columnSumOffset;
    private int dataOffset;

	private boolean DEBUG = false;

    //////////////////////////////////////////////////////////////////////////////
	// Constructors
	//////////////////////////////////////////////////////////////////////////////

	/** Prohibit public use of default constructor. */
	private YuConcordanceMatrix() {
    }

    /**
     * Returns the permuted data array
     *
     * @return The data
     */
    public short[] data() {
        short[] data = new short[size*size];

         IntStream.range(0, size*size).forEach(k -> {
             int i = k / size;
             int j = k % size;
             data[k] = pv(i,j);
        });

        return data;
    }
    public static YuConcordanceMatrix create(int[] data, int size) {
        final int N = size;

        YuConcordanceMatrix m = new YuConcordanceMatrix();
        m.size = N;

        // Set the offsets
        m.columnSumOffset = N;
        m.dataOffset = 2 * N;

        // calculate the concordance matrix
        m.data = new short[(2 * N) + (N * N)]; // 2 * N for the permutation index column and the column sum column

        IntConsumer lambda = gid -> {
            int i = gid / N;
            int j = gid % N;

            m.setV(i, j, (byte)data[i * size + j]);
        };

        // Calculate the concordance matrix either sequentially or in parallel based on job size
        IntStream.range(0,N*N).forEach(lambda);

        // initialize the base matrix with a naturally ordered permutation index, i.e. 0 .. N-1
        m.pi(new int[N]);

        // Initialize the natural ordering of pi values
        IntStream.range(0,N).forEach(i -> m.setPiAt(i, i));   // Cannot assign pi value with a JTP device because of the even/divisible by 64 bug

        return m;
    }

	public static YuConcordanceMatrix create(double[] x, double[] y) {
		final int N = x.length;

		YuConcordanceMatrix m = new YuConcordanceMatrix();
		m.size = N;

        // Set the offsets
        m.columnSumOffset = N;
        m.dataOffset = 2*N;

		// calculate the concordance matrix
		m.data = new short[(2*N) + (N * N)]; // 2 * N for the permutation index column and the column sum column

        IntConsumer lambda = gid -> {
            int i = gid / N;
            int j = gid % N;

            if (i != j)  {
                double product = ((y[i] - y[j]) * (x[i] - x[j]));

                if (product > 0) {
                    m.setV(i, j, (byte) 1);
                } else if (product < 0) {
                    m.setV(i, j, (byte) -1);
                }
            }
        };

        // Calculate the concordance matrix either sequentially or in parallel based on job size
        IntStream.range(0,N*N).forEach(lambda);

        // initialize the base matrix with a naturally ordered permutation index, i.e. 0 .. N-1
        m.pi(new int[N]);

        // Initialize the natural ordering of pi values
        IntStream.range(0,N).forEach(i -> m.setPiAt(i, i));   // Cannot assign pi value with a JTP device because of the even/divisible by 64 bug

		return m;
	}

	//////////////////////////////////////////////////////////////////////////////
	// API
	//////////////////////////////////////////////////////////////////////////////

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
        return IntStream.range(0,size).map(i -> this.data[i]).toArray();
	}

    public void pi(int[] pi) {
        for(int i = 0; i < size; i++){
            this.data[i] = (short)pi[i];
        }
    }

    public int piVal(int i) {
        return this.data[i];
    }

    public void setPiAt(int i, int value) {
        this.data[i] = (short)value;
    }

	public int[] cumulativeSums(int j, int high, int low) {
		int[] sums = new int[low - high + 1];
        int[] tmpCol = new int[low+1];


        for(int i = 0; i < tmpCol.length; i++) {
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
		for (int i = high+1; i <= low; i++) {
			sums[++sumIndex] = sums[sumIndex - 1] + tmpCol[i];
		}

		return sums;
	}

    public int matrixSum(int end) {
        int sum = 0;
        columnSums(end);
        for(int j = 0; j <= end; j++) {
            sum += this.data[columnSumOffset + j];
        }
        return sum;
    }



    /** Returns the sum of all columns from 0 to index, inclusive. */
    public void columnSums(int index) {
        for (int i = 0; i <= index; i++) {
            short sum = 0;
            for (int j = 0; j <= index; j++) {
                sum += pv(i,j);
            }
            this.data[columnSumOffset + i] = sum;
        }
    }

	/** Create a new concordance matrix from the permuted index */
	public YuConcordanceMatrix permute(int[] pi) {
		this.pi(pi);
		return this;
	}

	// Calculate a Kendall tau coefficient which is the average value of all off-diagonal values within a
	// concordance matrix. Since the matrix is reflective, the lower triangle can be used instead. The
	// tau path scores are generated by using partial matrices from 2 to the size of the concordance.
	// This implementation is O(n).
	public double[] tauPathScores() {
		double[] tauPathScores = new double[this.size() - 1];
		double sum = 0;
		double nosOffDiagonals = 0;

		for (int i = 1; i < this.size(); i++) {
			for (int j = 0; j < i; j++) {
				sum += pv(i, j);
				nosOffDiagonals++;
			}
			double score = sum / nosOffDiagonals;
			tauPathScores[i - 1] = score;
		}
		return tauPathScores;
	}

    /** Find the column(s) with the minimum sum */
    public Vector<Integer> tieList(int index) {
        Vector<Integer> listOfTies = new Vector<>();

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
                listOfTies.addElement((int)data[j]);
            }
        }

        if (DEBUG) System.out.format("tie list: %s\n", listOfTies.toString());
        return listOfTies;
    }



	/** Create a printable table of the matrix values */
	public String toTable() {
		StringBuffer sb = new StringBuffer("CONCORDANCE MATRIX\n");
		for (int i = 0; i < this.size; i++) {
			for (int j = 0; j < this.size; j++) {
				short v = pv(i, j);
				if (v == -1)
					sb.append(" " + v);
				else
					sb.append("  " + v);
			}
			sb.append('\n');
		}
        sb.append("pi:         " + Arrays.toString(pi()) + "\n");
        sb.append(("Column sums: "));
        IntStream.range(0, size).forEach(j -> sb.append(this.data[columnSumOffset + j] + ", "));
		return sb.toString();
	}

    /**
     * Create a printable table of the matrix values
     */
    public String toTable(int colId) {
        StringBuffer sb = new StringBuffer("CONCORDANCE MATRIX\n");
        for (int i = 0; i <= colId; i++) {
            for (int j = 0; j <= colId; j++) {
                short v = pv(i, j);
                if (v == -1)
                    sb.append(" " + v);
                else
                    sb.append("  " + v);
            }
            sb.append('\n');
        }
        sb.append(("Column sums: "));
        IntStream.rangeClosed(0, colId).forEach(j -> sb.append(this.data[columnSumOffset + j] + ", "));
        sb.append("\npi:         " + Arrays.toString(pi()) + "\n");

        return sb.toString();
    }

    /**
     * Prints the column sum vector up to the stage "stage"
     */
    public void printColumnSums(int stage, int forwardStage) {
        int extra = forwardStage - stage;
        System.out.format("Column Sums [%d]: ", stage);
        for(int j = 0; j <= stage + extra; j++) {
            if (j == stage) System.out.print(this.data[columnSumOffset + j] + "*, ");
            else System.out.print(this.data[columnSumOffset + j] + ", ");

        }
        System.out.println("");
    }
}
