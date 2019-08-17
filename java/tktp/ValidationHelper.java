package tktp;

import java.util.stream.*;
import java.util.List;
import java.util.Arrays;

public final class ValidationHelper {
  public static void main(String[] a) {

    if (a.length != 3) {
      System.out.println("Please provide an implementation and two equal-length vectors as input:\ne.g. java -cp tktp.jar tptp.ValidationHelper YU 1,2,3,4 4,3,2,1");
      System.exit(1);
    }

    double[] x = Stream.of(a[1].split(","))
      .mapToDouble (e -> new Double(Double.parseDouble(e)))
      .toArray();

    double[] y = Stream.of(a[2].split(","))
      .mapToDouble (e -> new Double(Double.parseDouble(e)))
      .toArray();

    if (x.length != y.length) {
      System.out.println("Please provide equal-length vectors.");
      System.exit(1);
    }

    int[] pi = null;
    switch (a[0]) {
      case "FastBCS2": pi = FastBCS2.getPi(x, y).pi();
      case "FastBCS": pi = FastBCS.getPi(x, y);
    }

    System.out.println(Arrays.toString(pi).replace("[", "").replace("]", ""));
  }
}
