package salsa.mdsaschisq;

import java.io.Serializable;

public class ChisqFirstandSecond implements Serializable {
    public double chisq;
    public double[] first;
    public double[] second;
    public int VectorSize;
    public int SecondSize;

    public ChisqFirstandSecond(int NumberParms) {
        VectorSize = NumberParms;
        first = new double[NumberParms];
        SecondSize = (NumberParms * (1 + NumberParms)) / 2;
        second = new double[SecondSize];
        chisq = 0.0;
        for (int iparm = 0; iparm < NumberParms; iparm++) {
            first[iparm] = 0.0;
        }
        for (int iparm = 0; iparm < SecondSize; iparm++) {
            second[iparm] = 0.0;
        }
    }

    public static void initializeToZero(ChisqFirstandSecond instance) {
        instance.chisq = 0.0;
        for (int iparm = 0; iparm < instance.VectorSize; iparm++) {
            instance.first[iparm] = 0.0;
        }
        for (int iparm = 0; iparm < instance.SecondSize; iparm++) {
            instance.second[iparm] = 0.0;
        }
    }
}