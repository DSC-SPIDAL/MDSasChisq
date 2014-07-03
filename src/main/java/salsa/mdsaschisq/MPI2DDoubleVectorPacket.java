package salsa.mdsaschisq;

import java.io.Serializable;

/**
 * Simple MPI Serializable packet for single double dimensional array
 */
public class MPI2DDoubleVectorPacket implements Serializable {
    public int FirstPoint;
    public int NumberofPoints;
    public double[][] Marray;

    public MPI2DDoubleVectorPacket(int Maxlength, int VectorDimension) {
        Marray = new double[Maxlength][];
        for (int LongIndex = 0; LongIndex < Maxlength; LongIndex++) {
            Marray[LongIndex] = new double[VectorDimension];
        }
    }

    public static void CopyMPI2DDoubleVectorPacket(MPI2DDoubleVectorPacket ToObject,
                                                   MPI2DDoubleVectorPacket FromObject) {
        ToObject.FirstPoint = FromObject.FirstPoint;
        ToObject.NumberofPoints = FromObject.NumberofPoints;
        SALSABLAS.CopyVector(ToObject.Marray, FromObject.Marray, 0, FromObject.NumberofPoints);

    }
}