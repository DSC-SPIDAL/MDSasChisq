package salsa.mdsaschisq;

import java.io.Serializable;

public class MPI1DStringVectorPacket implements Serializable {
    public int FirstPoint;
    public int NumberofPoints;
    public String[] Marray;

    public MPI1DStringVectorPacket(int Maxlength, int VectorDimension) {
        Marray = new String[Maxlength];
    }

    public static void CopyMPI1DStringVectorPacket(MPI1DStringVectorPacket ToObject,
                                                   MPI1DStringVectorPacket FromObject) {
        ToObject.FirstPoint = FromObject.FirstPoint;
        ToObject.NumberofPoints = FromObject.NumberofPoints;
        SALSABLAS.CopyVector(ToObject.Marray, FromObject.Marray, 0, FromObject.NumberofPoints);
    }
}