package salsa.mdsaschisq;

import mpi.MPI;

import java.io.Serializable;
import java.nio.ByteBuffer;

/**
 * Simple MPI Serializable packet for single double dimensional array
 */
public class MPI2DDoubleVectorPacket implements Serializable {
    private final int shift = Double.BYTES;
    private final int extent;
    private final int mArrayDim0Length; // length of dimension 0
    private final int mArrayDim1Length; // length of dimension 1
    ByteBuffer buffer; // data buffer
    private static int firstPointOffset = 0;
    private static int numberOfPointsOffset = Integer.BYTES;
    private static int mArrayDim0LengthOffset = 2*Integer.BYTES;
    private static int mArrayDim1LengthOffset = 3*Integer.BYTES;
    private static int mArrayOffset = 4*Integer.BYTES;

    public MPI2DDoubleVectorPacket(int maxLength, int vectorDimension) {
        mArrayDim0Length = maxLength;
        mArrayDim1Length = vectorDimension;
        extent = mArrayDim0Length * mArrayDim1Length * Double.BYTES + 4*Integer.BYTES;
        buffer = MPI.newByteBuffer(extent);
        buffer.putInt(mArrayDim0LengthOffset, mArrayDim0Length);
        buffer.putInt(mArrayDim1LengthOffset, mArrayDim1Length);
    }

    public int getMArrayDim0Length (){
        return mArrayDim0Length;
    }

    public void setMArrayDim0Length (int length){
        buffer.putInt(mArrayDim0LengthOffset, length);
    }

    public int getMArrayDim1Length (){
        return mArrayDim1Length;
    }

    public void setMArrayDim1Length (int length){
        buffer.putInt(mArrayDim1LengthOffset, length);
    }

    public int getFirstPoint(){
        return buffer.getInt(firstPointOffset);
    }

    public void setFirstPoint(int value){
        buffer.putInt(firstPointOffset, value);
    }

    public int getNumberOfPoints(){
        return buffer.getInt(numberOfPointsOffset);
    }

    public void setNumberOfPoints(int value){
        buffer.putInt(numberOfPointsOffset, value);
    }





    public static void CopyMPI2DDoubleVectorPacket(MPI2DDoubleVectorPacket ToObject,
                                                   MPI2DDoubleVectorPacket FromObject) {
        ToObject.FirstPoint = FromObject.FirstPoint;
        ToObject.NumberofPoints = FromObject.NumberofPoints;
        SALSABLAS.CopyVector(ToObject.Marray, FromObject.Marray, 0, FromObject.NumberofPoints);

    }
}