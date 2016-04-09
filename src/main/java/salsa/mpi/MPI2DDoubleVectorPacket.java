package salsa.mpi;

import mpi.MPI;
import mpi.MPIException;
import salsa.mdsaschisq.SALSAUtility;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;

/**
 * Simple MPI Serializable packet for single double dimensional array
 */

public class MPI2DDoubleVectorPacket{
    private static final int firstPointOffset = 0;
    private static final int numberOfPointsOffset = Integer.BYTES;
    private static final int mArrayOffset = 2*Integer.BYTES;

    ByteBuffer buffer;
    int extent;
    private int vectorDimension;

    public MPI2DDoubleVectorPacket(int firstPoint, int numberOfPoints, int maxLength, int vectorDimension) throws MPIException {
        this.vectorDimension = vectorDimension;
        extent = (maxLength*vectorDimension)*Double.BYTES + 2*Integer.BYTES; // 2 is to add two integers to represent first point and number of points
        buffer = MPI.newByteBuffer(extent);
        buffer.limit(extent);
        buffer.putInt(firstPointOffset,firstPoint).putInt(numberOfPointsOffset, numberOfPoints);

    }

    public void copyToMArray(double[][] from, int numberOfPoints){
        if (from == null || numberOfPoints != buffer.getInt(numberOfPointsOffset) || from[0].length != vectorDimension){
            SALSAUtility.printAndThrowRuntimeException("Error while copying double[][] to mArray");
            return;
        }
        buffer.position(mArrayOffset).limit(extent);
        DoubleBuffer dbuff = buffer.asDoubleBuffer();
        for (double[] aFrom : from) {
            dbuff.put(aFrom);
        }
        buffer.flip();
    }

    public void copyMArrayTo(double[][] to, int startIndex){
        try {
            int rank = MPI.COMM_WORLD.getRank();
            // TODO - Debugs
            System.out.println(
                "Rank=" + rank + " to=[" + to.length + "][" + to[0].length
                + "]" + " startIdx=" + startIndex + " numberOfPoints=" + buffer.getInt(numberOfPointsOffset));
        }
        catch (MPIException e) {
            e.printStackTrace();
        }


        int numberOfPoints = buffer.getInt(numberOfPointsOffset);
        buffer.position(mArrayOffset).limit(extent);
        DoubleBuffer dbuff = buffer.asDoubleBuffer();
        for (int i = 0; i < numberOfPoints; i++) {
            dbuff.get(to[i+startIndex],0,vectorDimension);
        }
    }

    public void setMArrayElementAt(int i, int j, double value){
        buffer.putDouble(i*vectorDimension+j, value);
    }

    public static void copy(MPI2DDoubleVectorPacket from, MPI2DDoubleVectorPacket to){
        from.buffer.position(0);
        to.buffer.position(0);
        to.buffer.put(from.buffer);
        to.buffer.flip();
        from.buffer.flip();
        /*int capacity = from.buffer.capacity();
        from.buffer.position(0).limit(capacity);
        to.buffer.position(0).limit(capacity);
        to.buffer.put(from.buffer);
        to.buffer.flip();
        from.buffer.flip();*/
    }

    public int getFirstPoint(){
        return buffer.getInt(firstPointOffset);
    }
}
