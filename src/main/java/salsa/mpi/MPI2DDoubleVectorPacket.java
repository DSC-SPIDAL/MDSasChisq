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
    private int extent;
    private int vectorDimension;

    public MPI2DDoubleVectorPacket(int firstPoint, int numberOfPoints, int maxLength, int vectorDimension) throws MPIException {
        this.vectorDimension = vectorDimension;
        extent = (maxLength*vectorDimension)*Double.BYTES + 2*Integer.BYTES; // 2 is to add two integers to represent first point and number of points
        buffer = MPI.newByteBuffer(extent);
        buffer.limit(extent);
        buffer.putInt(firstPointOffset,firstPoint).putInt(numberOfPointsOffset, numberOfPoints);

    }

    public void copyToMArray(double[][] from){
        if (from == null || from.length != buffer.getInt(numberOfPointsOffset) || from[0].length != vectorDimension){
            SALSAUtility.printAndThrowRuntimeException("Error while copying double[][] to mArray");
        }
        buffer.position(mArrayOffset).limit(extent);
        DoubleBuffer dbuff = buffer.asDoubleBuffer();
        // TODO - array to buffer copy - see if it's faster if HJ is used to copy elements individually in parallel than the sequential bulk method used here
        for (double[] aFrom : from) {
            dbuff.put(aFrom);
        }
        buffer.flip();
    }

    public double[][] getMArray(){
        int numberOfPoints = buffer.getInt(numberOfPointsOffset);
        double [][] array = new double[numberOfPoints][vectorDimension];
        buffer.position(mArrayOffset).limit(extent);
        DoubleBuffer dbuff = buffer.asDoubleBuffer();
        for (int i = 0; i < numberOfPoints; i++) {
           dbuff.get(array[i],0,vectorDimension);
        }
        return array;
    }

    public static void copy(MPI2DDoubleVectorPacket from, MPI2DDoubleVectorPacket to){
        int capacity = from.buffer.capacity();
        from.buffer.position(0).limit(capacity);
        to.buffer.position(0).limit(capacity);
        to.buffer.put(from.buffer);
        to.buffer.flip();
        from.buffer.flip();
    }




}
