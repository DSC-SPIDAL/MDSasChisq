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
//        buffer.position(mArrayOffset).limit(extent);
        buffer.position(mArrayOffset);
        DoubleBuffer dbuff = buffer.asDoubleBuffer();
        for (double[] aFrom : from) {
            dbuff.put(aFrom);
        }
        buffer.flip();

        // TODO - Debugs - sanity check
        int vecLength = from[0].length;
        double[] vec = new double[vecLength];
        buffer.position(mArrayOffset).limit(extent);
        DoubleBuffer dbuffer = buffer.asDoubleBuffer();
        for (int i = 0; i < numberOfPoints; ++i){
            dbuffer.position(i*vecLength);
            dbuffer.get(vec);
            for (int j = 0; j < vecLength; ++j){
                if (vec[j] != from[i][j]){
                    SALSAUtility.debugPrintCameHere(" in copyToMArray " + i  + "," + j + " elements do not match. from=" + from[i][j] + " to " + vec[j], -2);
                }
            }
        }
    }

    public void copyMArrayTo(double[][] to, int startIndex){
        int numberOfPoints = buffer.getInt(numberOfPointsOffset);
//        buffer.position(mArrayOffset).limit(extent);
        buffer.position(mArrayOffset);
        DoubleBuffer dbuff = buffer.asDoubleBuffer();

        // TODO - Debugs
        if (numberOfPoints+startIndex < 0){
            SALSAUtility.debugPrintCameHere(" +++ in copyMArrayTo numberOfPoints=" + numberOfPoints + " startIdx=" + startIndex, -2);
        }

        for (int i = 0; i < numberOfPoints; i++) {
            dbuff.get(to[i+startIndex],0,vectorDimension);
        }

        // TODO - Debugs - sanity check
        int vecLength = to[0].length;
        double[] vec = new double[vecLength];
        buffer.position(mArrayOffset).limit(extent);
        DoubleBuffer dbuffer = buffer.asDoubleBuffer();
        for (int i = 0; i < numberOfPoints; ++i){
            dbuffer.position(i*vecLength);
            dbuffer.get(vec);
            for (int j = 0; j < vecLength; ++j){
                if (vec[j] != to[i+startIndex][j]){
                    SALSAUtility.debugPrintCameHere(" in copyToMArray " + i  + "," + j + " elements do not match. to[i+startIdx]=" + to[i][j] + " buffer " + vec[j], -2);
                }
            }
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
    }

    public int getFirstPoint(){
        return buffer.getInt(firstPointOffset);
    }
}
