package salsa.mdsaschisq;

import edu.rice.hj.api.SuspendableException;
import mpi.MPI;
import mpi.MPIException;
import mpi.Struct;

import java.io.Serializable;
import java.nio.ByteBuffer;

import static edu.rice.hj.Module0.forallChunked;

/**
 * Simple MPI Serializable packet for single double dimensional array
 */

public class MPI2DDoubleVectorPacket{
    private MPI2DDoubleVectorPacketType.Data mpi2DDoubleVectorPacketTypeData;

    private MPI2DDoubleVectorPacket(int maxLength, int vectorDimension) throws MPIException {
        MPI2DDoubleVectorPacketType mpi2DDoubleVectorPacketType = new MPI2DDoubleVectorPacketType(maxLength, vectorDimension);
        ByteBuffer buffer = MPI.newByteBuffer(mpi2DDoubleVectorPacketType.getExtent());
        mpi2DDoubleVectorPacketTypeData = mpi2DDoubleVectorPacketType.getData(buffer);
    }

    public static MPI2DDoubleVectorPacketType.Data newMPI2DDoubleVectorPacket(int maxLength, int vectorDimension) throws MPIException {
        return new MPI2DDoubleVectorPacket(maxLength, vectorDimension).mpi2DDoubleVectorPacketTypeData;
    }



    public class MPI2DDoubleVectorPacketType extends Struct{
        private int maxLength;
        private int vectorDimension;
        private int firstPointOffset = addInt();
        private int numberOfPointsOffset = addInt();

        private int mArrayOffset = addDouble(maxLength * vectorDimension); // mArray represents a 2-dimensional double array flattened into a single dimension

        public MPI2DDoubleVectorPacketType(int maxLength, int vectorDimension) {
            this.maxLength = maxLength;
            this.vectorDimension = vectorDimension;
        }


        @Override
        protected Data newData() {
            return new Data();
        }

        public class Data extends Struct.Data{
            public int getFirstPoint(){return getInt(firstPointOffset);}
            public void setFirstPoint(int value){putInt(firstPointOffset, value);}

            public int getNumberOfPoints(){return getInt(numberOfPointsOffset);}
            public void setNumberOfPoints(int value){putInt(numberOfPointsOffset, value);}

            public double getMarrayElementAt(int i, int j){
                return getDouble(mArrayOffset, i*vectorDimension+j);
            }
            public void setMarrayElementAt(int i, int j, double value){
                putDouble(mArrayOffset, i*vectorDimension+j, value);
            }

            public void loadMarray(double [][] from, int startIndex, int totalSize){
                if (from == null || from.length == 0) return;

                if (vectorDimension != from[0].length) {
                    SALSAUtility.printAndThrowRuntimeException(
                            "Inconsistent Dimensions - expected " + vectorDimension + " received " + from[0].length);
                }



                if (SALSAUtility.sequentialBLAS) {
                    for (int LongIndex = 0; LongIndex < TotalSize; LongIndex++) {
                        System.arraycopy(from[LongIndex], 0, VectorC[LongIndex + StartIndex], 0, vectorDimension);
                    }
                    return;
                }

                // Parallel VectorC = from
                try {
                    forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
                    {
                        int indexlen = SALSAUtility.PointsperThread[threadIndex];
                        int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                        int endpoint = indexlen + beginpoint;
                        if (endpoint > TotalSize) {
                            endpoint = TotalSize;
                        }
                        if (threadIndex == (SALSAUtility.ThreadCount - 1)) {
                            endpoint = TotalSize;
                        }
                        for (int LongIndex = beginpoint; LongIndex < endpoint; LongIndex++) {
                            System.arraycopy(from[LongIndex], 0, VectorC[LongIndex + StartIndex], 0, vectorDimension);
                        }
                    });
                } catch (SuspendableException e) {
                    SALSAUtility.printAndThrowRuntimeException(e.getMessage());
                }
            }

        }
    }
}