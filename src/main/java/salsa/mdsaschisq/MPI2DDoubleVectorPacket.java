package salsa.mdsaschisq;

import mpi.MPI;
import mpi.Struct;

import java.io.Serializable;
import java.nio.ByteBuffer;

/**
 * Simple MPI Serializable packet for single double dimensional array
 */

public class MPI2DDoubleVectorPacket{
    private ByteBuffer buffer;
    private int maxLength;
    private int vectorDimension;
    public MPI2DDoubleVectorPacket(int maxLength, int vectorDimension){
        this.maxLength = maxLength;
        this.vectorDimension = vectorDimension;
    }

    public class MPI2DDoubleVectorPacketType extends Struct{
        private int firstPointOffset = addInt();
        private int numberOfPointsOffset = addInt();

        private int mArrayOffset = addDouble(maxLength * vectorDimension); // mArray represents a 2-dimensional double array flattened into a single dimension


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

            public void loadMarray(double [][] from, Data to){

            }

        }
    }
}