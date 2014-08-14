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

    public class VectorType extends Struct{
        int vectorOffset = addDou
        @Override
        protected Data newData() {
            return new Data()
        }

        public class Data extends Struct.Data{

        }
    }

    public class MPI2DDoubleVectorPacketType extends Struct{
        private int firstPointOffset = addInt();
        private int numberOfPointsOffset = addInt();

        private int mArrayOffset; // mArray represents a 2-dimensional


        @Override
        protected Data newData() {
            return new Data();
        }

        public class Data extends Struct.Data{

        }
    }
}