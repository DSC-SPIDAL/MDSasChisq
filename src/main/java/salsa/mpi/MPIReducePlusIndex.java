package salsa.mpi;

import mpi.Datatype;
import mpi.MPI;
import mpi.MPIException;
import mpi.User_Function;

import java.io.Serializable;

public class MPIReducePlusIndex implements Serializable
{
    private int index;
    private double value;

    public MPIReducePlusIndex(int index, double value)
    {
        this.index = index;
        this.value = value;
    }

    public int getIndex() {
        return index;
    }

    public double getValue() {
        return value;
    }

    private static final String MPI_USER_FUNCTION_MIN_WITH_INDEX_INVALID_DATA_TYPE = "MPI User Function - MinWithIndex: Invalid data type";
    private static final String MPI_USER_FUNCTION_MAX_WITH_INDEX_INVALID_DATA_TYPE = "MPI User Function - MaxWithIndex: Invalid data type";

    // Return Value and Index corresponding to minimum over all MPI Processes
    static mpi.Op minWithIndex = new mpi.Op(new User_Function() {
        @Override
        public void Call(Object inVec, int inOffset, Object outVec, int outOffset, int count, Datatype datatype) throws
                MPIException {
            if (datatype == MPI.OBJECT){
                Object [] inArr = (Object[]) inVec;
                Object [] outArr = (Object[]) outVec;
                int inDisplacement  = inOffset ;
                int outDisplacement = outOffset ;
                for (int i = 0; i < count; ++i, ++inDisplacement, ++outDisplacement){
                    MPIReducePlusIndex inValue = (MPIReducePlusIndex)inArr[inDisplacement];
                    MPIReducePlusIndex outValue = (MPIReducePlusIndex)outArr[outDisplacement];
                    if (outValue.getIndex() < 0 || inValue.getValue() < outValue.getValue())
                    {
                        outArr[outDisplacement] = inValue;
                    }
                }
            } else {
                System.out.println(
                        MPI_USER_FUNCTION_MIN_WITH_INDEX_INVALID_DATA_TYPE);
                try {
                    MPI.COMM_WORLD.Abort(1);
                } catch (MPIException ignored){}
            }

        }
    }, true);

    // Return Value and Index corresponding to maximum over all MPI Processes
    static mpi.Op maxWithIndex = new mpi.Op(new User_Function() {
        @Override
        public void Call(Object inVec, int inOffset, Object outVec, int outOffset, int count, Datatype datatype) throws MPIException {
            if (datatype == MPI.OBJECT){
                Object [] inArr = (Object[]) inVec;
                Object [] outArr = (Object[]) outVec;
                int inDisplacement  = inOffset ;
                int outDisplacement = outOffset ;
                for (int i = 0; i < count; ++i, ++inDisplacement, ++outDisplacement){
                    MPIReducePlusIndex inValue = (MPIReducePlusIndex)inArr[inDisplacement];
                    MPIReducePlusIndex outValue = (MPIReducePlusIndex)outArr[outDisplacement];
                    if (outValue.getIndex() < 0 || inValue.getValue() > outValue.getValue())
                    {
                        outArr[outDisplacement] = inValue;
                    }
                }
            } else {
                System.out.println(
                        MPI_USER_FUNCTION_MAX_WITH_INDEX_INVALID_DATA_TYPE);
                try {
                    MPI.COMM_WORLD.Abort(1);
                } catch (MPIException ignored){}
            }

        }
    }, true);

    public enum  Op{
        MIN_WITH_INDEX, MAX_WITH_INDEX
    }
}

