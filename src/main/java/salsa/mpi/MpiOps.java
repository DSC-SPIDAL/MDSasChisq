package salsa.mpi;

import com.google.common.base.Joiner;
import mpi.Intracomm;
import mpi.MPI;
import mpi.MPIException;
import mpi.Op;
import salsa.mdsaschisq.MPI1DStringVectorPacket;
import salsa.mdsaschisq.MPI2DDoubleVectorPacket;

import java.io.Serializable;

public class MpiOps {
    private int [] intBuff;
    private double [] doubleBuff;
    private boolean [] booleanBuff;
    private String [] stringBuff;
    private MPIReducePlusIndex [] mpiReducePlusIndexBuff;
    private Object [] objectSendBuff, objectRecvBuff;

    private Intracomm comm;
    private int size;

    public MpiOps(Intracomm comm) throws MPIException {
        intBuff = new int[1];
        doubleBuff = new double[1];
        booleanBuff = new boolean[1];
        stringBuff = new String[1];
        mpiReducePlusIndexBuff = new MPIReducePlusIndex[1];
        objectSendBuff = new Object[1];
        objectRecvBuff = new Object[1];
        this.comm = comm;
        size = comm.Size();
    }

    public MpiOps() throws MPIException {
        this(MPI.COMM_WORLD);
    }

    /* AllReduce */
    public int allReduce(int value, Op reduceOp) throws MPIException {
        return allReduce(value, reduceOp, comm);
    }
    public int allReduce(int value, Op reduceOp, Intracomm comm) throws MPIException {
        intBuff[0] = value;
        comm.Allreduce(intBuff,0,intBuff,0,1,MPI.INT,reduceOp);
        return intBuff[0];
    }

    public void allReduce(int [] values, Op reduceOp) throws MPIException{
        allReduce(values, reduceOp, comm);
    }

    public void allReduce(int [] values, Op reduceOp, Intracomm comm) throws MPIException {
        comm.Allreduce(values,0,values,0,values.length,MPI.INT,reduceOp);
    }

    public double allReduce(double value, Op reduceOp) throws MPIException {
        return allReduce(value, reduceOp, comm);
    }
    public double allReduce(double value, Op reduceOp, Intracomm comm) throws MPIException {
        doubleBuff[0] = value;
        comm.Allreduce(doubleBuff,0,doubleBuff,0,1,MPI.DOUBLE,reduceOp);
        return doubleBuff[0];
    }

    public void allReduce(double [] values, Op reduceOp) throws MPIException{
        allReduce(values, reduceOp, comm);
    }

    public void allReduce(double [] values, Op reduceOp, Intracomm comm) throws MPIException {
        comm.Allreduce(values,0,values,0,values.length,MPI.DOUBLE,reduceOp);
    }


    public boolean allReduce(boolean value, Op reduceOp) throws MPIException {
        return allReduce(value, reduceOp, comm);
    }
    public boolean allReduce(boolean value, Op reduceOp, Intracomm comm) throws MPIException {
        booleanBuff[0] = value;
        comm.Allreduce(booleanBuff,0,booleanBuff,0,1,MPI.BOOLEAN,reduceOp);
        return booleanBuff[0];
    }


    public String allReduce(String value) throws MPIException{
        return allReduce(value, comm);
    }

    public String allReduce(String value, Intracomm comm) throws MPIException {
        stringBuff[0] = value;
        String [] result = new String[comm.Size()];
        comm.Allgather(stringBuff,0,1,MPI.OBJECT,result,0,1,MPI.OBJECT);
        return Joiner.on("").join(result);
    }


    public MPIReducePlusIndex allReduce(MPIReducePlusIndex value, MPIReducePlusIndex.Op reduceOp) throws MPIException {
        return allReduce(value, reduceOp, comm);
    }

    public MPIReducePlusIndex allReduce(MPIReducePlusIndex value, MPIReducePlusIndex.Op reduceOp, Intracomm comm) throws MPIException {

        mpiReducePlusIndexBuff[0] = value;
        if (reduceOp == MPIReducePlusIndex.Op.MAX_WITH_INDEX) {
            comm.Allreduce(mpiReducePlusIndexBuff,0,mpiReducePlusIndexBuff,0,1,MPI.OBJECT,MPIReducePlusIndex.maxWithIndex);
        } else if (reduceOp == MPIReducePlusIndex.Op.MIN_WITH_INDEX){
            comm.Allreduce(mpiReducePlusIndexBuff,0,mpiReducePlusIndexBuff,0,1,MPI.OBJECT,MPIReducePlusIndex.minWithIndex);
        }
        return mpiReducePlusIndexBuff[0];
    }

    /* AllGather */
    public int[] allGather(int value) throws MPIException {
        int[] result = new int[size];
        allGather(value, result, comm);
        return result;
    }
    public void allGather(int value, int[] result) throws MPIException {
        allGather(value, result, comm);
    }
    public void allGather(int value, int[] result, Intracomm comm) throws MPIException {
        intBuff[0] = value;
        comm.Allgather(intBuff,0,1,MPI.INT,result,0,1,MPI.INT);
    }

    public double [] allGather (double value) throws MPIException {
        double [] result = new double[size];
        allGather(value, result, comm);
        return result;
    }

    public void allGather(double value, double [] result) throws MPIException {
        allGather(value, result, comm);
    }

    public void allGather(double value, double [] result, Intracomm comm) throws MPIException {
        doubleBuff[0] = value;
        comm.Allgather(doubleBuff,0,1,MPI.DOUBLE,result, 0, 1, MPI.DOUBLE);
    }


    /* Broadcast */
    public int broadcast(int value, int root) throws MPIException{
        return broadcast(value, root, comm);
    }

    public int broadcast(int value, int root, Intracomm comm) throws MPIException {
        intBuff[0] = value;
        comm.Bcast(intBuff,0,1,MPI.INT, root);
        return intBuff[0];
    }

    public void broadcast(int[] values, int root) throws MPIException {
        broadcast(values, root, comm);
    }

    public void broadcast(int[] values, int root, Intracomm comm) throws MPIException {
        comm.Bcast(values, 0, values.length, MPI.INT, root);
    }

    public double broadcast(double value, int root) throws MPIException {
        return broadcast(value, root, comm);
    }

    public double broadcast(double value, int root, Intracomm comm) throws MPIException {
        doubleBuff[0] = value;
        comm.Bcast(doubleBuff,0,1,MPI.DOUBLE,root);
        return doubleBuff[0];
    }

    public void broadcast(double[] values, int root) throws MPIException {
        broadcast(values, root, comm);
    }

    public void broadcast(double[] values, int root, Intracomm comm) throws MPIException {
        comm.Bcast(values, 0, values.length, MPI.DOUBLE, root);
    }

    public boolean broadcast(boolean value, int root) throws MPIException{
        return broadcast(value, root, comm);
    }

    public boolean broadcast(boolean value, int root, Intracomm comm) throws MPIException {
        booleanBuff[0] = value;
        comm.Bcast(booleanBuff,0,1,MPI.BOOLEAN,root);
        return booleanBuff[0];
    }

    public void broadcast (boolean[] values, int root) throws MPIException{
        broadcast(values, root, comm);
    }

    public void broadcast(boolean[] values, int root, Intracomm comm) throws MPIException{
        comm.Bcast(values, 0, values.length,MPI.BOOLEAN,root);
    }

    public MPI1DStringVectorPacket broadcast(MPI1DStringVectorPacket value, int root) throws MPIException{
        return broadcast(value, root, comm);
    }

    public MPI1DStringVectorPacket broadcast(MPI1DStringVectorPacket value, int root, Intracomm comm) throws MPIException{
        objectSendBuff[0] = value;
        comm.Bcast(objectSendBuff,0,1,MPI.OBJECT,root);
        return (MPI1DStringVectorPacket) objectSendBuff[0];
    }

    public MPI2DDoubleVectorPacket broadcast(MPI2DDoubleVectorPacket value, int root) throws MPIException{
        return broadcast(value, root, comm);
    }

    public MPI2DDoubleVectorPacket broadcast(MPI2DDoubleVectorPacket value, int root, Intracomm comm) throws MPIException{
        objectSendBuff[0] = value;
        comm.Bcast(objectSendBuff,0,1,MPI.OBJECT,root);
        return (MPI2DDoubleVectorPacket) objectSendBuff[0];
    }


    /* Sendrecv */
    public <T extends Serializable> T sendReceive(T sendValue, int dest, int destTag, int src, int srcTag){
        return sendReceive(sendValue, dest, destTag, src, srcTag,comm);
    }

    public <T extends Serializable> T sendReceive(T sendValue, int dest, int destTag, int src, int srcTag, Intracomm comm){
        objectSendBuff[0] = sendValue;
        comm.Sendrecv(objectSendBuff, 0, 1, MPI.OBJECT, dest, destTag, objectRecvBuff, 0, 1, MPI.OBJECT, src, srcTag);
        return (T) objectRecvBuff[0]; // unavoidable (without tapping into reflection) unchecked cast  due to how Java generics work
    }


    /* Send */
    public <T extends Serializable> void send(T value, int dest, int tag){
        send(value, dest, tag, comm);
    }

    public <T extends Serializable> void send(T value, int dest, int tag, Intracomm comm){
        objectSendBuff[0] = value;
        comm.Send(objectSendBuff,0,1,MPI.OBJECT,dest,tag);
    }


    /* Receive */
    public <T extends Serializable> T receive(int src, int tag){
        return receive(src, tag, comm);
    }

    public <T extends Serializable> T receive(int src, int tag, Intracomm comm){
        comm.Recv(objectRecvBuff,0,1,MPI.OBJECT,src,tag);
        return (T) objectRecvBuff[0]; // unavoidable (without tapping into reflection) unchecked cast  due to how Java generics work
    }


    /* Barrier */
    public void barrier() throws MPIException {
        barrier(MPI.COMM_WORLD);

    }
    public void barrier(Intracomm comm) throws MPIException {
        comm.Barrier();
    }



}
