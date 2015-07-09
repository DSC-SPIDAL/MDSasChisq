package salsa.mdsaschisq;

import edu.rice.hj.api.SuspendableException;
import mpi.MPI;
import mpi.MPIException;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;

public class GlobalReductions {

    public static class FindBoolOr {
        public double[] NumberofPoints;
        public int NumberofThreads;
        public boolean[] Orvalue;
        public double TotalNumberofPoints;
        public boolean TotalOr;

        public FindBoolOr(int NumThreads) {
            NumberofThreads = NumThreads;

            NumberofPoints = new double[NumThreads];
            Orvalue = new boolean[NumThreads];

            TotalNumberofPoints = 0.0;
            TotalOr = false;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                Orvalue[ThreadNo] = false;
            }
        }

        public final void addAPoint(int ThreadNo, boolean value) {
            NumberofPoints[ThreadNo] += 1.0;
            Orvalue[ThreadNo] = Orvalue[ThreadNo] || value;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                TotalOr = Orvalue[ThreadNo] || TotalOr;
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - boolean - or
                TotalOr = SALSAUtility.mpiOps.allReduce(TotalOr, MPI.LOR);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
        }
    }

    public static class FindIntSum {
        private int[] NumberofPoints;
        private int NumberofThreads;
        private int[] Intvalue;
        public int TotalNumberofPoints;
        public int TotalInt;

        public FindIntSum(int NumThreads) {
            NumberofThreads = NumThreads;

            NumberofPoints = new int[NumThreads];
            Intvalue = new int[NumThreads];

            TotalNumberofPoints = 0;
            TotalInt = 0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0;
                Intvalue[ThreadNo] = 0;
            }
        }

        public final void addAPoint(int ThreadNo, int value) {
            NumberofPoints[ThreadNo] += 1;
            Intvalue[ThreadNo] += value;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                TotalInt += Intvalue[ThreadNo];
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - int - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - int - sum
                TotalInt = SALSAUtility.mpiOps.allReduce(TotalInt, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
        }
    }

    /**
     * Used to do histograms
     * Must call startThread method at start of threads
     */
    public static class FindDoubleArraySum {
        public double[] NumberofPoints;
        public int NumberofThreads;
        public int NumberinSum;
        public double[][] Sum;
        public double TotalNumberofPoints;
        public double[] TotalSum;

        public FindDoubleArraySum(int NumThreads, int ArraySize) {
            NumberofThreads = NumThreads;
            NumberinSum = ArraySize;

            NumberofPoints = new double[NumThreads];
            TotalSum = new double[ArraySize];
            Sum = new double[NumThreads][];

            TotalNumberofPoints = 0.0;
            for (int loop = 0; loop < ArraySize; loop++) {
                TotalSum[loop] = 0.0;
            }
        }

        public final void startThread(int ThreadNo) {
            NumberofPoints[ThreadNo] = 0.0;
            Sum[ThreadNo] = new double[NumberinSum];
            for (int loop = 0; loop < NumberinSum; loop++) {
                Sum[ThreadNo][loop] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, int loopvalue) {
            if ((loopvalue < 0) || (loopvalue >= NumberinSum)) {
                return;
            }
            NumberofPoints[ThreadNo] += 1.0;
            Sum[ThreadNo][loopvalue] += 1.0;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                for (int loop = 0; loop < NumberinSum; loop++) {
                    TotalSum[loop] += Sum[ThreadNo][loop];
                }
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - int - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
                SALSAUtility.mpiOps.allReduce(TotalSum, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
        }
    }

    public static class FindDoubleMax {
        public double[] NumberofPoints;
        public int NumberofThreads;
        public double[] Maxvalue;
        public double TotalNumberofPoints;
        public double TotalMax;

        public FindDoubleMax(int NumThreads) {
            NumberofThreads = NumThreads;

            NumberofPoints = new double[NumThreads];
            Maxvalue = new double[NumThreads];

            TotalNumberofPoints = 0.0;
            TotalMax = 0.0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                Maxvalue[ThreadNo] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, double value) {
            NumberofPoints[ThreadNo] += 1.0;
            Maxvalue[ThreadNo] = Math.max(Maxvalue[ThreadNo], value);
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                TotalMax = Math.max(TotalMax, Maxvalue[ThreadNo]);
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double - max
                TotalMax = SALSAUtility.mpiOps.allReduce(TotalMax, MPI.MAX);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
        }
    }

    public static class FindMeanSigma {
        public double[] NumberofPoints;
        public int NumberofThreads;
        public double[] mean;
        public double[] square;
        public double TotalNumberofPoints;
        public double Totalmean;
        public double Totalsquare;
        public double Totalsigma;

        public FindMeanSigma(int NumThreads) {
            NumberofThreads = NumThreads;

            NumberofPoints = new double[NumThreads];
            mean = new double[NumThreads];
            square = new double[NumThreads];

            TotalNumberofPoints = 0.0;
            Totalmean = 0.0;
            Totalsquare = 0.0;
            Totalsigma = 0.0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                mean[ThreadNo] = 0.0;
                square[ThreadNo] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, double value1) {
            NumberofPoints[ThreadNo] += 1.0;
            mean[ThreadNo] += value1;
            square[ThreadNo] += value1 * value1;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                Totalmean += mean[ThreadNo];
                Totalsquare += square[ThreadNo];
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Totalmean = SALSAUtility.mpiOps.allReduce(Totalmean, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Totalsquare = SALSAUtility.mpiOps.allReduce(Totalsquare, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

            if (TotalNumberofPoints < 0.5) {
                return;
            }

            Totalmean = Totalmean / TotalNumberofPoints;
            Totalsquare = (Totalsquare / TotalNumberofPoints) - Totalmean * Totalmean;
            Totalsigma = Math.sqrt(Math.max(0.0, Totalsquare));
        }
    }

    public static class FindDoubleSum {
        private double[] NumberofPoints;
        private int NumberofThreads;
        private double[] TotalinThread;
        public double TotalNumberofPoints;
        public double Total;

        public FindDoubleSum(int NumThreads) {
            NumberofThreads = NumThreads;

            NumberofPoints = new double[NumThreads];
            TotalinThread = new double[NumThreads];

            TotalNumberofPoints = 0.0;
            Total = 0.0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                TotalinThread[ThreadNo] = 0.0;
            }
        }

        public final void zero() {
            this.TotalNumberofPoints = 0.0;
            this.Total = 0.0;
            for (int ThreadNo = 0; ThreadNo < this.NumberofThreads; ThreadNo++) {
                this.NumberofPoints[ThreadNo] = 0.0;
                this.TotalinThread[ThreadNo] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, double value1) {
            NumberofPoints[ThreadNo] += 1.0;
            TotalinThread[ThreadNo] += value1;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                Total += TotalinThread[ThreadNo];
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Total = SALSAUtility.mpiOps.allReduce(Total, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

        }
    }

    public static class FindDoubleMean {
        public double[] NumberofPoints;
        public int NumberofThreads;
        public double[] mean;
        public double TotalNumberofPoints;
        public double Totalmean;

        public FindDoubleMean(int NumThreads) {
            NumberofThreads = NumThreads;

            NumberofPoints = new double[NumThreads];
            mean = new double[NumThreads];

            TotalNumberofPoints = 0.0;
            Totalmean = 0.0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                mean[ThreadNo] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, double value1) {
            NumberofPoints[ThreadNo] += 1.0;
            mean[ThreadNo] += value1;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                Totalmean += mean[ThreadNo];
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Totalmean = SALSAUtility.mpiOps.allReduce(Totalmean, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

            if (TotalNumberofPoints < 0.5) {
                return;
            }

            Totalmean = Totalmean / TotalNumberofPoints;
        }
    }

    public static class FindVectorIntSum {
        private int[] NumberofPoints;
        private int NumberofThreads;
        private int[][] VectorSum;
        public int TotalNumberofPoints;
        public int[] TotalVectorSum;
        private int ArraySize;

        public FindVectorIntSum(int NumThreads, int NumberinArray) {
            NumberofThreads = NumThreads;
            ArraySize = NumberinArray;

            NumberofPoints = new int[NumThreads];
            TotalVectorSum = new int[ArraySize];
            VectorSum = new int[NumThreads][];

            TotalNumberofPoints = 0;
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                TotalVectorSum[ArrayLoop] = 0;
            }
        }

        public final void startThread(int ThreadNo) {
            NumberofPoints[ThreadNo] = 0;
            VectorSum[ThreadNo] = new int[ArraySize];
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                VectorSum[ThreadNo][ArrayLoop] = 0;
            }
        }

        public final void addAPoint(int ThreadNo, int[] value1) {
            NumberofPoints[ThreadNo] += 1;
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
            }
        }

        public final void addAPoint(int ThreadNo, int value1, int position) {
            NumberofPoints[ThreadNo] += 1;
            VectorSum[ThreadNo][position] += value1;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                    TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                }
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - int - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - int[] - sum
                SALSAUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
        }

    }

    public static class FindVectorDoubleMax {
        private double[] NumberofPoints;
        private int NumberofThreads;
        private double[][] VectorMax;
        public double TotalNumberofPoints;
        public double[] TotalVectorMax;
        private int ArraySize;

        public FindVectorDoubleMax(int NumThreads, int NumberinArray) {
            NumberofThreads = NumThreads;
            ArraySize = NumberinArray;

            NumberofPoints = new double[NumThreads];
            TotalVectorMax = new double[ArraySize];
            VectorMax = new double[NumThreads][];

            TotalNumberofPoints = 0.0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                VectorMax[ThreadNo] = new double[ArraySize];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                    VectorMax[ThreadNo][ArrayLoop] = -1.0E10;
                }
            }
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                TotalVectorMax[ArrayLoop] = -1.0E10;
            }
        }

        public final void addAPoint(int ThreadNo, double[] value) {
            NumberofPoints[ThreadNo] += 1.0;
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                VectorMax[ThreadNo][ArrayLoop] = Math.max(VectorMax[ThreadNo][ArrayLoop], value[ArrayLoop]);
            }
        }

        public final void addAPoint(int ThreadNo, double value1, int position) {
            NumberofPoints[ThreadNo] += 1.0;
            VectorMax[ThreadNo][position] = Math.max(VectorMax[ThreadNo][position], value1);
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                    TotalVectorMax[ArrayLoop] = Math.max(TotalVectorMax[ArrayLoop], VectorMax[ThreadNo][ArrayLoop]);
                }
            }

            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - max
                SALSAUtility.mpiOps.allReduce(TotalVectorMax, MPI.MAX);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

        }

    }

    public static class FindVectorDoubleSum {
        private double[] NumberofPoints;
        private int NumberofThreads;
        private double[][] VectorSum;
        public double TotalNumberofPoints;
        public double[] TotalVectorSum;
        private int ArraySize;

        public FindVectorDoubleSum(int NumThreads, int NumberinArray) {
            NumberofThreads = NumThreads;
            ArraySize = NumberinArray;

            NumberofPoints = new double[NumThreads];
            TotalVectorSum = new double[ArraySize];
            VectorSum = new double[NumThreads][];

            TotalNumberofPoints = 0.0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                VectorSum[ThreadNo] = new double[ArraySize];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                    VectorSum[ThreadNo][ArrayLoop] = 0.0;
                }
            }
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                TotalVectorSum[ArrayLoop] = 0.0;
            }
        }

        public final void zero() {
            this.TotalNumberofPoints = 0.0;
            for (int ThreadNo = 0; ThreadNo < this.NumberofThreads; ThreadNo++) {
                this.NumberofPoints[ThreadNo] = 0.0;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                    this.VectorSum[ThreadNo][ArrayLoop] = 0.0;
                }
            }
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                this.TotalVectorSum[ArrayLoop] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, double[] value1) {
            NumberofPoints[ThreadNo] += 1.0;
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
            }
        }

        public final void addAPoint(int ThreadNo, double value1, int position) {
            NumberofPoints[ThreadNo] += 1.0;
            VectorSum[ThreadNo][position] += value1;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                    TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                }
            }

            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
                SALSAUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

        }

    }

    public static class FindVectorDoubleSum2 {
        private double[] NumberofPoints;
        private int NumberofThreads;
        private double[][] VectorSum;
        public double TotalNumberofPoints;
        public double[] TotalVectorSum;
        private int ArraySize;
        private Range[] ParallelArrayRanges;

        public FindVectorDoubleSum2(int NumThreads, int NumberinArray) {
            NumberofThreads = NumThreads;
            ArraySize = NumberinArray;

            NumberofPoints = new double[NumThreads];
            TotalVectorSum = new double[ArraySize];
            VectorSum = new double[NumThreads][];

            TotalNumberofPoints = 0.0;
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                TotalVectorSum[ArrayLoop] = 0.0;
            }

            ParallelArrayRanges = RangePartitioner.Partition(NumberinArray, SALSAUtility.ThreadCount);
        }

        public final void startThread(int ThreadNo) {
            NumberofPoints[ThreadNo] = 0.0;
            VectorSum[ThreadNo] = new double[ArraySize];
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                VectorSum[ThreadNo][ArrayLoop] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, double[] value1) {
            NumberofPoints[ThreadNo] += 1.0;
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
            }
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            // Note - parallel for

            launchHabaneroApp(() ->forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) -> {
                int beginindex = ParallelArrayRanges[threadIndex].getStartIndex();
                int indexlength = ParallelArrayRanges[threadIndex].getLength();
                for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++) {
                    TotalVectorSum[ArrayLoop] = 0.0;
                    for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                        TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                    }
                }
            }));

            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
                SALSAUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

        }

    }

    public static class FindVectorDoubleSum3 {
        private double[] NumberofPoints;
        private int NumberofThreads;
        private double[][] VectorSum;
        public double TotalNumberofPoints;
        public double[] TotalVectorSum;
        private int ArraySize;
        private int ArraySize1;
        private int ArraySize2;
        private Range[] ParallelArrayRanges;

        public FindVectorDoubleSum3(int NumThreads, int NumberinArray1, int NumberinArray2) {
            NumberofThreads = NumThreads;
            ArraySize = NumberinArray1 * NumberinArray2;
            ArraySize1 = NumberinArray1;
            ArraySize2 = NumberinArray2;

            NumberofPoints = new double[NumThreads];
            TotalVectorSum = new double[ArraySize];

            VectorSum = new double[NumThreads][];

            TotalNumberofPoints = 0.0;

            ParallelArrayRanges = RangePartitioner.Partition(ArraySize, SALSAUtility.ThreadCount);
        }

        public final void startThread(int ThreadNo) {
            NumberofPoints[ThreadNo] = 0.0;
            VectorSum[ThreadNo] = new double[ArraySize];
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                VectorSum[ThreadNo][ArrayLoop] = 0.0;
            }

        }

        public final void addAPoint(int ThreadNo, double[] value1, int index2) {
            NumberofPoints[ThreadNo] += 1.0;
            int additive = index2 * ArraySize1;
            for (int ArrayLoop1 = 0; ArrayLoop1 < ArraySize1; ArrayLoop1++) {
                VectorSum[ThreadNo][ArrayLoop1 + additive] += value1[ArrayLoop1];
            }
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            SALSAUtility.StartSubTimer(SALSAUtility.ThreadTiming);
            // Note - parallel for

            launchHabaneroApp(() ->forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) -> {
                int beginindex = ParallelArrayRanges[threadIndex].getStartIndex();
                int indexlength = ParallelArrayRanges[threadIndex].getLength();
                for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++) {
                    double tmp = 0.0;
                    for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                        tmp += VectorSum[ThreadNo][ArrayLoop];
                    }
                    TotalVectorSum[ArrayLoop] = tmp;
                }
            }));

            SALSAUtility.StopSubTimer(SALSAUtility.ThreadTiming);

            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                int bigsize = TotalVectorSum.length;
                if (bigsize <= 4096) {
                    // Note - MPI Call - Allreduce - double[] - sum
                    SALSAUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                } else {
                    double[] buffer = new double[4096];
                    int start = 0;
                    while (start < bigsize) {
                        int whatsLeft = Math.min(bigsize - start, 4096);
                        System.arraycopy(TotalVectorSum, start, buffer, 0, whatsLeft);
                        // Note - MPI Call - Allreduce - double[] - sum
                        SALSAUtility.mpiOps.allReduce(buffer, MPI.SUM);
                        System.arraycopy(buffer, 0, TotalVectorSum, start, whatsLeft);
                        start += whatsLeft;
                    }

                }
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
        }

    }

    public static class FindIndirectVectorDoubleSum {
        private double[] NumberofPoints;
        private int NumberofThreads;
        private double[][] VectorSum;
        public double TotalNumberofPoints;
        public double[] TotalVectorSum;
        private int ArraySize;
        private Range[] ParallelArrayRanges;

        public FindIndirectVectorDoubleSum(int NumThreads, int NumberinArray) {
            NumberofThreads = NumThreads;
            ArraySize = NumberinArray;

            NumberofPoints = new double[NumThreads];
            TotalVectorSum = new double[ArraySize];
            VectorSum = new double[NumThreads][];

            TotalNumberofPoints = 0.0;

            ParallelArrayRanges = RangePartitioner.Partition(NumberinArray, SALSAUtility.ThreadCount);
        }

        public final void startThread(int ThreadNo) {
            NumberofPoints[ThreadNo] = 0.0;
            VectorSum[ThreadNo] = new double[ArraySize];
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                VectorSum[ThreadNo][ArrayLoop] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, int NumberLocations, int[] location, double[] value1) {
            if (NumberLocations <= 0) {
                return;
            }
            NumberofPoints[ThreadNo] += 1.0;
            for (int ArrayLoop = 0; ArrayLoop < NumberLocations; ArrayLoop++) {
                VectorSum[ThreadNo][location[ArrayLoop]] += value1[ArrayLoop];
            }
        }

        public final void sumOverThreadsAndMPI() throws MPIException {

            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
            }

            // Note - parallel for

            launchHabaneroApp(() ->forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) -> {
                int beginindex = ParallelArrayRanges[threadIndex].getStartIndex();
                int indexlength = ParallelArrayRanges[threadIndex].getLength();
                for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++) {
                    TotalVectorSum[ArrayLoop] = 0.0;
                    for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                        TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                    }
                }
            }));


            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
                SALSAUtility.mpiOps.allReduce(TotalVectorSum, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

        }

    }

    public static class FindArrayMean {
        public double[] NumberofPoints;
        public int NumberofThreads;
        public double[][] mean;
        public double TotalNumberofPoints;
        public double[] Totalmean;
        public int ArraySize;

        public FindArrayMean(int NumThreads, int NumberinArray) {
            NumberofThreads = NumThreads;
            ArraySize = NumberinArray;

            NumberofPoints = new double[NumThreads];
            Totalmean = new double[ArraySize];
            mean = new double[NumThreads][];

            TotalNumberofPoints = 0.0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                mean[ThreadNo] = new double[ArraySize];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                    mean[ThreadNo][ArrayLoop] = 0.0;
                }
            }
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                Totalmean[ArrayLoop] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, double[] value1) {
            NumberofPoints[ThreadNo] += 1.0;
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                mean[ThreadNo][ArrayLoop] += value1[ArrayLoop];
            }
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                    Totalmean[ArrayLoop] += mean[ThreadNo][ArrayLoop];
                }
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double[] - sum
                SALSAUtility.mpiOps.allReduce(Totalmean, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

            if (TotalNumberofPoints < 0.5) {
                return;
            }
            for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++) {
                Totalmean[ArrayLoop] = Totalmean[ArrayLoop] / TotalNumberofPoints;
            }
        }
    }

    public static class FindMinorMaxValuewithIndex {
        public double[] NumberofPoints;
        public int NumberofThreads;
        public double[] MaxOrMinvalue;
        public double TotalNumberofPoints;
        public double TotalMaxOrMin;
        public int TotalIndexValue;
        public int[] IndexValue;
        public int MinMaxPointer; // =0 Min = 1 Max

        public FindMinorMaxValuewithIndex(int NumThreads, int UseMax) {
            NumberofThreads = NumThreads;
            MinMaxPointer = UseMax;

            NumberofPoints = new double[NumThreads];
            MaxOrMinvalue = new double[NumThreads];
            IndexValue = new int[NumThreads];

            TotalNumberofPoints = 0.0;
            TotalMaxOrMin = 0.0;
            TotalIndexValue = -1;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                MaxOrMinvalue[ThreadNo] = 0.0;
                IndexValue[ThreadNo] = -1;
            }
        }

        public final void addAPoint(int ThreadNo, int indexposition, double value) {
            NumberofPoints[ThreadNo] += 1.0;
            if (MinMaxPointer != 0) { // Max
                if ((IndexValue[ThreadNo] >= 0) && (MaxOrMinvalue[ThreadNo] > value)) {
                    return;
                }
            } else { // Min
                if ((IndexValue[ThreadNo] >= 0) && (MaxOrMinvalue[ThreadNo] <= value)) {
                    return;
                }
            }
            MaxOrMinvalue[ThreadNo] = value;
            IndexValue[ThreadNo] = indexposition;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                if (IndexValue[ThreadNo] < 0) {
                    continue;
                }

                TotalNumberofPoints += NumberofPoints[ThreadNo];
                if (MinMaxPointer != 0) {
                    if ((TotalIndexValue >= 0) && (TotalMaxOrMin > MaxOrMinvalue[ThreadNo])) {
                        continue;
                    }
                } else {
                    if ((TotalIndexValue >= 0) && (TotalMaxOrMin <= MaxOrMinvalue[ThreadNo])) {
                        continue;
                    }
                }

                TotalMaxOrMin = MaxOrMinvalue[ThreadNo];
                TotalIndexValue = IndexValue[ThreadNo];
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                if (MinMaxPointer != 0) {
                    // Note - MPI Call - Allreduce - MPIReducePlusIndex - max with index
                    salsa.mpi.MPIReducePlusIndex result = SALSAUtility.mpiOps.allReduce(
                            new salsa.mpi.MPIReducePlusIndex(TotalIndexValue, TotalMaxOrMin),
                            salsa.mpi.MPIReducePlusIndex.Op.MAX_WITH_INDEX);
                    TotalMaxOrMin = result.getValue();
                    TotalIndexValue = result.getIndex();
                } else {
                    // Note - MPI Call - Allreduce - MPIReducePlusIndex - min with index
                    salsa.mpi.MPIReducePlusIndex result = SALSAUtility.mpiOps.allReduce(
                            new salsa.mpi.MPIReducePlusIndex(TotalIndexValue, TotalMaxOrMin),
                            salsa.mpi.MPIReducePlusIndex.Op.MIN_WITH_INDEX);
                    TotalMaxOrMin = result.getValue();
                    TotalIndexValue = result.getIndex();
                }
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
        }
    }

    /**
     * Finds top LimitNumberStored points by minimizing value given
     * in addAPoint
     * Rsults returned in order of figure of merit
     * Store values and indices
     * Uses FindMinimumSet
     * Results TotalNumberofPoints OrderedMinValue OrderedIndexValue
     */
    public static class FindManyMinValuewithIndex {

        public double[] NumberofPoints;
        public int NumberofThreads;
        public int Numbertofind;
        public double[][] MinValuebythread;
        public int[][] IndexValuebythread;
        public int[] CurrentWorstbythread;

        public double TotalNumberofPoints;
        public double[] TotalMinValue;
        public int[] TotalIndexValue;
        public int TotalWorst;
        public double[] OrderedMinValue;
        public int[] OrderedIndexValue;

        public FindManyMinValuewithIndex(int NumThreads, int LimitNumberStored) {
            NumberofThreads = NumThreads;
            Numbertofind = LimitNumberStored;

            NumberofPoints = new double[NumThreads];
            MinValuebythread = new double[NumThreads][];
            IndexValuebythread = new int[NumThreads][];
            CurrentWorstbythread = new int[NumThreads];

            TotalNumberofPoints = 0.0;
            TotalMinValue = new double[LimitNumberStored];
            TotalIndexValue = new int[LimitNumberStored];
            OrderedMinValue = new double[LimitNumberStored];
            OrderedIndexValue = new int[LimitNumberStored];

            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                CurrentWorstbythread[ThreadNo] = -1;
                MinValuebythread[ThreadNo] = new double[LimitNumberStored];
                IndexValuebythread[ThreadNo] = new int[LimitNumberStored];
                for (int storeloop = 0; storeloop < LimitNumberStored; storeloop++) {
                    MinValuebythread[ThreadNo][storeloop] = -1.0;
                    IndexValuebythread[ThreadNo][storeloop] = -1;
                }
            }
        }

        public final void addAPoint(int ThreadNo, int indexposition, double value) {
            NumberofPoints[ThreadNo] += 1.0;
            tangible.RefObject<Integer> tempRef_Object =
                    new tangible.RefObject<>(CurrentWorstbythread[ThreadNo]);
            FindMinimumSet(value, indexposition, tempRef_Object, MinValuebythread[ThreadNo],
                           IndexValuebythread[ThreadNo], Numbertofind);
            CurrentWorstbythread[ThreadNo] = tempRef_Object.argValue;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {

            for (int storeloop = 0; storeloop < Numbertofind; storeloop++) {
                TotalMinValue[storeloop] = -1.0;
                TotalIndexValue[storeloop] = -1;
            }
            TotalWorst = -1;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                for (int storeloop = 0; storeloop < Numbertofind; storeloop++) {
                    if (IndexValuebythread[ThreadNo][storeloop] < 0) {
                        continue; // End this thread
                    }
                    tangible.RefObject<Integer> tempRef_TotalWorst = new tangible.RefObject<>(TotalWorst);
                    FindMinimumSet(MinValuebythread[ThreadNo][storeloop], IndexValuebythread[ThreadNo][storeloop],
                                   tempRef_TotalWorst, TotalMinValue, TotalIndexValue, Numbertofind);
                    TotalWorst = tempRef_TotalWorst.argValue;
                }
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
            // Sort in absolute order and accumulate over processes. This takes Numbertofindsteps
            for (int OrderLoop = 0; OrderLoop < Numbertofind; OrderLoop++) {
                int localindex = -1; // unset
                double localvalue = -1.0;
                int loopused = -1;
                for (int internalloop = 0; internalloop < Numbertofind; internalloop++) { // Find minimum
                    if (TotalIndexValue[internalloop] < 0) {
                        continue;
                    }
                    if ((localindex < 0) || (TotalMinValue[internalloop] < localvalue)) {
                        localindex = TotalIndexValue[internalloop];
                        localvalue = TotalMinValue[internalloop];
                        loopused = internalloop;
                    }
                }
                int oldlocalindex = localindex;
                if (SALSAUtility.MPI_Size > 1) {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    // Note - MPI Call - Allreduce - MPIReducePlusIndex - min with index
                    salsa.mpi.MPIReducePlusIndex result = SALSAUtility.mpiOps
                                                                      .allReduce(new salsa.mpi.MPIReducePlusIndex(
                                                                              localindex, localvalue),
                                                                                 salsa.mpi.MPIReducePlusIndex.Op
                                                                                         .MIN_WITH_INDEX);
                    localvalue = result.getValue();
                    localindex = result.getIndex();
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }

                OrderedMinValue[OrderLoop] = localvalue;
                OrderedIndexValue[OrderLoop] = localindex;
                if ((oldlocalindex >= 0) && (OrderedIndexValue[OrderLoop] == oldlocalindex)) {
                    TotalIndexValue[loopused] = -1;
                    TotalMinValue[loopused] = -1.0;
                }
            } // Loop over Order Loop
        }
    }

    /**
     * Support finding list of minimum values by inserting new point with value newvalue and index newindex into
     * lists SmallValues SmallIndices
     *  In SmallIndices negative values correspond to unset values
     *  NumberSmallOnes is total number wanted
     *  Currentcut is position in SmallValues, SmallIndices of largest min value
     */
    public static void FindMinimumSet(double newvalue, int newindex, tangible.RefObject<Integer> currentcut,
                                      double[] SmallValues, int[] SmallIndices, int NumberSmallones) {
        if (currentcut.argValue < 0) {
            currentcut.argValue = 0;
            SmallValues[0] = newvalue;
            SmallIndices[0] = newindex;
            return;
        }
        if (SmallIndices[NumberSmallones - 1] < 0) { // Not all positions are filled so add at next available
            // Reset currentcut if worst
            for (int ivalue = 0; ivalue < NumberSmallones; ivalue++) {
                if (SmallIndices[ivalue] < 0) {
                    SmallValues[ivalue] = newvalue;
                    SmallIndices[ivalue] = newindex;
                    if (SmallValues[ivalue] > SmallValues[currentcut.argValue]) {
                        currentcut.argValue = ivalue;
                    }
                    return;
                }
            }
        }
        if (newvalue >= SmallValues[currentcut.argValue]) {
            return;
        }

        // Replace currentcut position with new values and Reset new worst position
        SmallValues[currentcut.argValue] = newvalue;
        SmallIndices[currentcut.argValue] = newindex;
        double maxvalue = -1.0;
        for (int ivalue = 0; ivalue < NumberSmallones; ivalue++) {
            if (SmallIndices[ivalue] < 0) {
                continue;
            }
            if (SmallValues[ivalue] > maxvalue) {
                currentcut.argValue = ivalue;
                maxvalue = SmallValues[ivalue];
            }
        }
    }

    public static class FindCorrelation {
        public double[] NumberofPoints;
        public int NumberofThreads;
        public double[] mean1;
        public double[] mean2;
        public double[] square1;
        public double[] square2;
        public double[] cross12;
        public double TotalNumberofPoints;
        public double Totalmean1;
        public double Totalmean2;
        public double Totalsquare1;
        public double Totalsquare2;
        public double Totalcross12;
        public double Totalsigma1;
        public double Totalsigma2;

        public FindCorrelation(int NumThreads) {
            NumberofThreads = NumThreads;

            NumberofPoints = new double[NumThreads];
            mean1 = new double[NumThreads];
            mean2 = new double[NumThreads];
            square1 = new double[NumThreads];
            square2 = new double[NumThreads];
            cross12 = new double[NumThreads];

            TotalNumberofPoints = 0.0;
            Totalmean1 = 0.0;
            Totalmean2 = 0.0;
            Totalsquare1 = 0.0;
            Totalsquare2 = 0.0;
            Totalcross12 = 0.0;
            Totalsigma1 = 0.0;
            Totalsigma2 = 0.0;
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                NumberofPoints[ThreadNo] = 0.0;
                mean1[ThreadNo] = 0.0;
                mean2[ThreadNo] = 0.0;
                square1[ThreadNo] = 0.0;
                square2[ThreadNo] = 0.0;
                cross12[ThreadNo] = 0.0;
            }
        }

        public final void addAPoint(int ThreadNo, double value1, double value2) {
            NumberofPoints[ThreadNo] += 1.0;
            mean1[ThreadNo] += value1;
            mean2[ThreadNo] += value2;
            square1[ThreadNo] += value1 * value1;
            square2[ThreadNo] += value2 * value2;
            cross12[ThreadNo] += value1 * value2;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++) {
                TotalNumberofPoints += NumberofPoints[ThreadNo];
                Totalmean1 += mean1[ThreadNo];
                Totalmean2 += mean2[ThreadNo];
                Totalsquare1 += square1[ThreadNo];
                Totalsquare2 += square2[ThreadNo];
                Totalcross12 += cross12[ThreadNo];
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Totalmean1 = SALSAUtility.mpiOps.allReduce(Totalmean1, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Totalmean2 = SALSAUtility.mpiOps.allReduce(Totalmean2, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Totalsquare1 = SALSAUtility.mpiOps.allReduce(Totalsquare1, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Totalsquare2 = SALSAUtility.mpiOps.allReduce(Totalsquare2, MPI.SUM);
                // Note - MPI Call - Allreduce - double - sum
                Totalcross12 = SALSAUtility.mpiOps.allReduce(Totalcross12, MPI.SUM);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }

            if (TotalNumberofPoints < 0.5) {
                return;
            }

            Totalmean1 = Totalmean1 / TotalNumberofPoints;
            Totalmean2 = Totalmean2 / TotalNumberofPoints;
            Totalsquare1 = (Totalsquare1 / TotalNumberofPoints) - Totalmean1 * Totalmean1;
            Totalsquare2 = (Totalsquare2 / TotalNumberofPoints) - Totalmean2 * Totalmean2;
            Totalcross12 = (Totalcross12 / TotalNumberofPoints) - Totalmean1 * Totalmean2;
            Totalsigma1 = Math.sqrt(Totalsquare1);
            Totalsigma2 = Math.sqrt(Totalsquare2);
            Totalcross12 = Totalcross12 / (Totalsigma1 * Totalsigma2);
        }

        public final void print(String label, String FPformat) {
            if ((SALSAUtility.DebugPrintOption == 0) || (SALSAUtility.MPI_Rank != 0)) {
                return;
            }
            SALSAUtility.SALSAPrint(1, label + " means " + String.format(FPformat, Totalmean1) + " " + String
                    .format(FPformat, Totalmean2) + " sigmas " + String.format(FPformat, Totalsigma1) + " " + String
                    .format(
                            FPformat, Totalsigma2) + " correl " + String.format(FPformat, Totalcross12));
        }
    }


    public static class Find2DDoubleArraySum {
        public double[] NumberOfPoints;
        public int NumberOfThreads;
        public int OuterDimension;
        public int InnerDimension;
        public double[][][] Sum; // NumberOfThreads x OuterDimension x InnerDimension
        public double TotalNumberofPoints;
        public double[][] TotalSum; // OuterDimension x InnerDimension

        public Find2DDoubleArraySum(int numThreads, int outerDimension, int innerDimension) {
            NumberOfThreads = numThreads;
            OuterDimension = outerDimension;
            InnerDimension = innerDimension;

            NumberOfPoints = new double[numThreads];
            TotalSum = new double[outerDimension][];
            Sum = new double[numThreads][][];

            TotalNumberofPoints = 0.0;
            for (int i = 0; i < outerDimension; ++i) {
                TotalSum[i] = new double[innerDimension];
                for (int j = 0; j < innerDimension; ++j) {
                    TotalSum[i][j] = 0.0;
                }
            }
        }

        public final void startThread(int threadNo) {
            NumberOfPoints[threadNo] = 0.0;
            Sum[threadNo] = new double[OuterDimension][];
            for (int i = 0; i < OuterDimension; ++i) {
                Sum[threadNo][i] = new double[InnerDimension];
                for (int j = 0; j < InnerDimension; ++j) {
                    Sum[threadNo][i][j] = 0.0;
                }
            }
        }

        public final void addAPoint(int threadNo, int row, int col) {
            if (row < 0 || row >= OuterDimension || col < 0 || col >= InnerDimension) {
                return;
            }
            NumberOfPoints[threadNo] += 1.0;
            Sum[threadNo][row][col] += 1.0;
        }

        public final void sumOverThreadsAndMPI() throws MPIException {
            for (int threadNo = 0; threadNo < NumberOfThreads; threadNo++) {
                TotalNumberofPoints += NumberOfPoints[threadNo];
                for (int i = 0; i < OuterDimension; ++i) {
                    for (int j = 0; j < InnerDimension; ++j) {
                        TotalSum[i][j] += Sum[threadNo][i][j];
                    }
                }
            }
            if (SALSAUtility.MPI_Size > 1) {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                // Note - MPI Call - Allreduce - double - sum
                TotalNumberofPoints = SALSAUtility.mpiOps.allReduce(TotalNumberofPoints, MPI.SUM);
                for (int i = 0; i < OuterDimension; ++i) {
                    // Note - MPI Call - Allreduce - double[] - sum
                    SALSAUtility.mpiOps.allReduce(TotalSum[i], MPI.SUM);
                }
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
            }
        }
    }
}