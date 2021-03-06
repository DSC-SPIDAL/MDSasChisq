package salsa.mdsaschisq;

import edu.rice.hj.api.SuspendableException;

import static edu.rice.hj.Module1.forallChunked;

public class SALSABLAS {

    // Zero one dimensional double Array
    public static void zrword(double[] TobeZeroed) {
        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = TobeZeroed.length;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                TobeZeroed[LongIndex] = 0.0;
            }
            return;
        }

        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    TobeZeroed[LongIndex] = 0.0;
                }
            }
            );
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

    }

    // Zero two dimensional double Array
    public static void zrword(double[][] TobeZeroed) {
        int LocalVectorDimension = TobeZeroed[0].length;

        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = TobeZeroed.length;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                    TobeZeroed[LongIndex][LocalVectorIndex] = 0.0;
                }
            }
            return;
        }

        // Parallel Initialization
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                        TobeZeroed[LongIndex][LocalVectorIndex] = 0.0;
                    }
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

    }

    // Zero Three dimensional double Array
    public static void zrword(double[][][] TobeZeroed) {
        int LocalVectorDimension1 = TobeZeroed[0].length;
        int LocalVectorDimension2 = TobeZeroed[0][0].length;

        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = TobeZeroed.length;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension1; LocalVectorIndex1++) {
                    for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension2; LocalVectorIndex2++) {
                        TobeZeroed[LongIndex][LocalVectorIndex1][LocalVectorIndex2] = 0.0;
                    }
                }
            }
            return;
        }

        // Parallel Initialization
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension1; LocalVectorIndex1++) {
                        for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension2; LocalVectorIndex2++) {
                            TobeZeroed[LongIndex][LocalVectorIndex1][LocalVectorIndex2] = 0.0;
                        }
                    }
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

    }

    // Zero one dimensional int Array
    public static void zrword(int[] TobeZeroed) {
        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = TobeZeroed.length;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                TobeZeroed[LongIndex] = 0;
            }
            return;
        }

        // Parallel Initialization
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    TobeZeroed[LongIndex] = 0;
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

    }

    // Zero two dimensional int Array
    public static void zrword(int[][] TobeZeroed) {
        int LocalVectorDimension = TobeZeroed[0].length;

        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = TobeZeroed.length;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                    TobeZeroed[LongIndex][LocalVectorIndex] = 0;
                }
            }
            return;
        }

        // Parallel Initialization
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                        TobeZeroed[LongIndex][LocalVectorIndex] = 0;
                    }
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }

    // Set one dimensional Boolean Array
    public static void zrword(boolean[] TobeZeroed, boolean boolvalue) {
        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = TobeZeroed.length;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                TobeZeroed[LongIndex] = boolvalue;
            }
            return;
        }

        // Parallel Initialization
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    TobeZeroed[LongIndex] = boolvalue;
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

    }

    // Set two dimensional Boolean Array
    public static void zrword(boolean[][] TobeZeroed, boolean boolvalue) {
        int LocalVectorDimension = TobeZeroed[0].length;

        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = TobeZeroed.length;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                    TobeZeroed[LongIndex][LocalVectorIndex] = boolvalue;
                }
            }
            return;
        }

        // Parallel Initialization
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                        TobeZeroed[LongIndex][LocalVectorIndex] = boolvalue;
                    }
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }

    public static void zrword(double[][][][] MatrixA) {
        int LocalVectorDimension0 = MatrixA[0][0].length;
        int LocalVectorDimension1 = MatrixA[0][0][0].length;
        int LongDimension0 = MatrixA.length;
        int LongDimension1 = MatrixA[0].length;

        if (!SALSAUtility.sequentialBLAS) {
            SALSAUtility.printAndThrowRuntimeException("zeroing a Matrix NOT defined for Decomposed Parameters");


        }

        for (int LongIndex1 = 0; LongIndex1 < LongDimension0; LongIndex1++) {
            for (int LongIndex2 = 0; LongIndex2 < LongDimension1; LongIndex2++) {
                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension0; LocalVectorIndex1++) {
                    for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension1; LocalVectorIndex2++) {
                        MatrixA[LongIndex1][LongIndex2][LocalVectorIndex1][LocalVectorIndex2] = 0.0;
                    }
                }
            }
        }
    }

    public static void LinearCombineVector(double[] VectorC, double CoeffAlpha, double[] VectorA, double CoeffBeta,
                                           double[] VectorB) {
        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = VectorC.length;

            if (CoeffBeta != 0.0) {
                for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                    VectorC[LongIndex] = CoeffAlpha * VectorA[LongIndex] + CoeffBeta * VectorB[LongIndex];
                }
            } else {
                for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                    VectorC[LongIndex] = CoeffAlpha * VectorA[LongIndex];
                }
            }
            return;
        }

        // Parallel VectorC = CoeffAlpha * VectorA + CoeffBeta * VectorB
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                if (CoeffBeta != 0.0) {
                    for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                        VectorC[LongIndex] = CoeffAlpha * VectorA[LongIndex] + CoeffBeta * VectorB[LongIndex];
                    }
                } else {
                    for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                        VectorC[LongIndex] = CoeffAlpha * VectorA[LongIndex];
                    }
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

    } // End LinearCombineVector -- One dimensional double

    public static void LinearCombineVector(double[][] VectorC, double CoeffAlpha, double[][] VectorA, double CoeffBeta,
                                           double[][] VectorB) {
        int LocalVectorDimension = VectorC[0].length;

        if (LocalVectorDimension != VectorA[0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LocalVectorDimension + " A " + VectorA[0].length);
        }

        if (LocalVectorDimension != VectorB[0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LocalVectorDimension + " B " + VectorB[0].length);
        }

        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = VectorC.length;

            if (CoeffBeta != 0.0) {
                for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                    for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                        VectorC[LongIndex][LocalVectorIndex] = CoeffAlpha * VectorA[LongIndex][LocalVectorIndex] +
                                CoeffBeta * VectorB[LongIndex][LocalVectorIndex];
                    }
                }
            } else {
                for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                    for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                        VectorC[LongIndex][LocalVectorIndex] = CoeffAlpha * VectorA[LongIndex][LocalVectorIndex];
                    }
                }
            }
            return;
        }

        // Parallel VectorC = CoeffAlpha * VectorA + CoeffBeta * VectorB
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                if (CoeffBeta != 0.0) {
                    for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                        for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                            VectorC[LongIndex][LocalVectorIndex] = CoeffAlpha * VectorA[LongIndex][LocalVectorIndex] +
                                    CoeffBeta * VectorB[LongIndex][LocalVectorIndex];
                        }
                    }
                } else {
                    for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                        for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                            VectorC[LongIndex][LocalVectorIndex] = CoeffAlpha * VectorA[LongIndex][LocalVectorIndex];
                        }
                    }
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }

    //  Copy TotalSize local vectors from VectorA to VectorC starting in position 0 of VectorC and position
    // StartIndex of VectorA
    public static void CopyVector(double[][] VectorC, double[][] VectorA, int StartIndex, int TotalSize) {
        int LocalVectorDimension = VectorC[0].length;

        if (LocalVectorDimension != VectorA[0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LocalVectorDimension + " A " + VectorA[0].length);
        }

        if (SALSAUtility.sequentialBLAS) {
            for (int LongIndex = 0; LongIndex < TotalSize; LongIndex++) {
                System.arraycopy(VectorA[LongIndex], 0, VectorC[LongIndex + StartIndex], 0, LocalVectorDimension);
            }
            return;
        }

        // Parallel VectorC = VectorA
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
                    System.arraycopy(VectorA[LongIndex], 0, VectorC[LongIndex + StartIndex], 0, LocalVectorDimension);
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }

    //  Copy TotalSize local vectors from VectorA to VectorC starting in position 0 of VectorC and position
    // StartIndex of VectorA
    public static void CopyVector(String[] VectorC, String[] VectorA, int StartIndex, int TotalSize) {
        if (SALSAUtility.sequentialBLAS) {
            System.arraycopy(VectorA, 0, VectorC, StartIndex, TotalSize);
            return;
        }

        // Parallel VectorC = VectorA
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
                System.arraycopy(VectorA, beginpoint, VectorC, beginpoint + StartIndex, endpoint - beginpoint);
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }

    //  Copy TotalSize local vectors from VectorA to VectorC starting in position 0 of VectorC and position
    // StartIndex of VectorA
    public static void CopyVector(double[][][] VectorC, double[][][] VectorA, int StartIndex, int TotalSize) {
        int LocalVectorDimension1 = VectorC[0].length;

        if (LocalVectorDimension1 != VectorA[0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LocalVectorDimension1 + " A " + VectorA[0].length);
        }

        int LocalVectorDimension2 = VectorC[0][0].length;

        if (LocalVectorDimension2 != VectorA[0][0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LocalVectorDimension2 + " A " + VectorA[0][0].length);
        }

        if (SALSAUtility.sequentialBLAS) {
            for (int LongIndex = 0; LongIndex < TotalSize; LongIndex++) {
                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension1; LocalVectorIndex1++) {
                    System.arraycopy(VectorA[LongIndex][LocalVectorIndex1], 0,
                                     VectorC[LongIndex + StartIndex][LocalVectorIndex1], 0, LocalVectorDimension2);
                }
            }
            return;
        }

        // Parallel VectorC = VectorA
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
                    for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension1; LocalVectorIndex1++) {
                        System.arraycopy(VectorA[LongIndex][LocalVectorIndex1], 0,
                                         VectorC[LongIndex + StartIndex][LocalVectorIndex1], 0, LocalVectorDimension2);
                    }
                }
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }

    //  Form Vector Dot Product -- One Dimensional
    public static double VectorScalarProduct(double[] VectorA, double[] VectorB) {
        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = VectorA.length;

            double Total = 0.0;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                Total += VectorA[LongIndex] * VectorB[LongIndex];
            }
            return Total;
        }


        // Parallel Scalar Product  = Sum(over i) VectorA[i] * VectorB[i]
        GlobalReductions.FindDoubleSum FindScalarProduct = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);

        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                double tmp = 0.0;
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    tmp += VectorA[LongIndex] * VectorB[LongIndex];
                }
                FindScalarProduct.addAPoint(threadIndex, tmp);
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

        FindScalarProduct.sumOverThreadsAndMPI();
        return FindScalarProduct.Total;
    }

    //  Form Vector Dot Product -- Two Dimensional
    public static double VectorScalarProduct(double[][] VectorA, double[][] VectorB) {
        int LocalVectorDimension = VectorB[0].length;

        if (LocalVectorDimension != VectorA[0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions B" + LocalVectorDimension + " A " + VectorA[0].length);
        }

        if (SALSAUtility.sequentialBLAS) {
            int LongDimension = VectorA.length;

            double Total = 0.0;

            for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++) {
                for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                    Total += VectorB[LongIndex][LocalVectorIndex] * VectorA[LongIndex][LocalVectorIndex];
                }
            }
            return Total;
        }

        // Parallel Scalar Product  = Sum(over i) VectorA[i] * VectorB[i]
        GlobalReductions.FindDoubleSum ScalarProduct = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);

        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
            {
                double tmp = 0.0;
                int indexlen = SALSAUtility.PointsperThread[threadIndex];
                int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                    for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++) {
                        tmp += VectorB[LongIndex][LocalVectorIndex] * VectorA[LongIndex][LocalVectorIndex];
                    }
                }
                ScalarProduct.addAPoint(threadIndex, tmp);
            });
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

        ScalarProduct.sumOverThreadsAndMPI();
        return ScalarProduct.Total;
    }

    public static void CopyMatrix(double[][][][] MatrixC, double[][][][] MatrixA) {
        int LocalVectorDimension0 = MatrixC[0][0].length;

        if (LocalVectorDimension0 != MatrixA[0][0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LocalVectorDimension0 + " A " + MatrixA[0][0].length);
        }

        int LocalVectorDimension1 = MatrixC[0][0][0].length;

        if (LocalVectorDimension1 != MatrixA[0][0][0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LocalVectorDimension1 + " A " + MatrixA[0][0][0].length);
        }
        int LongDimension0 = MatrixC.length;

        if (LongDimension0 != MatrixA.length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LongDimension0 + " A " + MatrixA.length);
        }
        int LongDimension1 = MatrixC[0].length;

        if (LongDimension1 != MatrixA[0].length) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Inconsistent Dimensions C" + LongDimension1 + " A " + MatrixA[0].length);
        }

        if (!SALSAUtility.sequentialBLAS) {
            SALSAUtility.printAndThrowRuntimeException("CopyMatrix NOT defined for Decomposed Parameters");


        }

        for (int LongIndex1 = 0; LongIndex1 < LongDimension0; LongIndex1++) {
            for (int LongIndex2 = 0; LongIndex2 < LongDimension1; LongIndex2++) {
                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension0; LocalVectorIndex1++) {
                    System.arraycopy(MatrixA[LongIndex1][LongIndex2][LocalVectorIndex1], 0,
                                     MatrixC[LongIndex1][LongIndex2][LocalVectorIndex1], 0, LocalVectorDimension1);
                }
            }
        }
    }
}