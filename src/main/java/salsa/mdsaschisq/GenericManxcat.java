package salsa.mdsaschisq;

import mpi.MPI;
import org.jblas.*;

//  Routines for Generic Manxcat with parallel computation of components for small number of parameters
public class GenericManxcat {
    public static ChisqFirstandSecond[] LocalAccum;
    public static ChisqFirstandSecond SummedoverThreadsAccum, SummedoverProcessesAccum;

    public static void SetupGenericManxcat() {
        LocalAccum = new ChisqFirstandSecond[SALSAUtility.ThreadCount];
        int NumberParms = Hotsun.npar;
        for (int ThreadNo = 0; ThreadNo < SALSAUtility.ThreadCount; ThreadNo++) {
            LocalAccum[ThreadNo] = new ChisqFirstandSecond(NumberParms);
        }
        SummedoverThreadsAccum = new ChisqFirstandSecond(NumberParms);
        SummedoverProcessesAccum = new ChisqFirstandSecond(NumberParms);

        Hotsun.Number_VectorParameters = Hotsun.npar;
        Hotsun.ParameterVectorDimension = 1;
        Hotsun.DecomposeParameters = false;
        SALSAUtility.sequentialBLAS = true;
        Hotsun.fullmatrixset = true;

    } // End SetupGenericManxcat()

    /**
     * Well Loved Function accum for threading
     * Accum called for each point in Chisq Sum with normation Chisq =  sum(points) value**2
     * deriv[i] is d(value)/d(parameter i)
     *
     * @param ThreadNo
     * @param value
     * @param deriv
     */
    public static void Accum(int ThreadNo, double value, double[] deriv) {

        LocalAccum[ThreadNo].chisq += value * value;
        int linearcount = 0;
        for (int iparm1 = 0; iparm1 < Hotsun.npar; iparm1++) {
            LocalAccum[ThreadNo].first[iparm1] += value * deriv[iparm1];
            for (int iparm2 = iparm1; iparm2 < Hotsun.npar; iparm2++) {
                LocalAccum[ThreadNo].second[linearcount] += deriv[iparm1] * deriv[iparm2];
                ++linearcount;
            }
        }
    }

    public static void AccumDebug(int index) {
        // SALSAUtility.SALSAPrint(0, index.ToString() + " Last two second derivs " + LocalAccum[0].second[26]
        // .ToString("e4") + " " + LocalAccum[0].second[27].ToString("e4"));
    }


    public static void ReInitializeAccum() {
        for (int ThreadNo = 0; ThreadNo < SALSAUtility.ThreadCount; ThreadNo++) {
            ChisqFirstandSecond.initializeToZero(LocalAccum[ThreadNo]);
        }
        ChisqFirstandSecond.initializeToZero(SummedoverThreadsAccum);
        ChisqFirstandSecond.initializeToZero(SummedoverProcessesAccum);
    }

    public static void AddupChisqContributions(Desertwind TotalSolution) {
        // Sum over Threads
        int NumberParms = SummedoverThreadsAccum.VectorSize;
        int SecondSize = SummedoverThreadsAccum.SecondSize;
        for (int ThreadNo = 0; ThreadNo < SALSAUtility.ThreadCount; ThreadNo++) {
            SummedoverThreadsAccum.chisq += LocalAccum[ThreadNo].chisq;
            for (int iparm = 0; iparm < NumberParms; iparm++) {
                SummedoverThreadsAccum.first[iparm] += LocalAccum[ThreadNo].first[iparm];
            }
            for (int iparm = 0; iparm < SecondSize; iparm++) {
                SummedoverThreadsAccum.second[iparm] += LocalAccum[ThreadNo].second[iparm];
            }
        }

        //  Sum over MPI Processes
        ChisqFirstandSecond UsethisAccum = SummedoverThreadsAccum;
        if (SALSAUtility.MPI_Size > 1) {
            UsethisAccum = SummedoverProcessesAccum;
            // Note - MPI Call - Allreduce - double - sum
            SummedoverProcessesAccum.chisq = SALSAUtility.mpiOps.allReduce(SummedoverThreadsAccum.chisq, MPI.SUM);
            // Note - MPI Call - Allreduce - double[] - sum
            SALSAUtility.mpiOps.allReduce(SummedoverThreadsAccum.first, MPI.SUM);
            // Note - MPI Call - Allreduce - double[] - sum
            SALSAUtility.mpiOps.allReduce(SummedoverThreadsAccum.second, MPI.SUM);
        }

        //  Put into Hotsun Form
        TotalSolution.Chisquared = UsethisAccum.chisq;
        Hotsun.zerocr = UsethisAccum.chisq;
        int linearcount = 0;
        for (int iparm1 = 0; iparm1 < NumberParms; iparm1++) {
            TotalSolution.first[iparm1][0] = UsethisAccum.first[iparm1];
            for (int iparm2 = iparm1; iparm2 < NumberParms; iparm2++) {
                double tmp = UsethisAccum.second[linearcount];
                linearcount++;
                TotalSolution.FullMatrix[iparm1][iparm2][0][0] = tmp;
                TotalSolution.FullMatrix[iparm2][iparm1][0][0] = tmp;
                TotalSolution.ExactFullMatrix[iparm1][iparm2][0][0] = tmp;
                TotalSolution.ExactFullMatrix[iparm2][iparm1][0][0] = tmp;
            }
        }
        // SALSAUtility.SALSAPrint(0, "Second Deriv " + TotalSolution.FullMatrix[5, 6][0, 0].ToString("E4"));
    }

    public static void FindQlimits(Desertwind Solution, tangible.RefObject<Double> Qhigh,
                                   tangible.RefObject<Double> Qlow, tangible.RefObject<Integer> ReasontoStop1,
                                   tangible.RefObject<Integer> ReasontoStop2) {
        if (Hotsun.FullSecondDerivative) {
            double[] traceAndNorm = FindTraceandNorm(Solution.ExactFullMatrix);
            Hotsun.ChisqMatrixTrace = traceAndNorm[0];
            Hotsun.ChisqMatrixNorm = traceAndNorm[1];
        } else {
            double[] traceAndNorm = FindTraceandNorm(Solution.FullMatrix);
            Hotsun.ChisqMatrixTrace = traceAndNorm[0];
            Hotsun.ChisqMatrixNorm = traceAndNorm[1];
        }

        double[][] ConventionalMatrix = new double[Hotsun.npar][Hotsun.npar];

        // Set Up Matrix to find eigenvalues
        // Scale but do NOT add Q
        double[][][][] Matrix;
        if (Hotsun.FullSecondDerivative) {
            Matrix = Solution.ExactFullMatrix;
        } else {
            Matrix = Solution.FullMatrix;
        }

        //  Set actual matrix
        for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++) {
            for (int GlobalIndex2 = 0; GlobalIndex2 < Hotsun.npar; GlobalIndex2++) {
                double MatrixElement = Matrix[GlobalIndex1][GlobalIndex2][0][0];
                if (Hotsun.UseDiagonalScaling) {
                    MatrixElement *= Hotsun.sqdginv[GlobalIndex1][0] * Hotsun.sqdginv[GlobalIndex2][0];
                }
                ConventionalMatrix[GlobalIndex1][GlobalIndex2] = MatrixElement;
            } // End GlobalIndex2

        } // End GlobalIndex1
        // Find Minimum and Maximum eigenvalue of ConventionalMatrix

        // Note - jblas - eigen values
        DoubleMatrix lmatrix = new DoubleMatrix(ConventionalMatrix);
        ComplexDoubleMatrix eigenValueDecomp = Eigen.eigenvalues(lmatrix);

        double minEigenValue = Double.MAX_VALUE;
        double maxEigenValue = Double.MIN_VALUE;

        //Assuming you want on the real part...
        for (int i = 0; i < eigenValueDecomp.length; i++) {
            ComplexDouble eigenValue = eigenValueDecomp.get(i);
            if (eigenValue.isReal()) {
                double realValue = eigenValue.real();
                minEigenValue = Math.min(minEigenValue, realValue);
                maxEigenValue = Math.max(maxEigenValue, realValue);
            }
        }
        ReasontoStop1.argValue = 1;
        ReasontoStop2.argValue = 1;
        Qlow.argValue = minEigenValue;
        Qhigh.argValue = maxEigenValue;
    }

    public static double[] FindTraceandNorm(double[][][][] Matrix) {
        double trace = 0.0;
        double norm = 0.0;
        for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++) {
            if (Hotsun.FixedParameter[GlobalIndex1][0]) {
                continue;
            }
            double RowNorm = 0.0;
            for (int GlobalIndex2 = 0; GlobalIndex2 < Hotsun.npar; GlobalIndex2++) {
                if (Hotsun.FixedParameter[GlobalIndex2][0]) {
                    continue;
                }
                double MatrixElement = Matrix[GlobalIndex1][GlobalIndex2][0][0];
                if (Hotsun.UseDiagonalScaling) {
                    MatrixElement *= Hotsun.sqdginv[GlobalIndex1][0] * Hotsun.sqdginv[GlobalIndex2][0];
                }
                if (GlobalIndex1 == GlobalIndex2) {
                    if (Hotsun.AddMarquardtQDynamically) {
                        MatrixElement += Hotsun.Q;
                    }
                    trace += Math.abs(MatrixElement);
                }
                RowNorm += Math.abs(MatrixElement);

            } // End GlobalIndex2

            if (norm < RowNorm) {
                norm = RowNorm;
            }

        }
        return new double[]{trace, norm};
    }

    public static boolean SolveMatrix(double[][] Answer, Desertwind Solution) {
        double[][] ConventionalMatrix = new double[Hotsun.npar][Hotsun.npar];
        double[][] ConventionalFirst = new double[Hotsun.npar][1];
        double[][] ConventionalAnswer = new double[Hotsun.npar][];
        for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++) {
            ConventionalAnswer[GlobalIndex1] = new double[1];
        }

        double[][][][] Matrix;
        if (Hotsun.FullSecondDerivative) {
            Matrix = Solution.ExactFullMatrix;
        } else {
            Matrix = Solution.FullMatrix;
        }

        //  Set actual matrix
        for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++) {
            if (Hotsun.FixedParameter[GlobalIndex1][0]) {
                ConventionalFirst[GlobalIndex1][0] = 0.0;
            } else {
                ConventionalFirst[GlobalIndex1][0] = Solution.first[GlobalIndex1][0] * Hotsun.sqdginv[GlobalIndex1][0];
            }

            for (int GlobalIndex2 = 0; GlobalIndex2 < Hotsun.npar; GlobalIndex2++) {
                double MatrixElement = Matrix[GlobalIndex1][GlobalIndex2][0][0];
                if (Hotsun.UseDiagonalScaling) {
                    MatrixElement *= Hotsun.sqdginv[GlobalIndex1][0] * Hotsun.sqdginv[GlobalIndex2][0];
                }
                if (GlobalIndex1 == GlobalIndex2) {
                    if (Hotsun.AddMarquardtQDynamically) {
                        MatrixElement += Hotsun.Q;
                    }
                }
                ConventionalMatrix[GlobalIndex1][GlobalIndex2] = MatrixElement;
            }

        }

        // Form ConventionalMatrix-1 * ConventionalFirst

        // Note - jblas - solve
        DoubleMatrix cMatrix = new DoubleMatrix(ConventionalMatrix);
        DoubleMatrix RightHandSide = new DoubleMatrix(ConventionalFirst);
        ConventionalAnswer = Solve.solve(cMatrix, RightHandSide).toArray2();


        for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++) {
            Answer[GlobalIndex1][0] = ConventionalAnswer[GlobalIndex1][0] * Hotsun.sqdginv[GlobalIndex1][0];
        }
        return true;
    }

    /**
     * DistributedVector = Matrix (calculated using GlobalxVector) . GlobalVectoronRight
     * In generic case, Matrix is fully calculated and does not use GlobalxVector
     *
     * @param DistributedVector
     * @param Solution
     * @param useexact
     * @param GlobalxVector
     * @param GlobalVectoronRight
     */
    public static void GlobalMatrixVectorProduct(double[][] DistributedVector, Desertwind Solution, boolean useexact,
                                                 double[][] GlobalxVector, double[][] GlobalVectoronRight) {

        double[][][][] Matrix;
        if (useexact) {
            Matrix = Solution.ExactFullMatrix;
        } else {
            Matrix = Solution.FullMatrix;
        }

        for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++) {
            if (Hotsun.FixedParameter[GlobalIndex1][0]) {
                DistributedVector[GlobalIndex1][0] = 0.0;
            } else {
                double tmp = 0.0;
                for (int GlobalIndex2 = 0; GlobalIndex2 < Hotsun.npar; GlobalIndex2++) {
                    if (Hotsun.FixedParameter[GlobalIndex2][0]) {
                        continue;
                    }
                    double MatrixElement = Matrix[GlobalIndex1][GlobalIndex2][0][0];
                    if (Hotsun.UseDiagonalScaling) {
                        MatrixElement *= Hotsun.sqdginv[GlobalIndex1][0] * Hotsun.sqdginv[GlobalIndex2][0];
                    }
                    if (Hotsun.AddMarquardtQDynamically && (GlobalIndex1 == GlobalIndex2)) {
                        MatrixElement += Hotsun.Q;
                    }
                    tmp += MatrixElement * GlobalVectoronRight[GlobalIndex2][0];
                }
                DistributedVector[GlobalIndex1][0] = tmp;
            }
        }
    }

}