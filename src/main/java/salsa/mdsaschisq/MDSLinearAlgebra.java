package salsa.mdsaschisq;

import Manxcat.*;
import SALSALibrary.*;

//C# TO JAVA CONVERTER NOTE: There is no Java equivalent to C# namespace aliases:
//using UserRoutines = salsa.mdsaschisq.ManxcatMDS;


public class MDSLinearAlgebra
{
	public static double[][] DistributedNewIteratedVector = null; // LHS of Iteration Equation
	public static double[][] DistributedOldIteratedVector; // RHS of Iteration Equation
	public static double[][] GlobalOldIteratedVector; // Global version of RHS of Iteration Equation
	public static double[][] DistributedCGVector_R; // Vector r in Conjugate Gradient Method
	public static double[][] DistributedCGVector_Q; // Vector q in Conjugate Gradient Method
	public static double[][] DistributedCGVector_P; // Vector p in Conjugate Gradient Method
	public static double[][] GlobalCGVector_P; // Global version of Vector p in Conjugate Gradient Method
	public static double[][] GlobalCGVector_R; // Global version of Vector R in Conjugate Gradient Method (only used for debug)

	public static void Initialize()
	{
		int LocalNumberofPoints = SALSAUtility.PointCount_Process;
		DistributedOldIteratedVector = new double[LocalNumberofPoints][];
		DistributedNewIteratedVector = new double[LocalNumberofPoints][];
		DistributedCGVector_R = new double[LocalNumberofPoints][];

		int GlobalNumberofPoints = SALSAUtility.PointCount_Global;
		GlobalOldIteratedVector = new double[GlobalNumberofPoints][];
		GlobalCGVector_R = new double[GlobalNumberofPoints][];

		for (int LongIndex = 0; LongIndex < LocalNumberofPoints; LongIndex++)
		{
			DistributedOldIteratedVector[LongIndex] = new double[Hotsun.ParameterVectorDimension];
			DistributedNewIteratedVector[LongIndex] = new double[Hotsun.ParameterVectorDimension];
			DistributedCGVector_R[LongIndex] = new double[Hotsun.ParameterVectorDimension];
		}
		for (int LongIndex = 0; LongIndex < GlobalNumberofPoints; LongIndex++)
		{
			GlobalOldIteratedVector[LongIndex] = new double[Hotsun.ParameterVectorDimension];
			GlobalCGVector_R[LongIndex] = new double[Hotsun.ParameterVectorDimension];
		}

		// re-use Initialized arrays for Conjugate Gradient
		DistributedCGVector_Q = DistributedNewIteratedVector;
		DistributedCGVector_P = DistributedOldIteratedVector;
		GlobalCGVector_P = GlobalOldIteratedVector;

	} // End Initialize() in MDSLinearAlgebra

	//  MaxIndicator = 0 Find Maximum Eigenvalue of Matrix
	//  MaxIndicator = 1 Find Maximum Eigenvalue of (Maximizr - Matrix)
	public static int PowerIterate(Desertwind Solution, int MaxIndicator, double Maximizer, tangible.RefObject<Double> PowerEigenvalue)
	{
		if (DistributedNewIteratedVector == null)
		{
			Initialize();
		}
		Hotsun.UseDiagonalScaling = true;
		Hotsun.AddMarquardtQDynamically = false;

		//  Set up Initial Power Vectors
		SetInitialPowerVector(DistributedNewIteratedVector);
		double AdotA = SALSABLAS.VectorScalarProduct(DistributedNewIteratedVector, DistributedNewIteratedVector);

		int somethingtodo = 0;
		int PowerIterationCount = 0;
		PowerEigenvalue.argValue = -1.0;
		while (true)
		{ // Iterate over Power Multiplications by Chisq Matrix

			//  Normalize Current A and move from New to Old

			double OldNorm = 1.0 / Math.sqrt(AdotA);
			SALSABLAS.LinearCombineVector(DistributedOldIteratedVector, OldNorm, DistributedNewIteratedVector, 0.0, DistributedNewIteratedVector);

			//  Make a Global Vector of DistributedOldIteratedVector
			ManxcatCentral.MakeVectorGlobal(DistributedOldIteratedVector, GlobalOldIteratedVector);

			//  Form Chisq Matrix Product with Old Vector
			if (Hotsun.FullSecondDerivative)
			{
				ManxcatMDS.GlobalMatrixVectorProduct(DistributedNewIteratedVector, Solution, true, Hotsun.GlobalParameter, GlobalOldIteratedVector);
			}
			else
			{
				ManxcatMDS.GlobalMatrixVectorProduct(DistributedNewIteratedVector, Solution, false, Hotsun.GlobalParameter, GlobalOldIteratedVector);
			}

			//  Correct case MaxIndicator = 1
			if (MaxIndicator > 0)
			{
				SALSABLAS.LinearCombineVector(DistributedNewIteratedVector, -1.0, DistributedNewIteratedVector, Maximizer, DistributedOldIteratedVector);
			}
			SALSABLAS.LinearCombineVector(DistributedNewIteratedVector, 1.0, DistributedNewIteratedVector, Hotsun.addonforQcomputation, DistributedOldIteratedVector);

			//  Form Scalar Products
			double NewEigenvalue = SALSABLAS.VectorScalarProduct(DistributedNewIteratedVector, DistributedOldIteratedVector);
			AdotA = SALSABLAS.VectorScalarProduct(DistributedNewIteratedVector, DistributedNewIteratedVector);

			++PowerIterationCount;
			somethingtodo = -1;
			if (SALSAUtility.MPI_Rank == 0)
			{
				if (PowerIterationCount > 10 && (NewEigenvalue > 0.0))
				{ // Arbitary criteria for starting +tests
					somethingtodo = 0;
					double scaleit = 1.0;
					if (MaxIndicator > 0)
					{
						scaleit = Hotsun.extraprecision;
					}
					if (Math.abs(NewEigenvalue - PowerEigenvalue.argValue) > PowerEigenvalue.argValue * scaleit * Hotsun.eigenvaluechange)
					{
						++somethingtodo;
					}
					double delta = AdotA - NewEigenvalue * NewEigenvalue; // (Ax- Eigenvalue*Axold)**2
					if (Math.abs(delta) > NewEigenvalue * NewEigenvalue * scaleit * Hotsun.eigenvectorchange)
					{
						++somethingtodo;
					}
				}
			}
			PowerEigenvalue.argValue = NewEigenvalue;

			if (SALSAUtility.MPI_Size > 1)
			{
				SALSAUtility.StartSubTimer(SALSAUtility.MPIBROADCASTTiming);
				tangible.RefObject<Integer> tempRef_somethingtodo = new tangible.RefObject<Integer>(somethingtodo);
				SALSAUtility.MPI_communicator.<Integer>Broadcast(tempRef_somethingtodo, 0);
				somethingtodo = tempRef_somethingtodo.argValue;
				SALSAUtility.StopSubTimer(SALSAUtility.MPIBROADCASTTiming);
			}
			if (PowerIterationCount >= Hotsun.PowerIterationLimit)
			{
				somethingtodo = -2;
				break;
			}
			if (somethingtodo == 0)
			{
				somethingtodo = PowerIterationCount;
				break;
			}

		} // End while over PowerIterationCounts

		return somethingtodo;

	} // End PowerIterate()

	//  Solve for Distributed Answer = (Matrix-1) First 
	//  Note RealMatrixSize takes account of forced zero parameters which does not crop up in explicit
	//  algorithm as implemented in Matrix Vector Multiplier
	public static boolean ConjugateGradientSolver(double[][] Answer, Desertwind Solution, boolean useexact, double[][] GlobalxVector, double[][] DistributedRHS, tangible.RefObject<Integer> NumberofIterations, int RealMatrixSize, double LimitonNormofR)
	{

		boolean matrixsuccess = true;

		//  Initialize
		SALSABLAS.zrword(Answer); // Zero Solution x called xshift in Manxcat and stored in Answer
		SALSABLAS.CopyVector(DistributedCGVector_R, DistributedRHS, 0, SALSAUtility.PointCount_Process); // Set R(0) = RHS
		double InitialRNorm = SALSABLAS.VectorScalarProduct(DistributedCGVector_R, DistributedCGVector_R);
		double RNormTest = InitialRNorm * LimitonNormofR; // Limit for Test on Iteration of Norm of R
		double[] SaveRNorm = new double[RealMatrixSize + 1];
		SaveRNorm[0] = InitialRNorm;
		int CountSteps = 0;

		double lastrho = 1.0;
		double currentrho = 1.0;

		//  Loop over Conjugate Gradient Steps
		while (true)
		{
			// Set value of rho
			lastrho = currentrho;
			currentrho = SALSABLAS.VectorScalarProduct(DistributedCGVector_R, DistributedCGVector_R);

			// Set Vector P
			++CountSteps;
			if (CountSteps == 1)
			{
				SALSABLAS.CopyVector(DistributedCGVector_P, DistributedCGVector_R, 0, SALSAUtility.PointCount_Process);
			}
			else
			{
				SALSABLAS.LinearCombineVector(DistributedCGVector_P, currentrho / lastrho, DistributedCGVector_P, 1.0, DistributedCGVector_R);
			}

			//  Make a Global Vector of DistributedCGVector_P
			ManxcatCentral.MakeVectorGlobal(DistributedCGVector_P, GlobalCGVector_P);

			//  Distributed Q = Matrix . Global P
			ManxcatMDS.GlobalMatrixVectorProduct(DistributedCGVector_Q, Solution, useexact, Hotsun.GlobalParameter, GlobalCGVector_P);

			//  New Answer is Old answer + (Current Rho / (P dot Q)) Vector P
			double PdotQ = SALSABLAS.VectorScalarProduct(DistributedCGVector_P, DistributedCGVector_Q);
			double alpha = currentrho / PdotQ;
			SALSABLAS.LinearCombineVector(Answer, alpha, DistributedCGVector_P, 1.0, Answer);

			// New residual R = Old Residual - (Current Rho / (P dot Q)) Vector Q
			SALSABLAS.LinearCombineVector(DistributedCGVector_R, -alpha, DistributedCGVector_Q, 1.0, DistributedCGVector_R);

			//  See if we can or should End
			double CurrentRNorm = SALSABLAS.VectorScalarProduct(DistributedCGVector_R, DistributedCGVector_R);
			SaveRNorm[CountSteps] = CurrentRNorm;
			boolean TestRNorm = CurrentRNorm <= RNormTest;
			tangible.RefObject<Boolean> tempRef_TestRNorm = new tangible.RefObject<Boolean>(TestRNorm);
			SALSAUtility.SynchronizeMPIvariable(tempRef_TestRNorm);
			TestRNorm = tempRef_TestRNorm.argValue;
			if (TestRNorm)
			{
				matrixsuccess = true;
				break;
			}

			// Singular Matrix
			if (CountSteps >= RealMatrixSize)
			{
				matrixsuccess = false;
				ManxcatCentral.MakeVectorGlobal(DistributedCGVector_R, GlobalCGVector_R);
				if (SALSAUtility.MPI_Rank == 0)
				{
					SALSAUtility.SALSAPrint(0, " CG Failure after " + (new Integer(RealMatrixSize)).toString() + " Steps");
					String ListofNorms = "";
					int Normindex = CountSteps;
					for (int inorm = 0; inorm < 10; inorm++)
					{
						ListofNorms += " " + String.format("%0.4E", SaveRNorm[Normindex]);
						--Normindex;
						if (Normindex < 0)
						{
							break;
						}
					}
					SALSAUtility.SALSAPrint(0, "Last 10 Norms " + ListofNorms);

					String fname = ManxcatCentral.ResultDirectoryName + "\\BadCGVector" + (new Integer(Hotsun.TotalCGFailures)).toString();
					StreamWriter sw = new StreamWriter(fname, false, Encoding.UTF8);
					double fractionnorm = 1.0 / Math.sqrt(CurrentRNorm);
					try
					{
						for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
						{
							String Coordinates = "";
							int UsedPointIndex = SALSAUtility.NaivetoActualUsedOrder[GlobalPointIndex];
							double pointsize = 0.0;
							for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
							{
								Coordinates += String.format("%0.4E", GlobalCGVector_R[UsedPointIndex][LocalVectorIndex]) + "\t";
								pointsize += GlobalCGVector_R[UsedPointIndex][LocalVectorIndex] * GlobalCGVector_R[UsedPointIndex][LocalVectorIndex];
							}
							pointsize = Math.sqrt(pointsize) * fractionnorm;
							sw.WriteLine(String.format((GlobalPointIndex).toString() + "\t" + Coordinates + " " + String.format("%0.4E", pointsize)));
						}

						sw.Flush();
						sw.Close();
					}
					catch (RuntimeException e)
					{
						System.out.println("Failed writing data in CG Solver " + e);
						throw (e);
					}
				}
				break;
			}
		} // End Loop over Conjugate Gradient Steps

		NumberofIterations.argValue = CountSteps;
		return matrixsuccess;

	} // End ConjugateGradientSolver(Desertwind Solution)

	// Scale Distributed Vector by a Global Vector
	public static void DiagScaleVector(double[][] ScaledDistributedVector, double[][] OldDistributedVector, double[][] GlobalScalingVector)
	{
		int LocalVectorDimension = ScaledDistributedVector[0].GetLength(0);

		// Parallel Setting
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
					{
						ScaledDistributedVector[LongIndex][LocalVectorIndex] = 0.0;
					}
					else
					{
						ScaledDistributedVector[LongIndex][LocalVectorIndex] = OldDistributedVector[LongIndex][LocalVectorIndex] * GlobalScalingVector[GlobalIndex][LocalVectorIndex];
					}
				}
			}
		}
	   );

	} // // End DiagScaleVector


	// Set Initial Power Vector to be a random number with normalization 1
	// Be sure to set fixed variables to be zero
	public static void SetInitialPowerVector(double[][] TobeSet)
	{
		int LocalVectorDimension = TobeSet[0].GetLength(0);
		java.util.Random Randobject = new java.util.Random();

		// Parallel Setting
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
					{
						TobeSet[LongIndex][LocalVectorIndex] = 0.0;
					}
					else
					{
						TobeSet[LongIndex][LocalVectorIndex] = Randobject.nextDouble();
					}
				}
			}
		}
	   );

	} // // End SetInitialPowerVector(double[][] TobeSet)



	//  Calculate Matrix Global Vector product storing as a distributed vector
	public static void FindTraceandNorm(double[][][] MatrixDiagonals, double[][] GlobalxVector, tangible.RefObject<Double> Trace, tangible.RefObject<Double> Norm)
	{
		Trace.argValue = 0.0;
		Norm.argValue = 0.0;
		GlobalReductions.FindDoubleSum FindTrace = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);
		GlobalReductions.FindDoubleMax FindNorm = new GlobalReductions.FindDoubleMax(SALSAUtility.ThreadCount);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
//	Start Code to calculate Trace and Norm
			double LocalTrace = 0.0;
			double LocalNorm = 0.0;
			double WeightFunction1, WeightFunction2 = 0;
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex1 = beginpoint; LongIndex1 < indexlen + beginpoint; LongIndex1++)
			{
				int GlobalIndex1 = LongIndex1 + SALSAUtility.PointStart_Process;
				for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < Hotsun.ParameterVectorDimension; LocalVectorIndex1++)
				{
					if (Hotsun.FixedParameter[GlobalIndex1][LocalVectorIndex1])
					{
						continue;
					}
						double RowNorm = 0.0;
						for (int GlobalIndex2 = 0; GlobalIndex2 < SALSAUtility.PointCount_Global; GlobalIndex2++)
						{
							tangible.RefObject<Double> tempRef_WeightFunction1 = new tangible.RefObject<Double>(WeightFunction1);
							tangible.RefObject<Double> tempRef_WeightFunction2 = new tangible.RefObject<Double>(WeightFunction2);
							boolean tempVar = !ManxcatMDS.SetChisqWeights(GlobalIndex1, GlobalIndex2, tempRef_WeightFunction1, tempRef_WeightFunction2);
								WeightFunction1 = tempRef_WeightFunction1.argValue;
							WeightFunction2 = tempRef_WeightFunction2.argValue;
							if (tempVar)
							{
								continue;
							}
								double DistanceFudge1 = 1.0;
								double DistanceFudge2 = 1.0;
								double SquaredDistance = 0.0;
								double funcl = 0.0;
								if (GlobalIndex1 != GlobalIndex2)
								{
									for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < Hotsun.ParameterVectorDimension; LocalVectorIndex2++)
									{
										double tmp1 = GlobalxVector[GlobalIndex1][LocalVectorIndex2] - GlobalxVector[GlobalIndex2][LocalVectorIndex2];
										SquaredDistance += tmp1 * tmp1;
									}
									double ActualDistance = SquaredDistance;
									if (ManxcatMDS.DistanceFormula == 1)
									{
										SquaredDistance += ManxcatMDS.SQUAREMinimumDistance;
										ActualDistance = Math.sqrt(SquaredDistance);
										DistanceFudge1 = 0.5 / ActualDistance;
										DistanceFudge2 = DistanceFudge1 / SquaredDistance;
									}
									funcl = WeightFunction1 - WeightFunction2 * ActualDistance;
								}
								for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < Hotsun.ParameterVectorDimension; LocalVectorIndex2++)
								{
									if (Hotsun.FixedParameter[GlobalIndex2][LocalVectorIndex2])
									{
										continue;
									}
										double MatrixElement;
										if (GlobalIndex1 == GlobalIndex2)
										{
											MatrixElement = MatrixDiagonals[LongIndex1][LocalVectorIndex1][LocalVectorIndex2];
											if (Hotsun.UseDiagonalScaling)
											{
												MatrixElement *= Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex1] * Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex2];
											}
												if (Hotsun.AddMarquardtQDynamically && (LocalVectorIndex1 == LocalVectorIndex2))
												{
													MatrixElement += Hotsun.Q;
												}
													if (LocalVectorIndex1 == LocalVectorIndex2)
													{
														LocalTrace += Math.abs(MatrixElement);
													}
										}
										else
										{
											double correction = 0.0;
											double VectorCrossProduct = (GlobalxVector[GlobalIndex1][LocalVectorIndex1] - GlobalxVector[GlobalIndex2][LocalVectorIndex1]) * (GlobalxVector[GlobalIndex1][LocalVectorIndex2] - GlobalxVector[GlobalIndex2][LocalVectorIndex2]);
											MatrixElement = -8.0 * VectorCrossProduct * DistanceFudge1 * DistanceFudge1 * WeightFunction2 * WeightFunction2;
											if (Hotsun.FullSecondDerivative)
											{
												if ((ManxcatMDS.DistanceFormula == 2) && (LocalVectorIndex1 == LocalVectorIndex2))
												{
													correction = 4.0 * funcl * WeightFunction2;
												}
													if ((ManxcatMDS.DistanceFormula == 1) && (LocalVectorIndex1 == LocalVectorIndex2))
													{
														correction = 4.0 * funcl * WeightFunction2 * DistanceFudge1;
													}
														if (ManxcatMDS.DistanceFormula == 1)
														{
															correction += - 4.0 * funcl * VectorCrossProduct * WeightFunction2 * DistanceFudge2;
														}
															MatrixElement += correction;
											}
											if (Hotsun.UseDiagonalScaling)
											{
												MatrixElement *= Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex1] * Hotsun.sqdginv[GlobalIndex2][LocalVectorIndex2];
											}
										}
										RowNorm += Math.abs(MatrixElement);
								}
						}
						if (LocalNorm < RowNorm)
						{
							LocalNorm = RowNorm;
						}
				}
			}
			FindTrace.addapoint(ThreadNo, LocalTrace);
			FindNorm.addapoint(ThreadNo, LocalNorm);
		}
	   );

		FindTrace.sumoverthreadsandmpi();
		FindNorm.sumoverthreadsandmpi();

		Trace.argValue = FindTrace.Total;
		Norm.argValue = FindNorm.TotalMax;
		return;


	} // End FindTraceandNorm(double[][] MatrixDiagonals,double[][] GlobalxVector, out double Trace, out double Norm)



} // End Namespace MDS