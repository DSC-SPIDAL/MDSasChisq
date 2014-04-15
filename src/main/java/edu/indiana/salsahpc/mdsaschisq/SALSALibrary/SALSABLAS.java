package SALSALibrary;

public class SALSABLAS
{

	// Zero one dimensional double Array
	public static void zrword(double[] TobeZeroed)
	{
		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = TobeZeroed.length;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				TobeZeroed[LongIndex] = 0.0;
			}
			return;
		}

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				TobeZeroed[LongIndex] = 0.0;
			}
		}
	   );

	} // End zrword -- One dimensional double

	// Zero two dimensional double Array
	public static void zrword(double[][] TobeZeroed)
	{
		int LocalVectorDimension = TobeZeroed[0].GetLength(0);

		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = TobeZeroed.length;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					TobeZeroed[LongIndex][LocalVectorIndex] = 0.0;
				}
			}
			return;
		}

		// Parallel Initialization
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					TobeZeroed[LongIndex][LocalVectorIndex] = 0.0;
				}
			}
		}
	   );

	} // // End zrword -- Two dimensional double

	// Zero Three dimensional double Array
	public static void zrword(double[][][] TobeZeroed)
	{
		int LocalVectorDimension1 = TobeZeroed[0].GetLength(0);
		int LocalVectorDimension2 = TobeZeroed[0].GetLength(1);

		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = TobeZeroed.length;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension1; LocalVectorIndex1++)
				{
					for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension2; LocalVectorIndex2++)
					{
						TobeZeroed[LongIndex][LocalVectorIndex1][LocalVectorIndex2] = 0.0;
					}
				}
			}
			return;
		}

		// Parallel Initialization
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension1; LocalVectorIndex1++)
				{
					for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension2; LocalVectorIndex2++)
					{
						TobeZeroed[LongIndex][LocalVectorIndex1][LocalVectorIndex2] = 0.0;
					}
				}
			}
		}
	   );

	} // // End zrword -- Three dimensional double

	// Zero one dimensional int Array
	public static void zrword(int[] TobeZeroed)
	{
		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = TobeZeroed.length;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				TobeZeroed[LongIndex] = 0;
			}
			return;
		}

		// Parallel Initialization
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				TobeZeroed[LongIndex] = 0;
			}
		}
	   );

	} // End zrword -- One dimensional int

	// Zero two dimensional int Array
	public static void zrword(int[][] TobeZeroed)
	{
		int LocalVectorDimension = TobeZeroed[0].GetLength(0);

		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = TobeZeroed.length;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					TobeZeroed[LongIndex][LocalVectorIndex] = 0;
				}
			}
			return;
		}

		// Parallel Initialization
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					TobeZeroed[LongIndex][LocalVectorIndex] = 0;
				}
			}
		}
	   );

	} // // End zrword -- Two dimensional int

	// Set one dimensional Boolean Array
	public static void zrword(boolean[] TobeZeroed, boolean boolvalue)
	{
		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = TobeZeroed.length;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				TobeZeroed[LongIndex] = boolvalue;
			}
			return;
		}

		// Parallel Initialization
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				TobeZeroed[LongIndex] = boolvalue;
			}
		}
	   );

	} // End zrword -- One dimensional Boolean

	// Set two dimensional Boolean Array
	public static void zrword(boolean[][] TobeZeroed, boolean boolvalue)
	{
		int LocalVectorDimension = TobeZeroed[0].GetLength(0);

		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = TobeZeroed.length;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					TobeZeroed[LongIndex][LocalVectorIndex] = boolvalue;
				}
			}
			return;
		}

		// Parallel Initialization
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					TobeZeroed[LongIndex][LocalVectorIndex] = boolvalue;
				}
			}
		}
	   );

	} // // End zrword -- Two dimensional Boolean

	public static void zrword(double[][][][] MatrixA)
	{
		int LocalVectorDimension0 = MatrixA[0][0].GetLength(0);
		int LocalVectorDimension1 = MatrixA[0][0].GetLength(1);
		int LongDimension0 = MatrixA.length;
		int LongDimension1 = MatrixA[0].length;

		if (!SALSAUtility.sequentialBLAS)
		{
			RuntimeException e = SALSAUtility.SALSAError("zeroing a Matrix NOT defined for Decomposed Parameters");

			throw (e);
		}

		for (int LongIndex1 = 0; LongIndex1 < LongDimension0; LongIndex1++)
		{
			for (int LongIndex2 = 0; LongIndex2 < LongDimension1; LongIndex2++)
			{
				for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension0; LocalVectorIndex1++)
				{
					for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension1; LocalVectorIndex2++)
					{
						MatrixA[LongIndex1][LongIndex2][LocalVectorIndex1][LocalVectorIndex2] = 0.0;
					}
				}
			}
		}
	}

	public static void LinearCombineVector(double[] VectorC, double CoeffAlpha, double[] VectorA, double CoeffBeta, double[] VectorB)
	{
		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = VectorC.length;

			if (CoeffBeta != 0.0)
			{
				for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
				{
					VectorC[LongIndex] = CoeffAlpha * VectorA[LongIndex] + CoeffBeta * VectorB[LongIndex];
				}
			}
			else
			{
				for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
				{
					VectorC[LongIndex] = CoeffAlpha * VectorA[LongIndex];
				}
			}
			return;
		}

		// Parallel VectorC = CoeffAlpha * VectorA + CoeffBeta * VectorB
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			if (CoeffBeta != 0.0)
			{
				for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
				{
					VectorC[LongIndex] = CoeffAlpha * VectorA[LongIndex] + CoeffBeta * VectorB[LongIndex];
				}
			}
			else
			{
				for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
				{
					VectorC[LongIndex] = CoeffAlpha * VectorA[LongIndex];
				}
			}
		}
	   );

	} // End LinearCombineVector -- One dimensional double

	public static void LinearCombineVector(double[][] VectorC, double CoeffAlpha, double[][] VectorA, double CoeffBeta, double[][] VectorB)
	{
		int LocalVectorDimension = VectorC[0].GetLength(0);

		if (LocalVectorDimension != VectorA[0].GetLength(0))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LocalVectorDimension)).toString() + " A " + VectorA[0].GetLength(0).toString());
		}

		if (LocalVectorDimension != VectorB[0].GetLength(0))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LocalVectorDimension)).toString() + " B " + VectorB[0].GetLength(0).toString());
		}

		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = VectorC.length;

			if (CoeffBeta != 0.0)
			{
				for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
				{
					for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
					{
						VectorC[LongIndex][LocalVectorIndex] = CoeffAlpha * VectorA[LongIndex][LocalVectorIndex] + CoeffBeta * VectorB[LongIndex][LocalVectorIndex];
					}
				}
			}
			else
			{
				for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
				{
					for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
					{
						VectorC[LongIndex][LocalVectorIndex] = CoeffAlpha * VectorA[LongIndex][LocalVectorIndex];
					}
				}
			}
			return;
		}

		// Parallel VectorC = CoeffAlpha * VectorA + CoeffBeta * VectorB
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			if (CoeffBeta != 0.0)
			{
				for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
				{
					for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
					{
						VectorC[LongIndex][LocalVectorIndex] = CoeffAlpha * VectorA[LongIndex][LocalVectorIndex] + CoeffBeta * VectorB[LongIndex][LocalVectorIndex];
					}
				}
			}
			else
			{
				for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
				{
					for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
					{
						VectorC[LongIndex][LocalVectorIndex] = CoeffAlpha * VectorA[LongIndex][LocalVectorIndex];
					}
				}
			}
		}
	   );

	} // End LinearCombineVector -- TWo dimensional double

	//  Copy TotalSize local vectors from VectorA to VectorC starting in position 0 of VectorC and position StartIndex of VectorA
	public static void CopyVector(double[][] VectorC, double[][] VectorA, int StartIndex, int TotalSize)
	{
		int LocalVectorDimension = VectorC[0].GetLength(0);

		if (LocalVectorDimension != VectorA[0].GetLength(0))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LocalVectorDimension)).toString() + " A " + VectorA[0].GetLength(0).toString());
		}

		if (SALSAUtility.sequentialBLAS)
		{
			for (int LongIndex = 0; LongIndex < TotalSize; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					VectorC[LongIndex + StartIndex][LocalVectorIndex] = VectorA[LongIndex][LocalVectorIndex];
				}
			}
			return;
		}

		// Parallel VectorC = VectorA
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			int endpoint = indexlen + beginpoint;
			if (endpoint > TotalSize)
			{
				endpoint = TotalSize;
			}
				if (ThreadNo == (SALSAUtility.ThreadCount - 1))
				{
					endpoint = TotalSize;
				}
				for (int LongIndex = beginpoint; LongIndex < endpoint; LongIndex++)
				{
					for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
					{
						VectorC[LongIndex + StartIndex][LocalVectorIndex] = VectorA[LongIndex][LocalVectorIndex];
					}
				}
		}
	   );

	} // End CopyVector -- Two dimensional double

	//  Copy TotalSize local vectors from VectorA to VectorC starting in position 0 of VectorC and position StartIndex of VectorA
	public static void CopyVector(String[] VectorC, String[] VectorA, int StartIndex, int TotalSize)
	{
		if (SALSAUtility.sequentialBLAS)
		{
			for (int LongIndex = 0; LongIndex < TotalSize; LongIndex++)
			{
				VectorC[LongIndex + StartIndex] = VectorA[LongIndex];
			}
			return;
		}

		// Parallel VectorC = VectorA
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			int endpoint = indexlen + beginpoint;
			if (endpoint > TotalSize)
			{
				endpoint = TotalSize;
			}
				if (ThreadNo == (SALSAUtility.ThreadCount - 1))
				{
					endpoint = TotalSize;
				}
				for (int LongIndex = beginpoint; LongIndex < endpoint; LongIndex++)
				{
					VectorC[LongIndex + StartIndex] = VectorA[LongIndex];
				}
		}
	   );

	} // End CopyVector -- One dimensional String

	//  Copy TotalSize local vectors from VectorA to VectorC starting in position 0 of VectorC and position StartIndex of VectorA
	public static void CopyVector(double[][][] VectorC, double[][][] VectorA, int StartIndex, int TotalSize)
	{
		int LocalVectorDimension1 = VectorC[0].GetLength(0);

		if (LocalVectorDimension1 != VectorA[0].GetLength(0))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LocalVectorDimension1)).toString() + " A " + VectorA[0].GetLength(0).toString());
		}

		int LocalVectorDimension2 = VectorC[0].GetLength(1);

		if (LocalVectorDimension2 != VectorA[0].GetLength(1))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LocalVectorDimension2)).toString() + " A " + VectorA[0].GetLength(1).toString());
		}

		if (SALSAUtility.sequentialBLAS)
		{
			for (int LongIndex = 0; LongIndex < TotalSize; LongIndex++)
			{
				for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension1; LocalVectorIndex1++)
				{
					for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension2; LocalVectorIndex2++)
					{
						VectorC[LongIndex + StartIndex][LocalVectorIndex1][LocalVectorIndex2] = VectorA[LongIndex][LocalVectorIndex1][LocalVectorIndex2];
					}
				}
			}
			return;
		}

		// Parallel VectorC = VectorA
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			int endpoint = indexlen + beginpoint;
			if (endpoint > TotalSize)
			{
				endpoint = TotalSize;
			}
				if (ThreadNo == (SALSAUtility.ThreadCount - 1))
				{
					endpoint = TotalSize;
				}
				for (int LongIndex = beginpoint; LongIndex < endpoint; LongIndex++)
				{
					for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension1; LocalVectorIndex1++)
					{
						for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension2; LocalVectorIndex2++)
						{
							VectorC[LongIndex + StartIndex][LocalVectorIndex1][LocalVectorIndex2] = VectorA[LongIndex][LocalVectorIndex1][LocalVectorIndex2];
						}
					}
				}
		}
	   );

	} // End CopyVector -- Three dimensional double

	//  Form Vector Dot Product -- One Dimensional
	public static double VectorScalarProduct(double[] VectorA, double[] VectorB)
	{
		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = VectorA.length;

			double Total = 0.0;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				Total += VectorA[LongIndex] * VectorB[LongIndex];
			}
			return Total;
		}


		// Parallel Scalar Product  = Sum(over i) VectorA[i] * VectorB[i]
		GlobalReductions.FindDoubleSum FindScalarProduct = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			double tmp = 0.0;
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				tmp += VectorA[LongIndex] * VectorB[LongIndex];
			}
			FindScalarProduct.addapoint(ThreadNo, tmp);
		}
	   );

		FindScalarProduct.sumoverthreadsandmpi();
		return FindScalarProduct.Total;

	} // End VectorScalarProduct -- One dimensional double

	//  Form Vector Dot Product -- Two Dimensional
	public static double VectorScalarProduct(double[][] VectorA, double[][] VectorB)
	{
		int LocalVectorDimension = VectorB[0].GetLength(0);

		if (LocalVectorDimension != VectorA[0].GetLength(0))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions B" + (new Integer(LocalVectorDimension)).toString() + " A " + VectorA[0].GetLength(0).toString());
		}

		if (SALSAUtility.sequentialBLAS)
		{
			int LongDimension = VectorA.length;

			double Total = 0.0;

			for (int LongIndex = 0; LongIndex < LongDimension; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					Total += VectorB[LongIndex][LocalVectorIndex] * VectorA[LongIndex][LocalVectorIndex];
				}
			}
			return Total;
		}

		// Parallel Scalar Product  = Sum(over i) VectorA[i] * VectorB[i]
		GlobalReductions.FindDoubleSum ScalarProduct = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			double tmp = 0.0;
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
				{
					tmp += VectorB[LongIndex][LocalVectorIndex] * VectorA[LongIndex][LocalVectorIndex];
				}
			}
			ScalarProduct.addapoint(ThreadNo, tmp);
		}
	   );

		ScalarProduct.sumoverthreadsandmpi();
		return ScalarProduct.Total;

	} // End VectorScalarProduct -- Two dimensional double

	public static void CopyMatrix(double[][][][] MatrixC, double[][][][] MatrixA)
	{
		int LocalVectorDimension0 = MatrixC[0][0].GetLength(0);

		if (LocalVectorDimension0 != MatrixA[0][0].GetLength(0))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LocalVectorDimension0)).toString() + " A " + MatrixA[0][0].GetLength(0).toString());
		}

		int LocalVectorDimension1 = MatrixC[0][0].GetLength(1);

		if (LocalVectorDimension1 != MatrixA[0][0].GetLength(1))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LocalVectorDimension1)).toString() + " A " + MatrixA[0][0].GetLength(1).toString());
		}
		int LongDimension0 = MatrixC.length;

		if (LongDimension0 != MatrixA.length)
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LongDimension0)).toString() + " A " + MatrixA.length.toString());
		}
		int LongDimension1 = MatrixC[0].length;

		if (LongDimension1 != MatrixA[0].length)
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions C" + (new Integer(LongDimension1)).toString() + " A " + MatrixA[0].length.toString());
		}

		if (!SALSAUtility.sequentialBLAS)
		{
			RuntimeException e = SALSAUtility.SALSAError("CopyMatrix NOT defined for Decomposed Parameters");

			throw (e);
		}

		for (int LongIndex1 = 0; LongIndex1 < LongDimension0; LongIndex1++)
		{
			for (int LongIndex2 = 0; LongIndex2 < LongDimension1; LongIndex2++)
			{
				for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < LocalVectorDimension0; LocalVectorIndex1++)
				{
					for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < LocalVectorDimension1; LocalVectorIndex2++)
					{
						MatrixC[LongIndex1][LongIndex2][LocalVectorIndex1][LocalVectorIndex2] = MatrixA[LongIndex1][LongIndex2][LocalVectorIndex1][LocalVectorIndex2];
					}
				}
			}
		}

	} // End CopyMatrix(double[,][,] MatrixC, double[,][,] MatrixA)

} // End SALSABLAS
 // End namespace SALSALibrary