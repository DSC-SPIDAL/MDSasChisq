package salsa.mdsaschisq;

// Basic data for vector based parameters
// func not stored as in general too big
public class Desertwind
{

	public double[][][] DiagonalofMatrix; // Diagonal of Matrix as a matrix in VectorDimension for usual Chisq approximation; HALF Second Derivative
	public double[][][] ExactDiagonalofMatrix; // Diagonal of Matrix as a matrix in VectorDimension without usual Chisq approximation; HALF Second Derivative
	public double[][][][] FullMatrix; // Full matrix
	public double[][][][] ExactFullMatrix; // Full matrix
	public double[][] param; // Parameter values (See xshift. These are incremented by MINUS xshift from previous solution)
	public double[][] first; // HALF First Derivatives
	public double[][] xshift; // Minus predicted change in parameter values -- This is Stored with Solution that was made by it
	public double Chisquared; // Value of Chisq
	public int IterationCalculated; // Iteration solution calculated on
	public double[][] func = null; // Chisq residuals

	public Desertwind(int NumberofPoints, int VectorDimension)
	{
		param = new double[NumberofPoints][];
		first = new double[NumberofPoints][];
		xshift = new double[NumberofPoints][];

		if (Hotsun.fullmatrixset)
		{
			FullMatrix = new double[NumberofPoints][NumberofPoints][][];
			ExactFullMatrix = new double[NumberofPoints][NumberofPoints][][];
		}
		else
		{
			DiagonalofMatrix = new double[NumberofPoints][][];
			ExactDiagonalofMatrix = new double[NumberofPoints][][];
		}

		for (int LongIndex = 0; LongIndex < NumberofPoints; LongIndex++)
		{
			if (Hotsun.fullmatrixset)
			{
				for (int OtherIndex = 0; OtherIndex < NumberofPoints; OtherIndex++)
				{
					FullMatrix[LongIndex][OtherIndex] = new double[VectorDimension][VectorDimension];
					ExactFullMatrix[LongIndex][OtherIndex] = new double[VectorDimension][VectorDimension];
				}
			}
			else
			{
				DiagonalofMatrix[LongIndex] = new double[VectorDimension][VectorDimension];
				ExactDiagonalofMatrix[LongIndex] = new double[VectorDimension][VectorDimension];
			}
			param[LongIndex] = new double[VectorDimension];
			first[LongIndex] = new double[VectorDimension];
			xshift[LongIndex] = new double[VectorDimension];
		}

		if (Hotsun.funcset == false)
		{
			func = null;
		}
		Chisquared = 0.0;
	}
} // end Desertwind // end Namespace Manxcat