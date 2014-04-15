package Manxcat;

import MPI.*;
import SALSALibrary.*;
import MathNet.Numerics.LinearAlgebra.*;

public class ChisqFirstandSecond implements Serializable
{
	public double chisq;
	public double[] first;
	public double[] second;
	public int VectorSize;
	public int SecondSize;

	public ChisqFirstandSecond(int NumberParms)
	{
		VectorSize = NumberParms;
		first = new double[NumberParms];
		SecondSize = (NumberParms * (1 + NumberParms)) / 2;
		second = new double[SecondSize];
		chisq = 0.0;
		for (int iparm = 0; iparm < NumberParms; iparm++)
		{
			first[iparm] = 0.0;
		}
		for (int iparm = 0; iparm < SecondSize; iparm++)
		{
			second[iparm] = 0.0;
		}
	}

	public static void Zeromember(tangible.RefObject<ChisqFirstandSecond> Instance)
	{
		Instance.argValue.chisq = 0.0;
		for (int iparm = 0; iparm < Instance.argValue.VectorSize; iparm++)
		{
			Instance.argValue.first[iparm] = 0.0;
		}
		for (int iparm = 0; iparm < Instance.argValue.SecondSize; iparm++)
		{
			Instance.argValue.second[iparm] = 0.0;
		}
	}
}