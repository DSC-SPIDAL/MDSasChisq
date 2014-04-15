package salsa.mdsaschisq;

import HPC.Utilities.*;
import Manxcat.*;
import Salsa.Core.Configuration.Sections.*;
import MPI.*;

public class MPIReducePlusIndex implements Serializable
{
	public int index;
	public double value;

	public MPIReducePlusIndex(int _index, double _value)
	{
		index = _index;
		value = _value;
	}
} // End MPIMinPlusIndex
 // end Namespace SALSALibrary