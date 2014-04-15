package salsa.mdsaschisq;

public class SALSADataPointProperties implements Serializable
{
	public double x = 0.0; // x value
	public double y = 0.0; // y value
	public double z = 0.0; // z value
	public double xerr = 0.0; // x error value
	public double yerr = 0.0; // y error value
	public double zerr = 0.0; // z error value
	public int family1 = -1; // Index into Family 1 -- all these indices UNDEFINED if negative
	public int family2 = -1; // Index into Family 2
	public int cluster = -1; // Index into Cluster
	public int group = -1; // Index into Undefined Group
	public String familylabel1 = ""; // Label of this family member
	public String familylabel2 = ""; // Label of this family member
	public String clusterlabel = ""; // Label of this cluster member
	public String grouplabel = ""; // Label of this group member
	public String pointlabel = ""; // Label of this point
	public int FixedorVaried = 0; // = 0 Ignored = 1 Varied = 2 Fixed
	public int PointType = 0; // = 0 Normal, > 0 Center
	public int LocalPointNumber = -1; // Point Number in this File
	public int OriginalPointNumber = -1; // Point Number in Original analysis
	public boolean valuesset = false; // If true x y z set in colon section
	public boolean errorsset = false; // If true x y z errors set in colon section
	public String source = ""; // Source run producing  position values

	public final SALSADataPointProperties ShallowCopy()
	{
		return (SALSADataPointProperties)this.clone();
	}

} // End SALSADataPointProperties