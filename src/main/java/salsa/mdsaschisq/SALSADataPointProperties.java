package salsa.mdsaschisq;

import java.io.Serializable;

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
        SALSADataPointProperties copy = new SALSADataPointProperties();
        copy.x = this.x;
        copy.y = this.y;
        copy.z = this.z;
        copy.xerr = this.xerr;
        copy.yerr = this.yerr;
        copy.zerr = this.zerr;
        copy.family1 = this.family1;
        copy.family2 = this.family2;
        copy.cluster = this.cluster;
        copy.group = this.group;
        copy.familylabel1 = this.familylabel1;
        copy.familylabel2 = this.familylabel2;
        copy.clusterlabel = this.clusterlabel;
        copy.grouplabel = this.grouplabel;
        copy.pointlabel = this.pointlabel;
        copy.pointlabel = this.pointlabel;
        copy.FixedorVaried = this.FixedorVaried;
        copy.PointType = this.PointType;
        copy.LocalPointNumber = this.LocalPointNumber;
        copy.OriginalPointNumber = this.OriginalPointNumber;
        copy.valuesset = this.valuesset;
        copy.errorsset = this.errorsset;
        copy.source = this.source;
		return copy;
	}

} // End SALSADataPointProperties