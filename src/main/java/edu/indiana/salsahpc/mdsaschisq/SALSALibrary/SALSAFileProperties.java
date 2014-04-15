package SALSALibrary;

//    OLD Comments on input Data
//    public string DataFileName = " "; // Input data file name
//    public string DataLabelsFileName = "";   // File with application specific labels for this data
//    public string ReducedVectorOutputFileName = "";  // File in ResultDirectoryName holding ReducedVector Output
//    public string ClusterDirectory = "";    // Directory to find Cluster Output
//    public int DataPoints = 0;  // Number of Original Data Points
//    public int ProcessingOption = 0; // Control Requested Processing
//    public int DistanceProcessingOption = 0; // Control Input data Processing ( = 2 SQUARE Input Values)
//    public int InitializationOption = 0;    // Control Initialization =0 None =1 Old Model =2 New File Type
//    public string InitializationFileName = "";  // File
//    public bool write2Das3D = true; // for 2D MDS write output for 3D viewer
//    public int LocalVectorDimension = 3;    // Final MDS dimension (2 or 3 expected)
//    public string selectedvariedpoints = "";  // List of varied points specified by point number or cluster number 
//    public string VariedPointCriterion = "all"; // If "all" or "rest" all the non-fixed points varied; if "allinfile" choose all in file and ignore selectedvariedpoints property
//    // if "originalpoint" select by Original Point; = "family1", "family2", "cluster", "group" select by this division
//    public string selectedvariedpointfile = "";  // File containing mapping of groups to original values 
//    public string selectedfixedpoints = "";  // List of fixed points specified by point number or group number 
//    public string FixedPointCriterion = "none"; // If "none" there are no fixed points; if "allinfile" choose all in file and ignore selectedfixedpoints property
//    // if "originalpoint" select by Original Point; = "family1", "family2", "cluster", "group" select by this division
//    public string selectedfixedpointfile = "";  // File containing coordinates of fixed points

//    public int RotationOption = 0;  // Control Rotation 
//    public string RotationLabelsFileName = "";   // File with application specific labels for data to be rotated

//    public int Chisqnorm = 0;   // Method of Normalizing Chisq
//    public int DistanceFormula = 0;   // Method of Calculating Distance
//    public int FullSecondDerivativeOption = 0; // Use Full second Derivative Calculation if > 0
//    public double MinimumDistance = -0.01;   // Minimum Distance for Chisq Normalization
//    public double FunctionErrorCalcMultiplier = 10.0; // Multiplicative Constant for Function (square this for Chisq) Error Calculation
//    public double ChisqPrintConstant = 1.0; // Multiplicative Constant for Chisq Printing
//    public int maxit = 60;  // Maximum Number of Manxcat Iterations
//    public int nbadgo = 6;  // Limit on Iterations with no progress
//    public double ChisqChangePerPoint = 0.001;  // Stop when this is Chisq change per point
//    public double FletcherRho = 0.25;   // Fletcher's parameter rho
//    public double FletcherSigma = 0.75; // Fletcher's Parameter sigma
//    public double Omega = 1.25; //  Parameter Omega deciding if to do Linear Search if OmegaOption = 1
//    public int OmegaOption = 0; // If nonzero do linear searches if promising
//    public double QHighInitialFactor = 0.01; // Initial call uses this fraction of Qhigh for Q (added 2009)
//    public double QgoodReductionFactor = 0.5;    // Set to zero to make Q = 0. for good solutions
//    public double addonforQcomputation = 2.0;   // Add this to make Power Method Safer
//    public double extraprecision = 0.05; //  Extra precision for Qlow
//    public int QLimitscalculationInterval = 1; //  recalculate Qhigh Qlow Trace and Norm after this number of iterations
//    public int InitialSteepestDescents = 0;  // Start minimization with this number of Steepest Descent Iterations
//    public double TimeCutmillisec = -1.0;  // if positive, Stop when time used in milliseconds greater than timmax
//    public bool derivtest = false; // If true, do a derivative 

//    public double CGResidualLimit = 0.00001;    // Stop when Norm of Residual is less than this
//    public int PowerIterationLimit = 200;    // Limit on Power Iterations
//    public double eigenvaluechange = 0.001; // convergence test in power method for eigenvalues
//    public double eigenvectorchange = 0.001; // convergence test in power method for eigenvector

//    public string TimingOutputFileName = " "; // Timing Output file name
//    public string BaseResultDirectoryName = "";  // Base Directory holding all Runs
//    public string ResultDirectoryExtension = ""; // 
//    public string ControlDirectoryName = ""; // Directory holding Control Information for this RunSet
//    public string RunSetLabel = "";  //  Label for RunSet
//    public int RunNumber = 0;    // Unique number of run

//    public string pattern = ""; // Parallel Pattern
//    public int ThreadCount = 1;    // Number of Threads
//    public int NodeCount = 1;    // Number of Nodes
//    public int MPIperNodeCount = 1;    // Number of MPI Processes per Node
//    public int MPIIOStrategy = 0;   // Strategy for I/O in MPI

//    public int HistogramBinCount = 100; // Bin Count for Histograms
//    public string Comment = ""; // User Comment
//    public string Extradata1 = "";  // Use in special ways
//    public string Extradata2 = "";  // Use in special ways
//    public string Extradata3 = "";  // Use in special ways
//    public string Extradata4 = "";  // Use in special ways
//    public int ExtraOption1 = 0;    // use in various ways
//    public int DebugPrintOption = 1; // Control Debug Printing (= 0 None, = 1 Full, ==2 Summary)
//    public bool ConsoleDebugOutput = false; // If true send debug output to Console

public class SALSAFileProperties implements Serializable
{
	public int LocalVectorDimension = 3; // Vector dimension of Mapped Points
	public int ClusterStartIndex = 0; // Cluster Indices start at 0 or 1
	public int OriginalPointStartIndex = 0; // Original Point Indices start at 0 or 1
	public int LocalPointStartIndex = 0; // File Point Indices start at 0 or 1
	public int FileGenerationType = 0; // 0 Clustering; 1 Grouping; 2 MDS Original; 3 MDS Incremental; 4 MDS Fixed; 5 MDS Rotation
	public String FamilyName1 = ""; // Name of Family 1
	public String FamilyName2 = ""; // Name of Family 2
	public String GroupName = "Undefined Group"; // Name of Group
	public String ClusterName = ""; // Name of Cluster -- typically Run Label
	public String Comment = ""; // General Comment
	public int NumberOriginalPoints = 0; // Number of Original Points
	public int NumberPointsinFile = 0; // Number of Points Written
	public String RotationParameters = ""; // List of Rotation parameters separated by commas
	public int NumberRotationParameters = 0; // Number of Rotation Parameters

	public final SALSAFileProperties ShallowCopy()
	{
		return (SALSAFileProperties) this.clone();
	}

} // End SALSAFileProperties