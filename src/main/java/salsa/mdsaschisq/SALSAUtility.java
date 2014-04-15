package salsa.mdsaschisq;

import HPC.Utilities.*;
import Manxcat.*;
import Salsa.Core.Configuration.Sections.*;
import MPI.*;

//C# TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if USE_UINT16
//C# TO JAVA CONVERTER NOTE: There is no Java equivalent to C# namespace aliases:
//using TDistance = System.UInt16;
//#elif USE_INT16
//C# TO JAVA CONVERTER NOTE: There is no Java equivalent to C# namespace aliases:
//using System.UInt16 = System.Int16;
//#else
//C# TO JAVA CONVERTER NOTE: There is no Java equivalent to C# namespace aliases:
//using System.Int16 = System.Double;
//#endif

public class SALSAUtility
{
	public static ParallelOptions _parallelOptions;

	public static int checkerboard = 0; // If 0 store full matrix, If 1 lower triangular, If 2 use load balanced checkerboard pattern
	public static MPI.Environment MPI_Environment; // MPI Environment
	public static Intracommunicator MPI_communicator = null; //MPI communicator
	public static int MPI_Rank = 0; // Rank of process
	public static int MPI_Size = 1; // Number of MPI Processes
	public static String ParallelPattern = ""; // Basic Parallel Pattern
	public static String PatternLabel = ""; // One line Label for Output
	public static int DistanceProcessingOption = 0; // Control Processing of Distance on Input if 0 or 1, leave distance unchanged; if =2 square input distance
	public static boolean CalcFixedCrossFixed = true; // If true calculate fixed-fixed terms in Chisq. Must be false if StoredDistanceOption =3;
	public static int chisqcomponent = -1; // -1 ignore, =1 calculate Varied-Varied only, =2 Calculate Varied Fixed, =3 Calculate Fixed Varied, =4 Calculate Fixed Fixed
	public static int StoredDistanceOption = 2; // Specify how distance data stored in memory =1 Illegal, =2 Used by Used (Used = Varied+Fixed), =3 Varied by Used
	public static int DiskDistanceOption = 2; // Specify whats in distance data on disk =1 Original by Original, =2 Used by Used (Used = Varied+Fixed), =3 Varied by Used
	public static double DistanceCut = 0.96; // Ignore if negative. If positive undefine all distances larger than this
	public static int LinkCut = 5; // Delete all points with this # Links <= this (Must have at least 3 for possibly well determined problem)
	public static double AllowedDeletedFraction = 0.25; // Veto run if Deleted Points (from missing distances) greater than this fraction
	public static int TransformMethod = 0; // Define transformation Method if nonzero; 11 Special Blast centered at one
	public static double TransformParameter = 0.125; // Define Parameter to control transformation. It is power of (1-d) for TransformMethod 11
	public static double UndefinedDistanceValue = -1.0; // If positive replace undefined distances by this value; use -1.0 if want to unset these distances
	public static double[] DistanceWeightCuts; // List of Distance Cuts for Weights. These define upper limits of distance bins with "infinity" as upper and 0 of course as lower
	public static int NumberDistanceWeightCuts = 0; // Number of Distance Weight Cuts. This is one less than number of bins
	public static double[] ActualWeightCuts; // List of Weights averaging to 1. Current alkgorithm makes each distance bin have equal total weight so bins with lowest number of points have largest weight

	public static int PointCount_Global = 0; // Total number of used points summed over all threads and processes; same as Hotsun.Number_Vectors
	// This includes varied and fixed points but not ignored points-- it is = NumberFixedPoints + NumberVariedPoints
	// It can include deleted points if some removed 
	public static int PointCount_Process = 0; // Total number of points summed over all threads in this process
	public static int PointCount_Largest = 0; // Largest number of points in all processes
	public static int PointStart_Process = 0; //    First data point in this process
	public static int VariedPointStart_Process = 0; //  First Varied Point (on global count) in Process
	public static int VariedPointCount_Process = 0; //  Number of Varied points in process

	public static int NumberFixedPoints = 0; // Number of Fixed Points
	public static int NumberVariedPoints = 0; // Number of Varied Points
	public static int NumberOriginalPoints = 0; // Total Number of Points specified in Input
	public static int SALSASHIFT = 100; // Allow room for flags in OriginalPointDisposition -- Should be 1 or higher
	public static int[] OriginalPointDisposition; // = 0 Original Point Not Used; = SALSASHIFT + n This is Varied Point n; = -SALSASHIFT -m This is Fixed Point m
	public static int[] FixedPointOriginal; // Value of Original Point Index Corresponding to this fixed point
	public static int[] VariedPointOriginal; // Value of Original Point Index Corresponding to this varied point
	public static int[] UsedPointtoOriginalPointMap; // Value of Original Point Index Corresponding to this Used point
	public static int[] OriginalPointtoUsedPointMap; // Value of Used Point Index Corresponding to this Original point; -1 if OriginalPoint not used
	public static int[] ActualtoNaiveUsedOrder; // Value of original used position of reordered (for load balancing) used points
	public static int[] NaivetoActualUsedOrder; // Value of actual used position for original used position

	public static boolean Usedreordered = false; // If true the used are reordered

	public static SALSADataPointProperties[] GlobalPointProperties; // Properties of USED Point System
	public static SALSAFileProperties GlobalFileProperties; // Properties of Current USED data set

	public static String[] OriginalDataLabels; // Labels of Current ORIGINAL data set
	public static SALSADataPointProperties[] OriginalPointProperties; // Properties of ORIGINAL Point System
	public static SALSAFileProperties OriginalFileProperties; // Properties of Original data set

	//  Within a job data points will be divided into MPI_Size parts -- each part is assigned to a separate MPI Process
	public static int[] PointsperProcess = null; //how many data points each process will take care
	public static int[][] PointsperThreadperProcess = null; // Number of data points in each process-thread

	//	Within a process, data points will be divided into ThreadCount segments( each thread a segment), the array keep the size of each segment
	public static int[] PointsperThread = null; //how many data points each thread will take care
	public static int[] StartPointperThread = null; //the starting point that a thread will take care

	public static int ThreadCount = 1; // maximum number of parallel threads in a process
	public static int NodeCount = 1; // maximum number of separate nodes in run
	public static int MPIperNodeCount = 1; // Number of MPI processes per node
	public static int MPIIOStrategy = 0; // MPI I/O Strategy
	public static boolean sequentialBLAS = false; // If true calculate BLAS sequentially using imput lengths; if false use SALSAParallelism and do using Threads/MPI

	public static double[][] PointDistances; // Point Distances
	public static int MatrixBreakFactor = 1;
	public static int[] ArrayDivision;
	public static int[] Indexsubtraction;
	public static int DivisionSize;

	public static java.util.Date startTime = new java.util.Date(0);
	public static java.util.Date endTime = new java.util.Date(0);

	// Todo. revert to HiPerfTimer if necessary
//        public static HiPerfTimer PreciseTimer;   //    Hold Precise Timing
//        public static HiPerfTimer[] SubTimers;   // Timing Objects
	public static HiPerfTimer PreciseTimer;
	public static HiPerfTimer[] SubTimers;

	public static int NumberofSubTimings = 0; // Number of subtimings
	public static double HPDuration = 0.0; // Time with Precision
	public static double[] SubDurations; // Hold partial timing
	public static String[] SubTimingNames; //  Labels of partial timing
	public static boolean[] SubTimingEnable;
	public static int MPIREDUCETiming = -1;
	public static int MPIREDUCETiming1 = -1;
	public static int MPISENDRECEIVETiming = -1;
	public static int MPIBROADCASTTiming = -1;
	public static int ThreadTiming = -1;

	/* Parameters for density graph creation */
	public static double Xmaxbound = 1.5; // bounding max x value
	public static double Ymaxbound = 1.5; // bounding max y value
	public static int Xres = 50; // resolution of x axis
	public static int Yres = 50; // resolution of y axis
	public static double Alpha = 2;
	public static double Pcutf = 0.85;
	public static boolean Normalize = true;
	public static String ClusterFile = "";
	public static java.util.HashSet<Integer> SelectedClusters; // Selected clusters to plot density graph
	public static boolean IsClustersSelected = false;

	/* Parameters for html creation */
	public static String ManxcatRunName = "";
	public static String ManxcatRunDescription = "";


	//  These are general parameters for C# codes
	public static java.util.ArrayList CosmicOutput = new java.util.ArrayList(1000); // Monitoring Output
	public static boolean ConsoleDebugOutput = true; // If true send Monitoring output to console
	public static int DebugPrintOption = 2; // Control Printing (= 0 None, ==1 Summary, = 2 Full)

	//  Set up Parallel Threading
	public static void SetupParallelOptions()
	{
		SALSAUtility._parallelOptions = new ParallelOptions();
		SALSAUtility._parallelOptions.MaxDegreeOfParallelism = SALSAUtility.ThreadCount;
	}

	public static ParallelOptions getParallelOptions()
	{
		return _parallelOptions;
	}


	public static void SetupDistanceWeights()
	{
		// Set up Distance Weights
		if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(ManxcatCentral.config.DistanceWeightsCuts))
		{
			char[] sep = new char[] {','};
			String[] weightCuts = ManxcatCentral.config.DistanceWeightsCuts.trim().split(java.util.regex.Pattern.quote(sep.toString()), -1);
			NumberDistanceWeightCuts = weightCuts.length;
			DistanceWeightCuts = new double[NumberDistanceWeightCuts];
			ActualWeightCuts = new double[NumberDistanceWeightCuts + 1];
			for (int i = 0; i < NumberDistanceWeightCuts; i++)
			{
				DistanceWeightCuts[i] = Double.parseDouble(weightCuts[i]);
			}
		}
	}

	public static RuntimeException SALSAError(String message)
	{
		System.out.println("SALSA Error " + message);
		RuntimeException e = new RuntimeException(message);
		return e;

	} // end SALSAError

	// PrintOption = 0 Essential Printout
	// PrintOption = 1 Summary Printout
	// PrintOption = 2 Only if full print out requested
	public static void SALSAPrint(int PrintOption, String StufftoPrint)
	{
		if (MPI_Rank != 0)
		{
			return;
		}
		if (DebugPrintOption < PrintOption)
		{
			return;
		}
		CosmicOutput.add(StufftoPrint);

		if (ConsoleDebugOutput)
		{
			System.out.println(StufftoPrint);
		}
		return;
	} // End SALSAPrint

	public static void SALSAGracefulend(String directory, tangible.RefObject<Integer> readinstruction)
	{
		readinstruction.argValue = -1;
		if (SALSAUtility.MPI_Rank != 0)
		{
			return;
		}

		String filename = directory + "\\GracefulEnd.txt";
		if (!(new java.io.File(filename)).isFile())
		{
			return;
		}

		try
		{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//			using (StreamReader sr = File.OpenText(filename))
			StreamReader sr = File.OpenText(filename);
			try
			{
				// Read first line of file
				String inputLineStr;
				if ((inputLineStr = sr.ReadLine()) == null)
				{
					sr.Close();
					return;
				}
				if (inputLineStr.length() > 0)
				{
					inputLineStr = tangible.DotNetToJavaStringHelper.trim(inputLineStr, new char[] {' ', '\t'});
					String[] split = inputLineStr.split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
					if ((split.length > 0) && (split[0].length() > 0))
					{
						readinstruction.argValue = Integer.parseInt(split[0]);
						SALSAPrint(0, " Termination Instruction " + (new Integer(readinstruction.argValue)).toString());
					}
				}
				sr.Close();
			}
			finally
			{
				sr.dispose();
			}
		}
		catch (RuntimeException e)
		{
			SALSAUtility.SALSAError(" Failure reading " + filename + " " + e.toString());
			throw (e);
		}
	}

	public static void SALSAStatus(String directory, String Message)
	{
		if (SALSAUtility.MPI_Rank != 0)
		{
			return;
		}
		String statusfilename = directory + "\\Status.txt";

		try
		{
			StreamWriter sw = null;

			if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(statusfilename))
			{
				sw = new StreamWriter(statusfilename, false, Encoding.UTF8);
			}

			if (sw != null)
			{
				sw.WriteLine(Message);
			}

			sw.Flush();
			sw.Close();
		}
		catch (RuntimeException e)
		{
			System.out.println("Failed writing status data" + e);
			throw (e);
		}
	} // End SALSAStatus

	// todo: saliya - Possible Improvement: handle comments in this function to use the List instead of the string
	//                you will have to search for the places where comment is manipulated in other
	//                places as well
	public static void SALSAUpdateMetaData(ManxcatSection Configuration)
	{
		String Comment = Configuration.Comment;
		for (int i = 0; i < CosmicOutput.size(); i++)
		{
			if (i == 0)
			{
				continue;
			}
			Comment += "\n" + CosmicOutput.get(i);
		}
		Comment = Comment.replace(':', '-');
		Configuration.Comment = Comment;
	} // end SALSAUpdateMetaData(string Comment)

	public static void InitializeTiming(int InputNumberofTimers)
	{
		NumberofSubTimings = InputNumberofTimers;
		SubTimers = new HiPerfTimer[NumberofSubTimings]; // Timing Objects

		SubDurations = new double[NumberofSubTimings]; // Hold partial timing
		SubTimingEnable = new boolean[NumberofSubTimings];
		SubTimingNames = new String[NumberofSubTimings];

		for (int itimer = 0; itimer < NumberofSubTimings; itimer++)
		{
			SubTimers[itimer] = new HiPerfTimer();
			SubDurations[itimer] = 0.0;
			SubTimingEnable[itimer] = true;
		}

		PreciseTimer = new HiPerfTimer();
		PreciseTimer.Start();

		startTime = new java.util.Date(); // Set the initial time.
	} // End InitializeTiming

	public static void SetUpSubTimer(int TimingIndex, String TimingLabel)
	{
		if (TimingIndex >= NumberofSubTimings)
		{
			SALSAError("Error in Timing Index " + (new Integer(TimingIndex)).toString() + " Max " + (new Integer(NumberofSubTimings)).toString());
			return;
		}
		SubTimingNames[TimingIndex] = TimingLabel;
		SubTimingEnable[TimingIndex] = true;
		SubDurations[TimingIndex] = 0.0;
	} // End SetUpSubTimer

	public static void SetUpMPISubTimers(int StartTimingIndex, String MPISetLabel)
	{
		if ((StartTimingIndex + 4) >= NumberofSubTimings)
		{
			SALSAError("Error in  MPI Timing Index " + (new Integer(StartTimingIndex + 2)).toString() + " Max " + (new Integer(NumberofSubTimings)).toString());
			return;
		}
		SetUpSubTimer(StartTimingIndex, MPISetLabel + "MPI Reduce");
		SetUpSubTimer(StartTimingIndex + 1, MPISetLabel + "MPI SR");
		SetUpSubTimer(StartTimingIndex + 2, MPISetLabel + "MPI Bcast");
		SetUpSubTimer(StartTimingIndex + 3, MPISetLabel + "MPI Global Reductions");
		SetUpSubTimer(StartTimingIndex + 4, MPISetLabel + "Thread Global Reductions");

		if ((!MPISetLabel.equals("")) && (!MPISetLabel.equals("Lib ")))
		{
			return;
		}
		MPIREDUCETiming = StartTimingIndex;
		MPISENDRECEIVETiming = StartTimingIndex + 1;
		MPIBROADCASTTiming = StartTimingIndex + 2;
		MPIREDUCETiming1 = StartTimingIndex + 3;
		ThreadTiming = StartTimingIndex + 4;

	} // End SetUpMPISubTimers

	public static void InterimTiming()
	{
		PreciseTimer.Stop();
		// Todo. revert to HiPerfTimer if necessary
//            HPDuration += PreciseTimer.Duration;
		HPDuration += PreciseTimer.getDuration() * 0.001;
		PreciseTimer.Start();
	} // end EndTiming

	public static void EndTiming()
	{
		endTime = new java.util.Date(); // Set the final time.
		PreciseTimer.Stop();
		HPDuration += PreciseTimer.getDuration() * 0.001;
	} // end EndTiming

	public static void StartSubTimer(int TimingIndex)
	{
		if (TimingIndex < 0)
		{
			return;
		}

		if (SubTimingEnable[TimingIndex])
		{
			SubTimers[TimingIndex].Start();
		}
	} // End StartSubTimer

	public static double StopSubTimer(int TimingIndex)
	{
		if (TimingIndex < 0)
		{
			return 0.0;
		}

		if (SubTimingEnable[TimingIndex])
		{
			SubTimers[TimingIndex].Stop();
			double tmp = SubTimers[TimingIndex].getDuration() * 0.001;
			SubDurations[TimingIndex] += tmp;
			return tmp;
		}
		return 0.0;

	} // End StopSubTimer


	// Synchronize a boolean across MPI processes
	public static void SynchronizeMPIvariable(tangible.RefObject<Boolean> cosmicboolean)
	{
		if (SALSAUtility.MPI_Size > 1)
		{
			SALSAUtility.StartSubTimer(SALSAUtility.MPIBROADCASTTiming);
			SALSAUtility.MPI_communicator.<Boolean>Broadcast(cosmicboolean, 0);
			SALSAUtility.StopSubTimer(SALSAUtility.MPIBROADCASTTiming);
		}
		return;
	} // End synchronizeboolean(bool cosmicboolean)

	// Synchronize a Double across MPI processes
	public static void SynchronizeMPIvariable(tangible.RefObject<Double> cosmicdouble)
	{
		if (SALSAUtility.MPI_Size > 1)
		{
			SALSAUtility.StartSubTimer(SALSAUtility.MPIBROADCASTTiming);
			SALSAUtility.MPI_communicator.<Double>Broadcast(cosmicdouble, 0);
			SALSAUtility.StopSubTimer(SALSAUtility.MPIBROADCASTTiming);
		}
		return;
	} // End synchronizeboolean(double cosmicdouble)
	public static void SynchronizeMPIvariable(tangible.RefObject<Integer> cosmicint)
	{
		if (SALSAUtility.MPI_Size > 1)
		{
			SALSAUtility.StartSubTimer(SALSAUtility.MPIBROADCASTTiming);
			SALSAUtility.MPI_communicator.<Integer>Broadcast(cosmicint, 0);
			SALSAUtility.StopSubTimer(SALSAUtility.MPIBROADCASTTiming);
		}
		return;
	} // End synchronizeboolean(int cosmicint)

	public static MPIReducePlusIndex MinwithIndex(MPIReducePlusIndex one, MPIReducePlusIndex two)
	{
		if (one.index < 0)
		{
			return two;
		}
		if (two.index < 0)
		{
			return one;
		}
		if (one.value < two.value)
		{
			return one;
		}
		else
		{
			return two;
		}
	}

	public static MPIReducePlusIndex MaxwithIndex(MPIReducePlusIndex one, MPIReducePlusIndex two)
	{
		if (one.index < 0)
		{
			return two;
		}
		if (two.index < 0)
		{
			return one;
		}
		if (one.value > two.value)
		{
			return one;
		}
		else
		{
			return two;
		}
	}

	// Return Value and Index corresponding to minimum over all MPI Processes
	//  Replace Input values by minimum values
	public static void AllReduceMinWithIndex(tangible.RefObject<Double> ProcessValue, tangible.RefObject<Integer> ProcessIndex)
	{
		if (SALSAUtility.MPI_Size > 1)
		{
			MPIReducePlusIndex LocalStructure = new MPIReducePlusIndex(ProcessIndex.argValue, ProcessValue.argValue);
			MPIReducePlusIndex TotalStructure = SALSAUtility.MPI_communicator.<MPIReducePlusIndex>Allreduce(LocalStructure, MinwithIndex);
			ProcessValue.argValue = TotalStructure.value;
			ProcessIndex.argValue = TotalStructure.index;
		}
	} // End AllReduceMinWithIndex

	// Return Value and Index corresponding to maximum over all MPI Processes
	//  Replace Input values by maximum values
	public static void AllReduceMaxWithIndex(tangible.RefObject<Double> ProcessValue, tangible.RefObject<Integer> ProcessIndex)
	{
		if (SALSAUtility.MPI_Size > 1)
		{
			MPIReducePlusIndex LocalStructure = new MPIReducePlusIndex(ProcessIndex.argValue, ProcessValue.argValue);
			MPIReducePlusIndex TotalStructure = SALSAUtility.MPI_communicator.<MPIReducePlusIndex>Allreduce(LocalStructure, MaxwithIndex);
			ProcessValue.argValue = TotalStructure.value;
			ProcessIndex.argValue = TotalStructure.index;
		}
	} // End AllReduceMaxWithIndex

	// Write performance results of process 0 into a file
	public static void WriteTiming_Cluster(String fileName, String RunLabel, int RunNumber, String DataFileName, String ProcessorName)
	{
		try
		{
			// Print results
			StreamWriter sw = null;

			// New file flag
			boolean newfile = false;

			if (!(new java.io.File(fileName)).isFile())
			{
				newfile = true;
			}

			if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(fileName))
			{
				sw = new StreamWriter(fileName, true, Encoding.UTF8); //Append

				if (newfile)
				{
					// Formats must match those in following data write
					sw.Write(String.format("%1$-15s%2$-15s%3$-8s%4$-12s%5$-8s%6$-12s%7$-12s%8$-16s%9$-12s%10$-25s%11$-25s%12$-25s", "Duration(ms)", "HPDuration(ms)", "Thread#", "MPIperNode", "Node#", "Pt#/Process", "Pt#/Global", "Run Label", "Run Number", "DataFileName", "CurrentTime", "ProcessorName"));

					for (int subtimer = 0; subtimer < SALSAUtility.NumberofSubTimings; subtimer++)
					{
						if (!SALSAUtility.SubTimingEnable[subtimer])
						{
							continue;
						}

						sw.Write(String.format("%1$-15s%2$-8s", SALSAUtility.SubTimingNames[subtimer], "Ratio")); // Subset Time    Name
					}
					sw.WriteLine(" "); // Blank line
					sw.WriteLine(" "); // Blank line
				}
			}

			if (sw != null)
			{
				TimeSpan duration = endTime - startTime;

				// Format syntax (Variable Number , Spaces with - left justified and no dash right aligned
				sw.Write(String.format("%1$-15s%2$-15s%3$-8s%4$-12s%5$-8s%6$-12s%7$-12s%8$-16s%9$-12s%10$-25s%11$-25s%12$-25s", Math.round(duration.TotalMilliseconds), Math.round(SALSAUtility.HPDuration * Math.pow(10, 0)) / Math.pow(10, 0), SALSAUtility.ThreadCount, SALSAUtility.MPIperNodeCount, SALSAUtility.NodeCount, SALSAUtility.PointCount_Process, SALSAUtility.PointCount_Global, RunLabel, RunNumber, DataFileName, SALSAUtility.startTime.ToLocalTime(), ProcessorName)); // Processor Used -  Current Time -  Name of input data file -  Run Number -  Label for Run -  Global pointsn -  Local points -  Node# aquired -  Process# per node -  Thread# per process - High performance timer -  Total time
				// Todo. revert to HiPerfTimer if necessary
//                    Math.Round(SALSAUtility.HPDuration * 0.001, 0),                                                                                                        //High performance timer

				for (int subtimer = 0; subtimer < SALSAUtility.NumberofSubTimings; subtimer++)
				{
					if (!SALSAUtility.SubTimingEnable[subtimer])
					{
						continue;
					}

					// Todo. revert to HiPerfTimer if necessary
//                        double SubTime = SALSAUtility.SubDurations[subtimer] * 0.001;
					double SubTime = SALSAUtility.SubDurations[subtimer];

//                        double SubRatio = SubTime / (HPDuration * 0.001);
					double SubRatio = SubTime / (HPDuration);
					sw.Write(String.format("%1$-15s%2$-8s", Math.round(SubTime), String.format("%0.4f", SubRatio))); // Subset Time
				}
				sw.WriteLine(" "); // Blank line
			}

			sw.Flush();
			sw.Close();
		}
		catch (RuntimeException e)
		{
			System.out.println("Failed writing data" + e);
			throw (e);
		}
	}
} // End class SalsaUtility