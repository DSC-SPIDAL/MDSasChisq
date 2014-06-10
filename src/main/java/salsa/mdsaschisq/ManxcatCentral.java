package salsa.mdsaschisq;

import com.google.common.base.Optional;
import com.google.common.base.Strings;
import mpi.MPI;
import mpi.MPIException;
import org.apache.commons.cli.*;
import salsa.configuration.ConfigurationMgr;
import salsa.configuration.sections.MDSasChisqSection;

import java.io.File;
import java.io.IOException;
import java.nio.ByteOrder;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.regex.Pattern;

public class ManxcatCentral
{
    private static Options programOptions = new Options();
    static {
        programOptions.addOption(String.valueOf(ManxcatConstants.CMD_OPTION_SHORT_C),ManxcatConstants.CMD_OPTION_LONG_C, true,
                ManxcatConstants.CMD_OPTION_DESCRIPTION_C);
        programOptions.addOption(String.valueOf(ManxcatConstants.CMD_OPTION_SHORT_N),ManxcatConstants.CMD_OPTION_LONG_N,true,
                ManxcatConstants.CMD_OPTION_DESCRIPTION_N);
        programOptions.addOption(String.valueOf(ManxcatConstants.CMD_OPTION_SHORT_T),ManxcatConstants.CMD_OPTION_LONG_T,true,
                ManxcatConstants.CMD_OPTION_DESCRIPTION_T);
    }
    private static Format dateFormatter = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
	//  Input Variables               
	public static MDSasChisqSection config;
	public static String ConfigurationFileName = ""; // Configuration file name
	public static int ProcessingOption = 0; // Processing Option

	//  Utility variables
	public static MPI2DDoubleVectorPacket TogoDistributed2DDoubleVector; // Parameter Data to be sent on MPI transfers
	public static MPI2DDoubleVectorPacket FromAfar2DDoubleVector; // Data received by MPI transfers
	public static MPI1DStringVectorPacket TogoDistributed1DStringVector; // Parameter Data to be sent on MPI transfers
	public static MPI1DStringVectorPacket FromAfar1DStringVector; // Data received by MPI transfers

	public static MPI2DDoubleVectorPacket TogoDiagVector; // diag Data to be sent on MPI transfers
	public static MPI2DDoubleVectorPacket TogoSqDgInvVector; // sqdiag Data to be sent on MPI transfers

	public static boolean violat; // error return for Calcfg (= false i,mplies parameter values changed as violated boundary conditions)
	public static double ChisqFunctionCalcMultiplier = 10.0; // Multiplicative Constant for Function (square this for Chisq) Calculation)
	public static double ChisqPrintConstant = 1.0; // Extrs Multiplicative Constant for Chisq Printing (divided by number of degrees of freedom)
	public static String ResultDirectoryName = ""; // Full String for Result Directory
	public static String ActualOutputFileName = ""; // Full Output File namee

	//  Note Best solution alwasys stored with Scaling or Addition of Marquardt Factor if this requested
	public static boolean BestSolutionSet = false; // If true there is a best solutiom
	public static boolean BestxshiftSet = false; // If true there is a best xshift set
	public static boolean GlobalParameterSet = false; // If True Global parameters set for Current Solution
	public static int IterationtorecalculateQLimits = 0; // Iteration after and including which calculating Q limits warranted
	public static double LineFactorGuess = 1.0; // Initial Guess at Line Factor

    static void main(String[] args) {
        Optional<CommandLine> parserResult = parseCommandLineArguments(args, programOptions);
        if (!parserResult.isPresent()) {
            System.out.println(ManxcatConstants.ERR_PROGRAM_ARGUMENTS_PARSING_FAILED);
            new HelpFormatter().printHelp(ManxcatConstants.PROGRAM_NAME, programOptions);
            return;
        }

        CommandLine cmd = parserResult.get();
        if (!(cmd.hasOption(ManxcatConstants.CMD_OPTION_LONG_C) &&
                cmd.hasOption(ManxcatConstants.CMD_OPTION_LONG_N) &&
                cmd.hasOption(ManxcatConstants.CMD_OPTION_LONG_T))) {
            System.out.println(ManxcatConstants.ERR_INVALID_PROGRAM_ARGUMENTS);
            new HelpFormatter().printHelp(ManxcatConstants.PROGRAM_NAME, programOptions);
            return;
        }

        //  Read Metadata using this as source of other metadata
        ReadConfiguration(cmd);

        try {
            //  Set up Parallelism
            SALSAParallelism.SetupParallelism(args);
        } catch (MPIException e) {
            SALSAUtility.printException(e);
            return; // End program on error
        }

        // Define Points to be used
        SALSA_ProcessVariedandFixed.setupFixedandVaried();

        // Set up Decomposition of USED points
        SALSAParallelism.SetParallelDecomposition();

        //  Redistribute points so an equal number (upto 1) varied in each thread
        SALSA_ProcessVariedandFixed.RedistributePoints();

        try {
            // Set up Normalizations
            // This can be customized in set up routines
            ChisqPrintConstant = config.ChisqPrintConstant;
            ChisqFunctionCalcMultiplier = config.FunctionErrorCalcMultiplier;

            // Determines type of fit
            ProcessingOption = config.ProcessingOption;

            //  Set up Hotsum and DoubleStar initial parameters
            if (ProcessingOption <= 100) {
                ManxcatMDS.SetupHotsunforMDS();
            } else if (ProcessingOption <= 200) {
                RotateManxcatMDS.SetupHotsunforRotateMDS();
            }
            Hotsun.SetupManxcat();

            if (Hotsun.DecomposeParameters) { // Set up MPI if parallel parameter
                FromAfar2DDoubleVector = setup2DDoubleMPIPacket();
                TogoDistributed2DDoubleVector = setup2DDoubleMPIPacket();
                FromAfar1DStringVector = setup1DStringMPIPacket();
                TogoDistributed1DStringVector = setup1DStringMPIPacket();
                TogoDiagVector = setup2DDoubleMPIPacket();
                TogoSqDgInvVector = setup2DDoubleMPIPacket();
            }

            //  Set up Timing
            SALSAUtility.InitializeTiming(8);
            SALSAUtility.SetUpMPISubTimers(3, "");
            SALSAUtility.SetUpSubTimer(0, "MatrixSolve");
            SALSAUtility.SetUpSubTimer(1, "EigenSolve");
            SALSAUtility.SetUpSubTimer(2, "Calcfg");
            Hotsun.TotalTimeUsed = 0.0;

            // Set up general Manxcat Application structure
            Hotsun.FullSecondDerivative = false;

            //  Increase Run Number
            if ((ProcessingOption == 0) || (ProcessingOption == 100)) {
                ++ManxcatCentral.config.RunNumber;
            }

            //  Construct One Line Label
            config.Pattern = SALSAUtility.ParallelPattern;
            SALSAUtility.PatternLabel = String
                    .format("==== MDS %1$s ==== Option:%2$s PointCount_Global:%3$s ==== Run %4$s %5$s == File %6$s ==" +
                                    " %7$s ==== ",
                            SALSAUtility.ParallelPattern, ManxcatCentral.config.ProcessingOption,
                            SALSAUtility.PointCount_Global, ManxcatCentral.config.RunSetLabel,
                            ManxcatCentral.config.RunNumber, ManxcatCentral.config.DistanceMatrixFile,
                            SALSAUtility.startTime);
            SALSAUtility.SALSAPrint(0, SALSAUtility.PatternLabel);

            ManxcatCentral.config.Comment += (!Strings
                    .isNullOrEmpty(ManxcatCentral.config.Comment) ? "\n" : "") + SALSAUtility.PatternLabel;

            // Initial Processing Complete
            // Make certain all processes have processed original data before writing updated
            SALSAUtility.MPI_communicator.Barrier();

            // Set results directory
            String timestamp = dateFormatter.format(new Date());
            String pattern = String.format("%1$sx%2$sx%3$s", SALSAUtility.ThreadCount, SALSAUtility.MPIperNodeCount,
                                           SALSAUtility.NodeCount);

            ManxcatCentral.config.ResultDirectoryExtension = String
                    .format("%1s %2$s %3$s %4$s", ManxcatConstants.PROGRAM_NAME, ManxcatCentral.config.RunSetLabel,
                            timestamp, pattern);
            ManxcatCentral.ResultDirectoryName = ManxcatCentral.config.BaseResultDirectoryName + File
                    .pathSeparatorChar + ManxcatCentral.config.ResultDirectoryExtension;

            if (SALSAUtility.MPIIOStrategy > 0) {
                ManxcatCentral.ResultDirectoryName += "-Unit" + SALSAUtility.MPI_Size;
            }

            if ((SALSAUtility.MPIIOStrategy > 0) || (SALSAUtility.MPI_Rank == 0)) {
                if ((new java.io.File(ResultDirectoryName)).isDirectory()) {
                    SALSAUtility.SALSAPrint(0, "The directory " + ResultDirectoryName + " exists");
                } else {
                    File dir = new File(ResultDirectoryName);
                    boolean success = dir.mkdir();
                    SALSAUtility.SALSAPrint(0,
                                            "The directory " + ResultDirectoryName + " creation was " + (success ?
                                                    "successful" : "unsuccessful") + " at " + dateFormatter
                                                    .format(new Date(dir.lastModified())));
                }
            }

            // Make certain all processes have processed original data before writing updated
            SALSAUtility.MPI_communicator.Barrier();

            //  Write out Current Summary and Updated Control File
            if ((SALSAUtility.MPIIOStrategy > 0) || (SALSAUtility.MPI_Rank == 0)) {
                config.saveAs(config.SummaryOutputFileName);
            }


            // Make certain all processes have processed original data before writing updated
            SALSAUtility.MPI_communicator.Barrier();

            if (ProcessingOption <= 90) {
                /* Actually run application */
                ManxcatCentral.ManxcatControl(ManxcatMDS::Calcfg,
                                              ManxcatMDS::SetupMDSasChisq,
                                              ManxcatMDS::InitializeParameters,
                                              ManxcatMDS::SolveMatrix,
                                              ManxcatMDSWriteSolution::WriteMDS,
                                              ManxcatMDS::FindQlimits,
                                              ManxcatMDS::GlobalMatrixVectorProduct,
                                              ManxcatMDS::Sequel);

            } else if (ProcessingOption <= 100) {
                /* Special case of <= 100 where doing ManxcatAsChisq for real is omitted.
				 * Original distances are read in from file along with previously
				 * calculated coordinates specified by initialization file and density
				 * graphs along with html page are created */
                ManxcatMDS.SetupMDSasChisq();
                // fill hotsun.globalparam with initialization file
                ManxcatMDS.FillupHotsun();
                ManxcatMDS.Sequel();
            } else if (ProcessingOption <= 200) {
                ManxcatCentral.ManxcatControl(RotateManxcatMDS::Calcfg, RotateManxcatMDS::SetupRotateMDS,
                                              RotateManxcatMDS::InitializeParameters, GenericManxcat::SolveMatrix,
                                              RotateManxcatMDS::WriteRotateMDS, GenericManxcat::FindQlimits,
                                              GenericManxcat::GlobalMatrixVectorProduct, RotateManxcatMDS::Sequel);
            }
            SALSAUtility.MPI_communicator.Barrier();

            String TimingMessage = "\nTiming ";
            for (int subtimer = 0; subtimer < SALSAUtility.NumberofSubTimings; subtimer++) {
                if (!SALSAUtility.SubTimingEnable[subtimer]) {
                    continue;
                }
                double SubTime = SALSAUtility.SubDurations[subtimer];
                TimingMessage += " " + SALSAUtility.SubTimingNames[subtimer] + " " + String
                        .format("%.0f", SubTime) + " (" + String
                        .format("%.3f", SubTime / SALSAUtility.HPDuration) + ")";
            }
            SALSAUtility.SALSAPrint(0, TimingMessage);

            if ((SALSAUtility.MPIIOStrategy > 0) || (SALSAUtility.MPI_Rank == 0)) {
                MetaDataIO.WriteResults_Cluster(config.SummaryOutputFileName, SALSAUtility.CosmicOutput);
                SALSAUtility.WriteTiming_Cluster(config.TimingOutputFileName, config.RunSetLabel, config.RunNumber,
                                                 config.DistanceMatrixFile,
                                                 MPI.Get_processor_name());
                ManxcatHtmlUtility.WriteHTML();
            }
        } catch (IOException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        } finally {
            SALSAParallelism.TearDownParallelism();
        }
    }

    /**
     * Parse command line arguments
     * @param args Command line arguments
     * @param opts Command line options
     * @return An <code>Optional&lt;CommandLine&gt;</code> object
     */
    private static Optional<CommandLine> parseCommandLineArguments(String [] args, Options opts){

        CommandLineParser optParser = new GnuParser();

        try {
            return Optional.fromNullable(optParser.parse(opts, args));
        } catch (ParseException e) {
            System.out.println(e);
        }
        return Optional.fromNullable(null);
    }

	private static void ReadConfiguration(CommandLine cmd)
	{
        config = ConfigurationMgr.LoadConfiguration(cmd.getOptionValue(ManxcatConstants.CMD_OPTION_LONG_C)).mdSasChisqSection;
        SALSAUtility.NodeCount = Integer.parseInt(cmd.getOptionValue(ManxcatConstants.CMD_OPTION_LONG_N));
        SALSAUtility.ThreadCount = Integer.parseInt(cmd.getOptionValue(ManxcatConstants.CMD_OPTION_LONG_T));
        // Override section's node and thread counts with command line values if different
        if (config.getNodeCount() != SALSAUtility.NodeCount) {
            config.setNodeCount(SALSAUtility.NodeCount);
        }
        if (config.getThreadCount() != SALSAUtility.ThreadCount) {
            config.setThreadCount(SALSAUtility.ThreadCount);
        }

		SALSAUtility.NumberOriginalPoints = config.DataPoints;
		SALSAUtility.CalcFixedCrossFixed = config.CalcFixedCrossFixed;
		SALSAUtility.StoredDistanceOption = config.StoredDistanceOption;
		SALSAUtility.DiskDistanceOption = config.DiskDistanceOption;
		SALSAUtility.UndefinedDistanceValue = config.UndefinedDistanceValue;

		SALSAUtility.DistanceCut = config.DistanceCut;
		SALSAUtility.LinkCut = config.LinkCut;
		SALSAUtility.TransformMethod = config.TransformMethod;
		SALSAUtility.TransformParameter = config.TransformParameter;

		SALSAUtility.NodeCount = config.NodeCount; // Number of Nodes

		SALSAUtility.Xmaxbound = config.XmaxBound > 0 ? config.XmaxBound : 1.5;
		SALSAUtility.Ymaxbound = config.YmaxBound > 0 ? config.YmaxBound : 1.5;
		SALSAUtility.Xres = config.Xres > 0 ? config.Xres : 50;
		SALSAUtility.Yres = config.Yres > 0 ? config.Yres : 50;
		SALSAUtility.Alpha = config.Alpha > 0 ? config.Alpha : 2;
		SALSAUtility.Pcutf = config.Pcutf > 0 ? config.Pcutf : 0.85;
		SALSAUtility.Normalize = config.Normalize;
		SALSAUtility.ClusterFile = config.ClusterFile;
		SALSAUtility.SelectedClusters = new java.util.HashSet<Integer>();
		if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(config.SelectedClusters) && !"none".equals(config.SelectedClusters))
		{
            Pattern pat = Pattern.compile("[,]");
			String[] splits = pat.split(config.SelectedClusters.trim());
			for (String split : splits)
			{
				int cnum = Integer.parseInt(split);
				if (!SALSAUtility.SelectedClusters.contains(cnum))
				{
					SALSAUtility.SelectedClusters.add(cnum);
				}
			}
		}

		if (SALSAUtility.SelectedClusters.size() > 0 && !(new java.io.File(SALSAUtility.ClusterFile)).isFile())
		{
			throw new RuntimeException("Cluster file not specified to decide selected clusters");
		}

		SALSAUtility.IsClustersSelected = SALSAUtility.SelectedClusters.size() > 0;

		SALSAUtility.ManxcatRunName = !tangible.DotNetToJavaStringHelper.isNullOrEmpty(config.ManxcatRunName) ? config.ManxcatRunName : "Unspecified Run";
		SALSAUtility.ManxcatRunDescription = !tangible.DotNetToJavaStringHelper.isNullOrEmpty(config.ManxcatRunDescription) ? config.ManxcatRunDescription : "Description not specified";
        SALSAUtility.ConsoleDebugOutput = config.ConsoleDebugOutput;
        SALSAUtility.DebugPrintOption = config.DebugPrintOption;
        SALSAUtility.dataTypeSize = config.getDataTypeSize();
        SALSAUtility.endianness = config.isBigEndian() ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN;
	}

    public static void ManxcatControl(Hotsun.CalcfgSignature Calcfg, Hotsun.IntializeSignature InitializeApplication,
                                      Hotsun.InitializeParametersSignature InitializeParameters,
                                      Hotsun.SolveMatrixSignature SolveMatrix,
                                      Hotsun.WriteSolutionSignature WriteSolution,
                                      Hotsun.FindQlimitsSignature FindQlimits,
                                      Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct,
                                      Hotsun.SequelSignature Sequel) {
		//  Set up Specific Manxcat Applicatin
		InitializeApplication.invoke(); // Reads in Data

		while (Hotsun.InitializationLoopCount < Hotsun.InitializationLoops)
		{
			// Iterate over Initializations
			// First call to user routine
			Hotsun.numit = 0;
			InitializeParameters.invoke(Hotsun.CurrentSolution, Hotsun.InitializationLoopCount);
			ZeroSolution(Hotsun.CurrentSolution);
			MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
			GlobalParameterSet = true; // Set Indicator  that Global Parameters are set

			SALSAUtility.StartSubTimer(2);
			violat = Calcfg.invoke(Hotsun.CurrentSolution);
			SALSAUtility.StopSubTimer(2);

			Hotsun.tcalcfg = SALSAUtility.SubDurations[2];
			Hotsun.tsolve = 0.0;
			Hotsun.teigen = 0.0;
			Hotsun.NumberofCGIterations = -2;
			Hotsun.EigenvalueIndicator1 = 0;
			Hotsun.EigenvalueIndicator2 = 0;

			ManxcatCentral.launchQlimits(FindQlimits); // Find limits on Q

			Hotsun.ichsav = -1; // changed from = 1 in old Fortran which seems wrong
			Hotsun.chsave[0] = Hotsun.zerocr; // Unnecessary
			Hotsun.expchg = 0.0;
			Hotsun.delchi = 0.0;
			Hotsun.isaved = 0;
			Hotsun.bnderr = 0;
			Hotsun.igood = 0;
			Hotsun.CountToSD = 0;
			Hotsun.CountQSmall = 0;
			Hotsun.CountQLarge = 0;
			Hotsun.CountToRandom = 0;
			Hotsun.CountExercusion = 0;
			Hotsun.QgoodFactor = 1.0;

			Hotsun.Q = Hotsun.QHighInitialFactor * Hotsun.Qhigh;
			boolean QlessthanQlow = Hotsun.Q < Hotsun.Qlow;
			tangible.RefObject<Boolean> tempRef_QlessthanQlow = new tangible.RefObject<Boolean>(QlessthanQlow);
			SALSAUtility.SynchronizeMPIvariable(tempRef_QlessthanQlow);
			QlessthanQlow = tempRef_QlessthanQlow.argValue;

			if (QlessthanQlow)
			{
				Hotsun.Q = Hotsun.Qlow;
			}
			Hotsun.Qbest = Hotsun.Q;

			Hotsun.Qgood = -1.0;
			Hotsun.isdtry = 0;

			if (Hotsun.InitialSteepestDescents > 0)
			{
				Hotsun.isdtry = 1;
			}

			// Print Initial Trial
			SALSAUtility.SALSAStatus(ManxcatCentral.ResultDirectoryName, " Initial Chisq " + String.format("%.3f", ChisqPrintConstant * Hotsun.zerocr));

			if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
			{
				SALSAUtility.SALSAPrint(1,
                                        "\n-----------------------------------------------\n" + SALSAUtility.startTime
                                                                                                            .ToLocalTime() + " Iterations " + Hotsun.numit + " Chisq " + String


                                                .format("%.3f",
                                                        ChisqPrintConstant * Hotsun.zerocr) + " Parameters " + Hotsun
                                                .npar + " Data Points " + (new Long(
                                                Hotsun.ndata))
                                                .toString());
			}

			Hotsun.DecisionMethod = 0; //   Initial Conditions on Decisions
			Hotsun.DecisionMethod_1 = 0;
			Hotsun.DecisionMethod_2 = 0;
			Hotsun.DecisionMethod_3 = 0;
			Hotsun.DecisionMethod_4 = 0;
			Hotsun.DecisionMethod_5 = 0;
			Hotsun.DecisionMethod_LineSearch = 0;

			//  Save initial guess as best solution
			SaveBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
			BestSolutionSet = true;
			BestxshiftSet = false;

			int CoordinateWriteCount = 0;
			while (true)
			{
				// Iterate Chisq Solution
				//  Save diagonal elements removed as we no longer overwrite by adding Matquardt's good idea Q

				// Hotsun.PreviousSolution.param holds parameters of Starting Solution
				// while Hotsun.PreviousSolution.xshift was shift that was added to make this param
				CopySolution(Hotsun.PreviousSolution, Hotsun.CurrentSolution);

				//  Start code used each iteration whether or not last iteration succeeded
				Hotsun.materr = 0;

				boolean wefailed = false;
				Hotsun.DecisionMethod_LineSearch = 0;

				while (true)
				{
					// Solve Matrix and Find Good Q value if fails

					//  Process Options for Q getting Stuck
					if (Hotsun.CountQSmall < 0)
					{
						Hotsun.CountQSmall = 0;
					}

					if (Hotsun.CountQLarge < 0)
					{
						Hotsun.CountQLarge = 0;
					}
					double Qtest = 0.05 * Hotsun.Qhigh;
					Qtest = Math.max(2.0 * Hotsun.Qlow, Qtest);
					boolean Qisaboveaverage = Hotsun.Q > Qtest;
					tangible.RefObject<Boolean> tempRef_Qisaboveaverage = new tangible.RefObject<Boolean>(Qisaboveaverage);
					SALSAUtility.SynchronizeMPIvariable(tempRef_Qisaboveaverage);
					Qisaboveaverage = tempRef_Qisaboveaverage.argValue;

					if (Qisaboveaverage)
					{
						Hotsun.CountQSmall = 0;
						++Hotsun.CountQLarge;
					}
					else
					{
						Hotsun.CountQLarge = 0;
						++Hotsun.CountQSmall;
					}

					if (Hotsun.CountQSmall > Hotsun.ResetQLimit)
					{
						Hotsun.CountQSmall = -1;
						Hotsun.CountQLarge = 0;
						Hotsun.Q = Qtest * 4.0;
						Hotsun.DecisionMethod_5 = 1;
					}

					if (Hotsun.CountQLarge > Hotsun.ResetQLimit)
					{
						Hotsun.CountQLarge = -1;
						Hotsun.CountQSmall = 0;
						Hotsun.Q = Math.max(Hotsun.Qlow, 0.1 * Qtest);
						Hotsun.DecisionMethod_5 = 2;
					}

					// Process Steepest Descent Options
					if ((Hotsun.isdtry == 2) || (Hotsun.isdtry == 4))
					{
						Hotsun.isdtry = 0;
					}

					if ((Hotsun.isdtry == 1) || (Hotsun.isdtry == 3))
					{
						// Steepest Descent Solution  -- Q is NOT added
						Hotsun.AddMarquardtQDynamically = false;
						SteepestDescentSolution(Hotsun.CurrentSolution, new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
						++Hotsun.isdtry;
						Hotsun.Qsd = Hotsun.Q;
						Hotsun.materr = 0;
						break;
					}
					else
					{
						// Typical Full solution

						Hotsun.UseDiagonalScaling = Hotsun.UseDiagonalScalinginSolvers;

						Hotsun.FullSecondDerivative = false;

						if (config.FullSecondDerivativeOption == 1)
						{
							Hotsun.FullSecondDerivative = true;
						}

						if (Hotsun.UseDiagonalScaling)
						{
							SetupDiagonalScaling(Hotsun.CurrentSolution); // This does NOT change matrices
						}

						// Add in Marquardt Parameter (removed afterwards)
						Hotsun.AddMarquardtQDynamically = !Hotsun.AddMarquardtQExplicitly;
						AddSubtractMarquardtFactor(Hotsun.CurrentSolution, 1.0);

						//  Invert Modified Chisq Matrix
						SALSAUtility.StartSubTimer(0);
						boolean SolveMatrixSuccess = SolveMatrix.invoke(Hotsun.CurrentSolution.xshift, Hotsun.CurrentSolution);
						Hotsun.tsolve = SALSAUtility.StopSubTimer(0);
						AddSubtractMarquardtFactor(Hotsun.CurrentSolution, -1.0);

						if (SolveMatrixSuccess)
						{
							Hotsun.materr = 0;
							break;
						}
						else
						{
							// error in Matrix Solver -- either give up or throw away trial and increase Q
							if (Hotsun.materr > 1)
							{
								// Too many failures
								EndupManxcat(6, new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
								wefailed = true;
								break;
							}

							++Hotsun.materr;

							// First try in matrix failure is increasing Q
							Hotsun.Q = 4.0 * Hotsun.Q;

							boolean Qessentiallyzero = Hotsun.Q <= (2.0 * Hotsun.QHighInitialFactor);
							tangible.RefObject<Boolean> tempRef_Qessentiallyzero = new tangible.RefObject<Boolean>(Qessentiallyzero);
							SALSAUtility.SynchronizeMPIvariable(tempRef_Qessentiallyzero);
							Qessentiallyzero = tempRef_Qessentiallyzero.argValue;
							if (Qessentiallyzero)
							{
								Hotsun.Q = 1.0;
							}

							// Second try in matrix failure is steepest descent
							if (Hotsun.materr == 2)
							{
								Hotsun.isdtry = 3;
							}
							continue;
						}
					}
				} // End while Solving Matrix and finding good Q value if failure

				if (wefailed)
				{
					break;
				}

				// We only reach here if Matrix Solved Correctly by Steepest Descent or Marquardt
				// Old Manxcat used two step solution (LU and Solve)
				//  Here we just use direct Conjugate Gradient

				//  Calculate Predicted Chisq Change from this solution
				FindPredictedChisqChange(Hotsun.CurrentSolution, new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
				Hotsun.expchg = Hotsun.pred1 + Hotsun.pred2;

				if (Hotsun.pred2 > 0)
				{
					if (SALSAUtility.MPI_Rank == 0)
					{
						SALSAUtility.SALSAPrint(2, "\nIllegal pred2 " + String.format("%.3f", Hotsun.pred2 * ChisqPrintConstant) + " CurrentChisq " + String.format("%.3f", Hotsun.CurrentSolution.Chisquared * ChisqPrintConstant));
					}
				}

				// For derivative test avoid first iteration that could be special
				if (Hotsun.derivtest && (Hotsun.numit == 3))
				{
					DoDerivTest(new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
					Hotsun.derivtest = false;
				}
				Hotsun.doingderivtest = false;

				//  Set up parameters for next call to Calcfg
				//  This adds in estimated change to param and zeros first derivative, Second Derivative and Chisq (zerocr) in a  CurrentSolution
				// xshift was actually calculated from PREVIOUS Current Solution (now stored in Previous Solution)
				SALSABLAS.LinearCombineVector(Hotsun.CurrentSolution.param, 1.0, Hotsun.CurrentSolution.param, -1.0, Hotsun.CurrentSolution.xshift);
				MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
				GlobalParameterSet = true;
				ZeroSolution(Hotsun.CurrentSolution); // Does not vary param or xshift

				// Save Current best Chisq in circular list
				++Hotsun.ichsav;

				if (Hotsun.ichsav >= Hotsun.nbadgo)
				{
					Hotsun.ichsav = 0;
				}
				Hotsun.chsave[Hotsun.ichsav] = Hotsun.zeromn;

				Hotsun.numit++; // Increment Iteration count
				Hotsun.idata = 0;


				//  Call Calcfg to calculate Taylor expansion
				SALSAUtility.StartSubTimer(2);
				violat = Calcfg.invoke(Hotsun.CurrentSolution);
				SALSAUtility.StopSubTimer(2);

				//  Cope with boundary errors
				// In this case parameter values will have been changed by Calcfg
				if (violat)
				{
					Hotsun.bnderr++;
					Findxnorm(2); // Reset xshift and xnorm to reflect Calcfg resetting
				}
				else
				{
					Hotsun.bnderr = 0;
				}

				//   set change in chisq
				Hotsun.delchi = Hotsun.PreviousSolution.Chisquared - Hotsun.zerocr;

				//  Look at case where actual change significantly larger than expected
				// Only do search if last solution was also a success
				boolean LineSearchUsed = false;
				boolean randomized = false;
				boolean searchsuccess = false;
				double LineFactor = 1.0;
				double ExtraDecrease = 0.0;
				int LineIterations = 0;
				boolean Delchipositive = Hotsun.delchi >= 0.0;
				tangible.RefObject<Boolean> tempRef_Delchipositive = new tangible.RefObject<Boolean>(Delchipositive);
				SALSAUtility.SynchronizeMPIvariable(tempRef_Delchipositive);
				Delchipositive = tempRef_Delchipositive.argValue;

				if (Delchipositive)
				{
					if ((Hotsun.OmegaOption > 0) && (Hotsun.igood > 0))
					{
						boolean LookatOmegaOption = Hotsun.delchi > (Hotsun.Omega * Hotsun.expchg);
						tangible.RefObject<Boolean> tempRef_LookatOmegaOption = new tangible.RefObject<Boolean>(LookatOmegaOption);
						SALSAUtility.SynchronizeMPIvariable(tempRef_LookatOmegaOption);
						LookatOmegaOption = tempRef_LookatOmegaOption.argValue;

						if (LookatOmegaOption)
						{
							// Do a line search exploring larger shifts
							LineSearchUsed = true;
							LineFactorGuess = -1.0; // Explore larger solutions
							tangible.RefObject<Integer> tempRef_DecisionMethod_LineSearch = new tangible.RefObject<Integer>(Hotsun.DecisionMethod_LineSearch);
							tangible.RefObject<Double> tempRef_LineFactor = new tangible.RefObject<Double>(LineFactor);
							tangible.RefObject<Double> tempRef_ExtraDecrease = new tangible.RefObject<Double>(ExtraDecrease);
							tangible.RefObject<Integer> tempRef_LineIterations = new tangible.RefObject<Integer>(LineIterations);
							searchsuccess = LineSearch(Hotsun.CurrentSolution, Hotsun.PreviousSolution, tempRef_DecisionMethod_LineSearch, tempRef_LineFactor, tempRef_ExtraDecrease, tempRef_LineIterations, new Hotsun.Hotsun.CalcfgSignature(Calcfg), -1.0);
						Hotsun.DecisionMethod_LineSearch = tempRef_DecisionMethod_LineSearch.argValue;
						LineFactor = tempRef_LineFactor.argValue;
						ExtraDecrease = tempRef_ExtraDecrease.argValue;
						LineIterations = tempRef_LineIterations.argValue;
							Hotsun.DecisionMethod_LineSearch += 100;
							Hotsun.TotalSearchIterations += LineIterations;
						}
					}
				}
				else
				{
					if ((Hotsun.pred1 > 0.0) && (Hotsun.pred2 < 0.0))
					{
						double dchi1 = Hotsun.pred1 * ChisqPrintConstant;

						if (dchi1 < 0.5 * Hotsun.zerocr)
						{
							if (Hotsun.CountExercusion == 0)
							{
								++Hotsun.CountToRandom;
							}

							if ((Hotsun.CountToRandom >= Hotsun.RandomLimit) && ((Hotsun.maxit - Hotsun.numit) > 20))
							{
								Hotsun.CountToRandom = 0;
								randomized = true;
							}
							else
							{
								double dchi4 = Hotsun.delchi * ChisqPrintConstant;
								double alpha = dchi1 / (2.0 * (dchi1 - dchi4));
								LineFactorGuess = alpha;
								tangible.RefObject<Integer> tempRef_DecisionMethod_LineSearch2 = new tangible.RefObject<Integer>(Hotsun.DecisionMethod_LineSearch);
								tangible.RefObject<Double> tempRef_LineFactor2 = new tangible.RefObject<Double>(LineFactor);
								tangible.RefObject<Double> tempRef_ExtraDecrease2 = new tangible.RefObject<Double>(ExtraDecrease);
								tangible.RefObject<Integer> tempRef_LineIterations2 = new tangible.RefObject<Integer>(LineIterations);
								searchsuccess = LineSearch(Hotsun.CurrentSolution, Hotsun.PreviousSolution, tempRef_DecisionMethod_LineSearch2, tempRef_LineFactor2, tempRef_ExtraDecrease2, tempRef_LineIterations2, new Hotsun.Hotsun.CalcfgSignature(Calcfg), alpha);
							Hotsun.DecisionMethod_LineSearch = tempRef_DecisionMethod_LineSearch2.argValue;
							LineFactor = tempRef_LineFactor2.argValue;
							ExtraDecrease = tempRef_ExtraDecrease2.argValue;
							LineIterations = tempRef_LineIterations2.argValue;
								Hotsun.TotalSearchIterations += LineIterations;
								Hotsun.zerocr = Hotsun.CurrentSolution.Chisquared;
								Hotsun.delchi = Hotsun.PreviousSolution.Chisquared - Hotsun.zerocr;
								Delchipositive = Hotsun.delchi >= 0.0;
								LineSearchUsed = true;
							}
						}
					}
				}

				// Change Information
				double size1 = SALSABLAS.VectorScalarProduct(Hotsun.CurrentSolution.xshift, Hotsun.CurrentSolution.xshift);
				double size2 = SALSABLAS.VectorScalarProduct(Hotsun.CurrentSolution.first, Hotsun.CurrentSolution.first);
				double Changexshiftfirst = SALSABLAS.VectorScalarProduct(Hotsun.CurrentSolution.first, Hotsun.CurrentSolution.xshift);
				Changexshiftfirst = Changexshiftfirst / Math.sqrt(size1 * size2);
				double Changeinfirst = 0.0;
				double Changeinxshift = 0.0;
				double xshiftRatio = 0.0;
				double firstRatio = 0.0;

				// Changeinfirst is scalar product of first derivative at I and First Derivative at II / Size ( I * II)
				// firstratio is Ratio( Size first for II / Size first for I)
				// In most algorithms I is Best solution and II current solution
				// If exploring away from best I is previous solution
				double size3 = SALSABLAS.VectorScalarProduct(Hotsun.PreviousSolution.first, Hotsun.PreviousSolution.first);
				Changeinfirst = SALSABLAS.VectorScalarProduct(Hotsun.CurrentSolution.first, Hotsun.PreviousSolution.first);
				Changeinfirst = Changeinfirst / Math.sqrt(size2 * size3);
				firstRatio = Math.sqrt(size2 / size3);

				// Changeinxshift is scalar product of xshift at I and xshift at II / Size ( I * II)
				// xshiftratio is Ratio( Size xshift for II / Size xshift for I)
				// In most algorithms I is Best solution and II current solution
				// If exploring away from best I is previous solution
				boolean shiftchangeset = false;
				if (Hotsun.numit > 1)
				{
					shiftchangeset = true;
					double size4 = SALSABLAS.VectorScalarProduct(Hotsun.PreviousSolution.xshift, Hotsun.PreviousSolution.xshift);
					Changeinxshift = SALSABLAS.VectorScalarProduct(Hotsun.PreviousSolution.xshift, Hotsun.CurrentSolution.xshift);
					Changeinxshift = Changeinxshift / Math.sqrt(size1 * size4);
					xshiftRatio = Math.sqrt(size1 / size4);
				}

				//   Reset change in chisq
				Hotsun.delchi = Hotsun.PreviousSolution.Chisquared - Hotsun.zerocr;

				//  Timing Information
				// Hotsun.tsolve and Hotsun.teigen set in StopTimer
				Hotsun.tcalcfg = SALSAUtility.SubDurations[2];
				SALSAUtility.InterimTiming();
				Hotsun.TotalTimeUsed = SALSAUtility.HPDuration;

				//  Print Summary of this solution
				SALSAUtility.SALSAStatus(ManxcatCentral.ResultDirectoryName,
                                         "In Loop " + Hotsun.InitializationLoopCount + " Iteration " + Hotsun.numit + " Chisq " + String

                                                 .format("%.3f",
                                                         ChisqPrintConstant * Hotsun.zerocr) + " Best Chisq " + String

                                                 .format("%.3f", ChisqPrintConstant * Hotsun.BestSolution.Chisquared));

				if ((SALSAUtility.DebugPrintOption > 1) && (SALSAUtility.MPI_Rank == 0))
				{
					String shift1 = "Unset";
					String shift2 = "unset";
					if (shiftchangeset)
					{
						shift1 = String.format("%.4f", Changeinxshift);
						shift2 = String.format("%.4f", xshiftRatio);
					}
					SALSAUtility.SALSAPrint(2,
                                            "\nLoop " + Hotsun.InitializationLoopCount + " Iteration " + Hotsun.numit
                                                    + " Q " + String
                                                    .format("%.4E", Hotsun.Q) + " Qlow " + String
                                                    .format("%.3f",
                                                            Hotsun.Qlow) + " (" + Hotsun.EigenvalueIndicator2 + ") " +
                                                    "Qhigh " + String
                                                    .format("%.3f",
                                                            Hotsun.Qhigh) + " (" + Hotsun.EigenvalueIndicator1 + ") " +
                                                    "Qgood " + String
                                                    .format("%.4E", Hotsun.Qgood) + " Trace Q " + String
                                                    .format("%.3f", Hotsun.ChisqMatrixTrace) + " Norm Q " + String
                                                    .format("%.3f",
                                                            Hotsun.ChisqMatrixNorm) + "\nDecision " + Hotsun
                                                    .DecisionMethod + " 1:" + Hotsun.DecisionMethod_1 + " 2:" +
                                                    Hotsun.DecisionMethod_2 + " 3:" + Hotsun.DecisionMethod_3 + " 4:"
                                                    + Hotsun.DecisionMethod_4 + " 5:" + Hotsun.DecisionMethod_5 + " " +
                                                    "Counts Q Low " + Hotsun.CountQSmall + " Q High "
                                                    + Hotsun.CountQLarge + " SD Count " + Hotsun.CountToSD + " " +
                                                    "Randoms " +
                                                    Hotsun.CountToRandom + " Exercusions " + Hotsun.CountExercusion +
                                                    " SD Opt " + Hotsun.isdtry + " Solver Errs " + Hotsun.materr + " " +
                                                    "CG Iters " + Hotsun.NumberofCGIterations + "\nScalProd Shift" +
                                                    ".First " + String
                                                    .format("%.4f",
                                                            Changexshiftfirst) + " Change Shift from Prev DotProd " +
                                                    shift1 + " and Change Shift in norm " + shift2 + "\nChange in " +
                                                    "first deriv from Prev DotProd " + String


                                                    .format("%.4f",
                                                            Changeinfirst) + " and in first deriv norm " + String
                                                    .format("%.4f", firstRatio));

					String chisqmethod = "Marquardt's Algorithm";

					if (Hotsun.FullSecondDerivative)
					{
						chisqmethod = "Corrected Marquardt";
					}

					if ((Hotsun.isdtry == 2) || (Hotsun.isdtry == 4))
					{
						chisqmethod = "Steepest Descent";
						if (Hotsun.isdtry == 4)
						{
							chisqmethod += " from Singularity";
						}
						if (Hotsun.FullSecondDerivative)
						{
							chisqmethod = "Corrected Steep Desc";
						}
					}

					String chisqstatus = "OK Parameters";

					if (violat)
					{
						chisqstatus = "Out of Range Parameter";
					}
					double dchi1 = Hotsun.pred1 * ChisqPrintConstant;
					double dchi2 = Hotsun.pred2 * ChisqPrintConstant;
					double dchi3 = Hotsun.pred3 * ChisqPrintConstant;
					double dchi4 = Hotsun.delchi * ChisqPrintConstant;

					if (LineSearchUsed)
					{
						SALSAUtility.SALSAPrint(2, "Linesearch Chisq gains " + String
                                .format("%.3f", ChisqPrintConstant * ExtraDecrease) + " Shift Factor " + String
                                .format("%.3f", LineFactor) + " with guess " + String
                                .format("%.3f",
                                        LineFactorGuess) + " Iterations " + LineIterations + " Method " + Hotsun.DecisionMethod_LineSearch);
					}

					// Output parameter values if not too many
					String parametervalues = "";

					if ((!Hotsun.DecomposeParameters) && (Hotsun.npar < 20))
					{
						parametervalues = "\nParameters:";

						for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
						{
							for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
							{
								double tmp = Hotsun.CurrentSolution.param[LongIndex][LocalVectorIndex];
								parametervalues += " " + String.format("%.4E", tmp);
							}
						}
					}

					String labelCG = "Not Scaled ";
					if (Hotsun.NumberofCGIterations > 0)
					{
						labelCG = "Scaled ";
						Hotsun.tsolve = Hotsun.tsolve / Hotsun.NumberofCGIterations;
					}
					String labeleigen = "Not Scaled ";
					if ((Hotsun.EigenvalueIndicator1 > 0) && (Hotsun.EigenvalueIndicator2 > 0))
					{
						labeleigen = "Scaled ";
						Hotsun.teigen = Hotsun.teigen / (Hotsun.EigenvalueIndicator1 + Hotsun.EigenvalueIndicator2);
					}

					SALSAUtility.SALSAPrint(2, "Chisq " + String.format("%.3f", ChisqPrintConstant * Hotsun.zerocr) + " Best Chisq " + String.format("%.3f", ChisqPrintConstant * Hotsun.BestSolution.Chisquared) + " " + chisqmethod + " " + chisqstatus + " Calcfg Time " + String.format("%.0f", Hotsun.tcalcfg) + " Matrix Solve " + labelCG + String.format("%.1f", Hotsun.tsolve) + " Eigen " + labeleigen + String.format("%.1f", Hotsun.teigen) + " Total " + String.format("%.0f", Hotsun.TotalTimeUsed) + "\n1st Deriv Reduction " + String.format("%.4f", dchi1) + " Exp 2nd " + String.format("%.4f", dchi2) + " Other 2nd " + String.format("%.4f", dchi3) + " Actual Reduction " + String.format("%.4f", dchi4) + parametervalues);

					Hotsun.tsolve = 0.0;
					Hotsun.teigen = 0.0;
					Hotsun.tcalcfg = 0.0;
					Hotsun.NumberofCGIterations = -2;
					Hotsun.EigenvalueIndicator1 = 0;
					Hotsun.EigenvalueIndicator2 = 0;

				} // End Print Out

				// Save Results periodically
				CoordinateWriteCount++;
				if (CoordinateWriteCount == Hotsun.CoordinateWriteFrequency)
				{
					CoordinateWriteCount = 0;
					CopySolution(Hotsun.SearchSolution1, Hotsun.CurrentSolution);
					CopySolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
					double chisave = Hotsun.zerocr;
					Hotsun.zerocr = Hotsun.zeromn;

					SALSABLAS.zrword(Hotsun.perr);
					MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
					GlobalParameterSet = true; // Set Indicator  that Global Parameters are set

					WriteCoordinates(0, new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution));

					CopySolution(Hotsun.CurrentSolution, Hotsun.SearchSolution1);
					Hotsun.zerocr = chisave;
					MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
				}

				//  Examine Nature of Solution
				Hotsun.DecisionMethod = 1; // Default non Initialization
				Hotsun.DecisionMethod_1 = 0;
				Hotsun.DecisionMethod_2 = 0;
				Hotsun.DecisionMethod_3 = 0;
				Hotsun.DecisionMethod_4 = 0;
				Hotsun.DecisionMethod_5 = 0;

				if (Hotsun.CountExercusion > 0)
				{
					++Hotsun.CountExercusion;
					Hotsun.DecisionMethod_1 = 1;

					if (Hotsun.CountExercusion >= Hotsun.ExercusionLimit)
					{
						Hotsun.DecisionMethod_1 = 2;
						boolean CurrentSolutionBest = Hotsun.CurrentSolution.Chisquared < Hotsun.BestSolution.Chisquared;
						tangible.RefObject<Boolean> tempRef_CurrentSolutionBest = new tangible.RefObject<Boolean>(CurrentSolutionBest);
						SALSAUtility.SynchronizeMPIvariable(tempRef_CurrentSolutionBest);
						CurrentSolutionBest = tempRef_CurrentSolutionBest.argValue;

						if (!CurrentSolutionBest)
						{
							RestoreBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
							Hotsun.isaved = 2;
						}
						Hotsun.CountExercusion = 0;
						Hotsun.CountToRandom = 0;
					}
				}

				if (Delchipositive)
				{
					// Chisq Decreased - Success

					// If necessary, save current best parameters
					// We save true (diagonal of) second derivative matrix -- not one overwritten by Marquardt
					boolean CurrentSolutionBest = Hotsun.CurrentSolution.Chisquared < Hotsun.BestSolution.Chisquared;
					tangible.RefObject<Boolean> tempRef_CurrentSolutionBest2 = new tangible.RefObject<Boolean>(CurrentSolutionBest);
					SALSAUtility.SynchronizeMPIvariable(tempRef_CurrentSolutionBest2);
					CurrentSolutionBest = tempRef_CurrentSolutionBest2.argValue;

					if (CurrentSolutionBest)
					{
						SaveBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
						BestSolutionSet = true;
						BestxshiftSet = true;
						Hotsun.CountExercusion = 0;
					}

					Hotsun.DecisionMethod = 2;
					Hotsun.succ = true;
					Hotsun.isaved = 0;

					if (Hotsun.isdtry == 2)
					{
						Hotsun.DecisionMethod = 3;
					}

					if (Hotsun.isdtry == 4)
					{
						Hotsun.DecisionMethod = 4;
						boolean QsdTest = Hotsun.Q < (Hotsun.Qsd - 0.0000000001);
						tangible.RefObject<Boolean> tempRef_QsdTest = new tangible.RefObject<Boolean>(QsdTest);
						SALSAUtility.SynchronizeMPIvariable(tempRef_QsdTest);
						QsdTest = tempRef_QsdTest.argValue;

						if (QsdTest)
						{
							Hotsun.DecisionMethod = 5;
						}
					}

					if (LineSearchUsed)
					{
						Hotsun.DecisionMethod += 20;

						if (!searchsuccess)
						{
							Hotsun.DecisionMethod += 20;
						}
					}
				}
				else
				{
					// Chisq Increases -- Failure
					Hotsun.DecisionMethod = 11;
					Hotsun.isaved = 1;

					if (randomized)
					{
						Hotsun.DecisionMethod = 12;
						Hotsun.CountExercusion = 1;
					}
					Hotsun.succ = false;

					if (Hotsun.CountExercusion == 0)
					{
						RestoreBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
						Hotsun.isaved = 2;
					}
				}

				//  Decide if to give up ghost and call this the last Iteration
				int readInstruction = SALSAUtility.SALSAGracefulend(ManxcatCentral.ResultDirectoryName);
				readInstruction = SALSAUtility.synchronizeMPIVariable(readInstruction);
				if (readInstruction == 0)
				{
					Hotsun.InitializationLoops = Hotsun.InitializationLoopCount + 1;
					EndupManxcat(7, WriteSolution, Calcfg, GlobalMatrixVectorProduct);
					break;
				}
				else if (readInstruction > 0)
				{
					Hotsun.maxit = readInstruction;
				}

				boolean TimeExceeded = (Hotsun.timmax > 0.0) && (Hotsun.TotalTimeUsed > Hotsun.timmax);
				tangible.RefObject<Boolean> tempRef_TimeExceeded = new tangible.RefObject<Boolean>(TimeExceeded);
				SALSAUtility.SynchronizeMPIvariable(tempRef_TimeExceeded);
				TimeExceeded = tempRef_TimeExceeded.argValue;

				if (TimeExceeded)
				{
					// Time limit reached
					EndupManxcat(3, new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
					break;
				}

				if (Hotsun.numit >= Hotsun.maxit)
				{
					// Iteration limit reached
					EndupManxcat(1, new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
					break;
				}

				boolean SmallExpectedChisqChange = (Hotsun.expchg * ChisqPrintConstant) <= Hotsun.dellim;
				tangible.RefObject<Boolean> tempRef_SmallExpectedChisqChange = new tangible.RefObject<Boolean>(SmallExpectedChisqChange);
				SALSAUtility.SynchronizeMPIvariable(tempRef_SmallExpectedChisqChange);
				SmallExpectedChisqChange = tempRef_SmallExpectedChisqChange.argValue;

				if (SmallExpectedChisqChange)
				{
					// Expected change in chisq <= preset limit
					EndupManxcat(2, new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
					break;
				}

				if (Hotsun.bnderr > Hotsun.bnderrLimit)
				{
					// Boundary Value limit reached
					EndupManxcat(5, new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
					break;
				}

				// Stop if too little progress made in nbadggo steps
				// Do this by comparing current best chisq with that found nbadgo iterations before
				if (Hotsun.numit > Hotsun.nbadgo)
				{
					int itest = Hotsun.ichsav + 1;

					if (itest >= Hotsun.nbadgo)
					{
						itest = 0;
					}

					boolean TooLittleProgress = (Hotsun.chsave[itest] - Hotsun.zeromn) < Hotsun.dellim;
					tangible.RefObject<Boolean> tempRef_TooLittleProgress = new tangible.RefObject<Boolean>(TooLittleProgress);
					SALSAUtility.SynchronizeMPIvariable(tempRef_TooLittleProgress);
					TooLittleProgress = tempRef_TooLittleProgress.argValue;

					if (TooLittleProgress)
					{
						EndupManxcat(4, new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
						break;
					}
				}

				//  Another Iteration called for -- Decide on Strategy

				//  Reasonable Success -- Decrease Marquardt Parameter if going very well
				//  Namely actual change is larger that Fletcher's rho * Expected change estimated from Taylor expansion
				boolean GoodenoughChisqChange = Hotsun.delchi >= (Hotsun.rho * Hotsun.expchg);
				tangible.RefObject<Boolean> tempRef_GoodenoughChisqChange = new tangible.RefObject<Boolean>(GoodenoughChisqChange);
				SALSAUtility.SynchronizeMPIvariable(tempRef_GoodenoughChisqChange);
				GoodenoughChisqChange = tempRef_GoodenoughChisqChange.argValue;

				if (GoodenoughChisqChange)
				{
					launchQlimits(new Hotsun.Hotsun.FindQlimitsSignature(FindQlimits)); // Reset Q limits every now and then
					Hotsun.DecisionMethod_4 = 2;

					Hotsun.Qgood = Hotsun.Q;
					++Hotsun.igood;
					boolean Qessentiallyzero = Hotsun.Q < 0.0000000001;
					tangible.RefObject<Boolean> tempRef_Qessentiallyzero2 = new tangible.RefObject<Boolean>(Qessentiallyzero);
					SALSAUtility.SynchronizeMPIvariable(tempRef_Qessentiallyzero2);
					Qessentiallyzero = tempRef_Qessentiallyzero2.argValue;

					if (Qessentiallyzero)
					{
						Hotsun.DecisionMethod_2 = 1;
						continue; // Next Iteration
					}

					//  Modest change are those "good" iterations that change Chisq by less than Feltcher's Sigma *   Expected change estimated from Taylor expansion
					//  Good but modest chisq changes leave Q unchanged
					boolean ChisqChangemodest = Hotsun.delchi < (Hotsun.sigma * Hotsun.expchg);
					tangible.RefObject<Boolean> tempRef_ChisqChangemodest = new tangible.RefObject<Boolean>(ChisqChangemodest);
					SALSAUtility.SynchronizeMPIvariable(tempRef_ChisqChangemodest);
					ChisqChangemodest = tempRef_ChisqChangemodest.argValue;

					if (ChisqChangemodest)
					{
						Hotsun.DecisionMethod_2 = 2;
						++Hotsun.CountToSD;

						if (Hotsun.CountToSD >= Hotsun.SDLimit)
						{
							Hotsun.CountToSD = 0;

							if (Hotsun.isdtry == 0)
							{
								Hotsun.isdtry = 1;
							}

							Hotsun.DecisionMethod_2 = 3;
						}
						continue; // Next Iteration
					}

					//  This is case of really good Chisq change -- decrease Q
					Hotsun.CountToSD = 0;
					boolean Qalreadysmall = Hotsun.Q <= (Hotsun.Qscale * Hotsun.Qlow);
					tangible.RefObject<Boolean> tempRef_Qalreadysmall = new tangible.RefObject<Boolean>(Qalreadysmall);
					SALSAUtility.SynchronizeMPIvariable(tempRef_Qalreadysmall);
					Qalreadysmall = tempRef_Qalreadysmall.argValue;

					if (Qalreadysmall)
					{
						if (Hotsun.igood >= 3)
						{
							// Reduce Q below Qlow if really good convergence
							Hotsun.QgoodFactor *= Hotsun.QgoodReductionFactor;
							Hotsun.Q = Hotsun.Qlow * Hotsun.QgoodFactor;
							Hotsun.DecisionMethod_2 = 4;
						}
						else
						{
							Hotsun.DecisionMethod_2 = 5;
						}
						continue; // Next Iteration
					}
					else
					{
						Hotsun.DecisionMethod_2 = 6;
						Hotsun.Q = Math.max(Math.sqrt(Hotsun.Q * Hotsun.Qlow), Hotsun.Q * 0.33);
						continue; // Next Iteration
					}
				} // End case when reasonable success -- possibly decrease Q

					//  Unexpectedly small or negative change -- Increase Marquardt Parameter
					// Reset limits Qlow and Qhigh on allowed range of Q values if needed
				else
				{
					Hotsun.DecisionMethod_3 = 1;

					Hotsun.igood = 0; // Restart count of good solutions
					Hotsun.QgoodFactor = 1.0; // Restart factor for lowering Q below Qlow
					++Hotsun.CountToSD;

					if (Hotsun.CountToSD >= Hotsun.SDLimit)
					{
						Hotsun.CountToSD = 0;

						if (Hotsun.isdtry == 0)
						{
							Hotsun.isdtry = 1;
						}

						Hotsun.DecisionMethod_3 = 2;
						continue; // Next Iteration
					}
					boolean Qlookswrong = false;
					boolean qlimitsvalid = true;
					boolean QnearQhigh = false;

					if (Hotsun.CurrentSolution.IterationCalculated != Hotsun.IterationforQlimit)
					{
						qlimitsvalid = false;
					}
					boolean totaltest = qlimitsvalid || (Hotsun.numit < ManxcatCentral.IterationtorecalculateQLimits);

					if (!totaltest)
					{
						// Qlow and Qhigh need to be reset
						Hotsun.DecisionMethod_4 = 1;
						launchQlimits(new Hotsun.Hotsun.FindQlimitsSignature(FindQlimits));

					} // End case Qlow and Qhigh need to be reset

					QnearQhigh = Hotsun.Qhigh <= (Hotsun.Qscale * Hotsun.Q);
					tangible.RefObject<Boolean> tempRef_QnearQhigh = new tangible.RefObject<Boolean>(QnearQhigh);
					SALSAUtility.SynchronizeMPIvariable(tempRef_QnearQhigh);
					QnearQhigh = tempRef_QnearQhigh.argValue;

					if (QnearQhigh)
					{
						// Q is near enough Qhigh and it still didn't work very well! That's very disappointing

						Hotsun.DecisionMethod_3 = 3;

						if (Hotsun.isdtry != 0)
						{
							// Previous Round was not a full Chisq calculation
							Hotsun.Q *= 0.5;
							Qlookswrong = false;
							Hotsun.DecisionMethod_3 = 4;
							Hotsun.isdtry = 2;
						}
						else
						{
							// This will select Steepest Descent next time
							Qlookswrong = true;
						}
					} // End case Q is near enough Qhigh
					else
					{
						Qlookswrong = false;
						boolean Qgoodsmall = Hotsun.Qgood <= (Hotsun.Qscale * Hotsun.Q);
						tangible.RefObject<Boolean> tempRef_Qgoodsmall = new tangible.RefObject<Boolean>(Qgoodsmall);
						SALSAUtility.SynchronizeMPIvariable(tempRef_Qgoodsmall);
						Qgoodsmall = tempRef_Qgoodsmall.argValue;

						if (Qgoodsmall)
						{
							Hotsun.Q = Math.sqrt(Hotsun.Qlow * Hotsun.Qhigh);
							Hotsun.DecisionMethod_3 = 5;
						}
						else
						{
							Hotsun.Q = Math.sqrt(Hotsun.Q * Hotsun.Qgood);
							Hotsun.DecisionMethod_3 = 6;
						}
					}

					tangible.RefObject<Boolean> tempRef_Qlookswrong = new tangible.RefObject<Boolean>(Qlookswrong);
					SALSAUtility.SynchronizeMPIvariable(tempRef_Qlookswrong);
					Qlookswrong = tempRef_Qlookswrong.argValue;

					if (Qlookswrong)
					{
						Hotsun.DecisionMethod_3 += 10;

						if (Hotsun.isdtry == 0)
						{
							Hotsun.isdtry = 1;
						}

						Hotsun.Q = 0.5 * Hotsun.Qbest;
					}
				} // End Case when unexpectedly poor performance -- increase Q

				//  Another Iteration called for -- End decision on Strategy

			} // End Iterate Chisq Solution with while(true)

		} // End while over Initialization Iteration Loops

		// Perform follow-up
		Sequel.invoke();

		System.out.println("Completed");

		return;
	}

	// End Control

	// Explain Decision codes
	// DecisionMethod = 0  	Initial Conditions
	// DecisionMethod = 1	Default case if not initialization
	// DecisionMethod = 2	Default success -- chisq decreases
	// DecisionMethod = 3	chisq decreases: isdtry=2 (Steepest Descent used with no matrix failure) reset isdtry=0
	// DecisionMethod = 4	chisq decreases: isdtry=3 (Steepest Descent used with matrix failure) Q >=than previous steepest value
	// DecisionMethod = 5	chisq decreases: isdtry=3 (Steepest Descent used with matrix failure) Q < than previous steepest value; reset isdtry=04
	// DecisionMethod = 11	Default failure -- chisq increases
	// DecisionMethod = 12	Forced failure -- chisq increases and Exercusion started

	// DeltaDecisionMethod = 20	Success Chisq decreases and use Line Search Successfully
	// DeltaDecisionMethod = 40	Success Chisq decreases and use Line Search Unsuccessfully

	// DecisionMethod_1 = 0 No Exercusion
	// DecisionMethod_1 = 1 End Exercusion
	// DecisionMethod_1 = 2 Continue Exercusion

	// DecisionMethod_2 = 1	Good success as measured by rho and current Q essentially zero -- Final value
	// DecisionMethod_2 = 2	Good success as measured by rho and Modest success as measured by sigma -- Final Value
	// DecisionMethod_2 = 3	Good Success as measured by sigma and rho -- Q already small and Hotsun.igood >= 3
	// DecisionMethod_2 = 4	Good Success as measured by sigma and rho -- Q already small and Hotsun.igood < 3
	// DecisionMethod_2 = 5	Good Success as measured by sigma and rho -- Q not already small
	// DecisionMethod_2 = 6	Good success as measured by rho and Modest success as measured by sigma -- Try Steepest descent as happened too often

	// DecisionMethod_3 = 1 Case of Chisq increasing or very small change
	// DecisionMethod_3 = 2 Q looks too small -- Case of Chisq increasing or very small change
	// DecisionMethod_3 = 3 Q is near Qhigh Use Steepest Descent as not used last time -- Case of Chisq increasing or very small change
	// DecisionMethod_3 = 4 Q is near Qhigh Do not use Steepest Descent as used last time -- Case of Chisq increasing or very small change
	// DecisionMethod_3 = 5 Q NOT near enough Qhigh and current Qgood quite small -- Case of Chisq increasing or very small change
	// DecisionMethod_3 = 6 Q NOT near enough Qhigh and current Qgood quite large -- Case of Chisq increasing or very small change

	// DeltaDecisionMethod_3 = 10 Chose Steepest Descent - Case of Chisq increasing or very small change

	// DecisionMethod_4 = 1 Need to reset limits on Q -- Case of Chisq increasing or very small change
	// DecisionMethod-4 = 2 Need to reset limits on Q -- Every now and then

	// DecisionMethod_5 = 1 Q too small for too long so reset SET whether success or failure
	// DecisionMethod_5 = 2 Q too big for too long so reset SET whether success or failure

	// Iteration has ended -- Clean up in EndupManxcat
	//  ReasonToStop = 1 Iteration Limit
	//  ReasonToStop = 2 Expected Change in Chisq too small
	//  ReasonToStop = 3 Time exhausted
	//  ReasonToStop = 4 Progress too
	//  ReasonToStop = 5 Boundary value limit reached
	//  ReasonToStop = 6 Matrix Singular even though Q added

	public static void EndupManxcat(int ReasonToStop, Hotsun.WriteSolutionSignature WriteSolution, Hotsun.CalcfgSignature Calcfg, Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
	{

		if (Hotsun.isaved != 0)
		{
			// Copy residuals func from BestSolution to CurrentSolution -- This Logic NOT IMPLEMENTED in Current Manxcat
			Hotsun.zerocr = Hotsun.zeromn;
			if (Hotsun.isaved == 1)
			{
				RestoreBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
			}
		}
		// Hotsun.tsolve and Hotsun.teigen set in StopTimer
		Hotsun.tcalcfg = SALSAUtility.SubDurations[2];
		SALSAUtility.EndTiming();
		Hotsun.TotalTimeUsed = SALSAUtility.HPDuration;

		SALSAUtility.SALSAStatus(ManxcatCentral.ResultDirectoryName,
                                 "End Loop " + Hotsun.InitializationLoopCount + " Chisq " + String
                                         .format("%.3f", ChisqPrintConstant * Hotsun.zerocr));

		if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
		{

			SALSAUtility.SALSAPrint(1, "\n-------------------" + SALSAUtility.endTime
                                                                             .ToLocalTime() + "\nIterations " + Hotsun.numit + " Chisq " + String

                    .format("%.3f", ChisqPrintConstant * Hotsun.zerocr) + " Q " + String
                    .format("%.3f", Hotsun.Q) + " Qlow " + String.format("%.3f", Hotsun.Qlow) + " Qhigh " + String
                    .format("%.3f", Hotsun.Qhigh) + " Qgood " + String
                    .format("%.3f", Hotsun.Qgood) + " Trace Q " + String
                    .format("%.3f", Hotsun.ChisqMatrixTrace) + " Norm Q " + String
                    .format("%.3f", Hotsun.ChisqMatrixNorm));
			String bnderrlabel = "";

			if (Hotsun.bnderr != 0)
			{
				bnderrlabel = " Boundary Violations " + Hotsun.bnderr + " ";
			}
			String EndReason = "Correct Convergence";

			if (ReasonToStop == 1)
			{
				EndReason = "Iteration Limit";
			}

			if (ReasonToStop == 3)
			{
				EndReason = "Time Cut";
			}

			if (ReasonToStop == 4)
			{
				EndReason = "No Progress in " + Hotsun.nbadgo + " Iterations";
			}

			if (ReasonToStop == 5)
			{
				EndReason = "Boundary Violations in Calcfg";
			}

			if (ReasonToStop == 6)
			{
				EndReason = "Matrix Singular";
			}

			if (ReasonToStop == 7)
			{
				EndReason = "User End";
			}

			// Output parameter values if not too many
			String parametervalues = "";

			if ((!Hotsun.DecomposeParameters) && (Hotsun.npar < 20))
			{
				parametervalues = "\nParameters:";

				for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
				{
					for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
					{
						double tmp = Hotsun.CurrentSolution.param[LongIndex][LocalVectorIndex];
						parametervalues += " " + String.format("%.4E", tmp);
					}
				}
			}

			SALSAUtility.SALSAPrint(1,
                                    "Loop " + Hotsun.InitializationLoopCount + " " + bnderrlabel + EndReason + " CG Failures " + Hotsun.TotalCGFailures + " Iterations " + Hotsun.TotalCGIterations + " Eigenvalue Failures " + Hotsun.TotalPowerFailures + " Iterations " + Hotsun.TotalPowerIterations + " Search Iterations " + Hotsun.TotalSearchIterations + " Total Time " + String

                                            .format("%.1f", Hotsun.TotalTimeUsed) + " Calcfg Time " + String
                                            .format("%.1f", Hotsun.tcalcfg) + " Per Iter Solve " + String.format("%.1f",
                                                                                                                 SALSAUtility.SubDurations[0] / Hotsun.TotalCGIterations) + " Per Iter Eigen " + String
                                            .format("%.1f",
                                                    SALSAUtility.SubDurations[1] / Hotsun.TotalPowerIterations) + parametervalues);
			Hotsun.tsolve = 0.0;
			Hotsun.teigen = 0.0;
			Hotsun.tcalcfg = 0.0;
			Hotsun.NumberofCGIterations = -2;
			Hotsun.EigenvalueIndicator1 = 0;
			Hotsun.EigenvalueIndicator2 = 0;
		}

		// Output Chisqcomponent
		if (SALSAUtility.NumberFixedPoints > 0)
		{
			CalculateChisqComponents(new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct), Hotsun.CurrentSolution, ChisqPrintConstant);
		}

		//      Process Iteration over Initialization parameters
		if (Hotsun.InitializationLoopCount == 0)
		{
			CopySolution(Hotsun.BestLoopedSolution, Hotsun.CurrentSolution);
			Hotsun.BestChisqLoop = 0;
			if (Hotsun.InitializationLoopCount != (Hotsun.InitializationLoops - 1))
			{
				WriteOutParameters(new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
			}
		}
		else
		{
			if (Hotsun.CurrentSolution.Chisquared < Hotsun.BestLoopedSolution.Chisquared)
			{
				CopySolution(Hotsun.BestLoopedSolution, Hotsun.CurrentSolution);
				Hotsun.BestChisqLoop = Hotsun.InitializationLoopCount;

				if (Hotsun.InitializationLoopCount != (Hotsun.InitializationLoops - 1))
				{
					WriteOutParameters(new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));
				}
			}
		}
		Hotsun.InitLoopChisq[Hotsun.InitializationLoopCount] = Hotsun.CurrentSolution.Chisquared;
		Hotsun.InitializationLoopCount++;

		if (Hotsun.InitializationLoopCount < Hotsun.InitializationLoops)
		{
			return;
		}
		WriteOutParameters(new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution), new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct));

	} // End EndupMancat(int ReasonToStop)

	public static void WriteOutParameters(Hotsun.WriteSolutionSignature WriteSolution, Hotsun.CalcfgSignature Calcfg, Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
	{
		//  Set up really best Solution
		RestoreBestSolution(Hotsun.CurrentSolution, Hotsun.BestLoopedSolution);

		Hotsun.HotsunComment = ""; // Build final comment here

		//  Output Information on Initialization Loop
		if (Hotsun.InitializationLoopCount == Hotsun.InitializationLoops)
		{
			if (Hotsun.InitializationLoops > 1)
			{
				String firstcomment = "Best Initial Condition " + Hotsun.BestChisqLoop + " out of " + Hotsun
                        .InitializationLoops;
				if (SALSAUtility.MPI_Rank == 0)
				{
					SALSAUtility.SALSAPrint(1, firstcomment);
				}

				String ChisqList = "Chisq List";

				for (int localloop = 0; localloop < Hotsun.InitializationLoops; localloop++)
				{
					ChisqList += " " + String.format("%.3f", ChisqPrintConstant * Hotsun.InitLoopChisq[localloop]);
				}
				firstcomment += "\n" + ChisqList;
				if (SALSAUtility.MPI_Rank == 0)
				{
					SALSAUtility.SALSAPrint(1, ChisqList);
				}
				Hotsun.HotsunComment = firstcomment;
			}
			else
			{
				Hotsun.HotsunComment = String.format("%.3f", Hotsun.CurrentSolution.Chisquared * ChisqPrintConstant);
			}
			// Output Chisqcomponent
			if (SALSAUtility.NumberFixedPoints > 0)
			{
				CalculateChisqComponents(new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct), Hotsun.CurrentSolution, ChisqPrintConstant);
			}
		}

		//  Final Initialization Loop -- save results
		if (Hotsun.errcal > 0)
		{
			// In this version we assume that we cannot invert Chisq matrix so we estimate errors from Chisq matrix -- not its inverse
			CalculateParameterErrors(Hotsun.CurrentSolution);
			Hotsun.errcal = 1;
		}
		else
		{
			SALSABLAS.zrword(Hotsun.perr);
			Hotsun.errcal = 0;
		}
		WriteCoordinates(1, new Hotsun.Hotsun.WriteSolutionSignature(WriteSolution));
		return;

	} // End WriteOutParameters

	public static void WriteCoordinates(int Outputtype, Hotsun.WriteSolutionSignature WriteSolution)
	{
		String outputfiletype = "";
		if (Outputtype == 0)
		{
			outputfiletype = "RESTART";
		}

		//  Output Parameter Values and Errors
		if (!GlobalParameterSet)
		{
			GlobalParameterSet = true;
			MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
		}
		MakeVectorGlobal(Hotsun.perr, Hotsun.UtilityGlobalVector1);

		if (SALSAUtility.MPI_Rank == 0)
		{
			// Set Output File Name and write file with errors
			String labelfilename = ManxcatCentral.config.ReducedVectorOutputFileName;

            ActualOutputFileName = labelfilename;

			WriteSolution.invoke(ActualOutputFileName + outputfiletype, 0, Hotsun.GlobalParameter, Hotsun.UtilityGlobalVector1);

            int slashpos = labelfilename.lastIndexOf(java.io.File.separatorChar);
            String endpart = labelfilename.substring(slashpos + 1);
            labelfilename = labelfilename.substring(0, slashpos + 1);
            ActualOutputFileName = labelfilename + "SIMPLE" + endpart;
			WriteSolution.invoke(ActualOutputFileName + outputfiletype, 1, Hotsun.GlobalParameter, Hotsun.UtilityGlobalVector1);
		}
		return;
	}

	public static void launchQlimits(Hotsun.FindQlimitsSignature FindQlimits)
	{
		if ((Hotsun.numit == 0) || (Hotsun.numit >= ManxcatCentral.IterationtorecalculateQLimits))
		{
			ManxcatCentral.IterationtorecalculateQLimits = Hotsun.numit + Hotsun.QLimitscalculationInterval;
			Hotsun.FullSecondDerivative = false;

			if (ManxcatCentral.config.FullSecondDerivativeOption == 1)
			{
				Hotsun.FullSecondDerivative = true;
			}

			Hotsun.UseDiagonalScaling = Hotsun.UseDiagonalScalinginSolvers;

			if (Hotsun.UseDiagonalScaling)
			{
				SetupDiagonalScaling(Hotsun.CurrentSolution);
			}
			Hotsun.AddMarquardtQDynamically = false;

			SALSAUtility.StartSubTimer(1);
			tangible.RefObject<Double> tempRef_Qhigh = new tangible.RefObject<Double>(Hotsun.Qhigh);
			tangible.RefObject<Double> tempRef_Qlow = new tangible.RefObject<Double>(Hotsun.Qlow);
			tangible.RefObject<Integer> tempRef_EigenvalueIndicator1 = new tangible.RefObject<Integer>(Hotsun.EigenvalueIndicator1);
			tangible.RefObject<Integer> tempRef_EigenvalueIndicator2 = new tangible.RefObject<Integer>(Hotsun.EigenvalueIndicator2);
			FindQlimits.invoke(Hotsun.CurrentSolution, tempRef_Qhigh, tempRef_Qlow, tempRef_EigenvalueIndicator1, tempRef_EigenvalueIndicator2);
			Hotsun.Qhigh = tempRef_Qhigh.argValue;
			Hotsun.Qlow = tempRef_Qlow.argValue;
			Hotsun.EigenvalueIndicator1 = tempRef_EigenvalueIndicator1.argValue;
			Hotsun.EigenvalueIndicator2 = tempRef_EigenvalueIndicator2.argValue;
			Hotsun.teigen = SALSAUtility.StopSubTimer(1);

			Hotsun.IterationforQlimit = Hotsun.CurrentSolution.IterationCalculated;
		}
		return;

	} // End launchQlimits

	public static void CalculateChisqComponents(Hotsun.CalcfgSignature Calcfg, Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct, Desertwind FromSolution, double MultiplicationFactor)
	{
		int savecomponentflag = SALSAUtility.chisqcomponent;
		SALSAUtility.chisqcomponent = 1;
		double chisq1 = MultiplicationFactor * StandaAloneChisq(new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct), FromSolution);
		SALSAUtility.chisqcomponent = 2;
		double chisq2 = MultiplicationFactor * StandaAloneChisq(new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct), FromSolution);
		SALSAUtility.chisqcomponent = 3;
		double chisq3 = MultiplicationFactor * StandaAloneChisq(new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct), FromSolution);
		SALSAUtility.chisqcomponent = 4;
		double chisq4 = MultiplicationFactor * StandaAloneChisq(new Hotsun.Hotsun.CalcfgSignature(Calcfg), new Hotsun.Hotsun.GlobalMatrixVectorProductSignature(GlobalMatrixVectorProduct), FromSolution);
		String componentcomment = " V-V " + String.format("%.3f", chisq1) + " V-F " + String.format("%.3f", chisq2) + " F-V " + String.format("%.3f", chisq3) + " F-F " + String.format("%.3f", chisq4);
		if (SALSAUtility.MPI_Rank == 0)
		{
			SALSAUtility.SALSAPrint(1, componentcomment);
		}
		SALSAUtility.chisqcomponent = savecomponentflag;
		Hotsun.HotsunComment += "\n" + componentcomment;
		return;

	} // End CalculateChisqComponents()

	public static double StandaAloneChisq(Hotsun.CalcfgSignature Calcfg, Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct, Desertwind FromSolution)
	{
		double savechisq = Hotsun.zerocr;
		Hotsun.zerocr = 0.0;
		boolean SaveFullSecondDerivative = Hotsun.FullSecondDerivative;
		Hotsun.FullSecondDerivative = true;
		boolean SaveUseDiagonalScaling = Hotsun.UseDiagonalScaling;
		Hotsun.UseDiagonalScaling = false;
		boolean SaveAddMarquardtQDynamically = Hotsun.AddMarquardtQDynamically;
		Hotsun.AddMarquardtQDynamically = false;
		boolean saveCalcFixedCrossFixed = SALSAUtility.CalcFixedCrossFixed;
		SALSAUtility.CalcFixedCrossFixed = true;
		if (SALSAUtility.StoredDistanceOption == 3)
		{
			SALSAUtility.CalcFixedCrossFixed = false;
		}

		ZeroSolution(Hotsun.SearchSolution1);
		int Numberparms = Hotsun.Number_VectorParameters;
		if (Hotsun.DecomposeParameters)
		{
			Numberparms = SALSAUtility.PointCount_Process;
		}
		SALSABLAS.CopyVector(Hotsun.SearchSolution1.param, FromSolution.param, 0, Numberparms);
		MakeVectorGlobal(Hotsun.SearchSolution1.param, Hotsun.GlobalParameter);

		boolean localviolat = Calcfg.invoke(Hotsun.SearchSolution1);
		double toreturn = Hotsun.zerocr;

		SALSAUtility.CalcFixedCrossFixed = saveCalcFixedCrossFixed;
		Hotsun.FullSecondDerivative = SaveFullSecondDerivative;
		Hotsun.UseDiagonalScaling = SaveUseDiagonalScaling;
		Hotsun.AddMarquardtQDynamically = SaveAddMarquardtQDynamically;
		Hotsun.zerocr = savechisq;
		MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
		return toreturn;

	} // End standalone chisq calculation

	public static void DoDerivTest(Hotsun.CalcfgSignature Calcfg, Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
	{
		if (SALSAUtility.MPI_Size > 1)
		{
			return;
		}
		double savechisq = Hotsun.zerocr;

		int PointtoVary1 = Hotsun.Number_VectorParameters - 25;
		if (PointtoVary1 >= SALSAUtility.NumberVariedPoints)
		{
			PointtoVary1 = SALSAUtility.NumberVariedPoints - 5;
		}
		if (PointtoVary1 < 0)
		{
			PointtoVary1 = 5;
		}

		int IndextoVary1 = 2;
		if (IndextoVary1 > Hotsun.ParameterVectorDimension)
		{
			IndextoVary1 = 0;
		}

		int PointtoVary2 = 25;
		if (PointtoVary2 >= SALSAUtility.NumberVariedPoints)
		{
			PointtoVary2 = 4;
		}

		int IndextoVary2 = 2;
		if (IndextoVary2 > Hotsun.ParameterVectorDimension)
		{
			IndextoVary2 = 0;
		}
		if (Hotsun.npar < 10)
		{
			PointtoVary1 = 6;
			PointtoVary2 = 3;
		}
		int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[PointtoVary1];
		PointtoVary1 = SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex1];
		int OriginalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[PointtoVary2];
		PointtoVary2 = SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex2];

		SALSAUtility.SALSAPrint(1,
                                "Point1 " + PointtoVary1 + "," + IndextoVary1 + " Point2 " + PointtoVary2 + "," +
                                        "" + IndextoVary2);
		Hotsun.doingderivtest = true;
		boolean SaveFullSecondDerivative = Hotsun.FullSecondDerivative;
		Hotsun.FullSecondDerivative = true;
		boolean SaveUseDiagonalScaling = Hotsun.UseDiagonalScaling;
		Hotsun.UseDiagonalScaling = false;
		boolean SaveAddMarquardtQDynamically = Hotsun.AddMarquardtQDynamically;
		Hotsun.AddMarquardtQDynamically = false;

		ZeroSolution(Hotsun.SearchSolution1);
		int Numberparms = Hotsun.Number_VectorParameters;

		if (Hotsun.DecomposeParameters)
		{
			Numberparms = SALSAUtility.PointCount_Global;
		}
		SALSABLAS.CopyVector(Hotsun.SearchSolution1.param, Hotsun.CurrentSolution.param, 0, Numberparms);
		Hotsun.zerocr = 0.0;
		boolean localviolat = Calcfg.invoke(Hotsun.SearchSolution1);
		double startingpoint = Hotsun.zerocr;

		double variableincrement = 0.0001 * Math.abs(Hotsun.CurrentSolution.xshift[PointtoVary1][IndextoVary1]);

		double FirstDerivative1 = 2.0 * Hotsun.SearchSolution1.first[PointtoVary1][IndextoVary1];
		double FirstDerivative2 = 2.0 * Hotsun.SearchSolution1.first[PointtoVary2][IndextoVary2];
		double SecondDerivative11;

		if (Hotsun.fullmatrixset)
		{
			SecondDerivative11 = 2.0 * Hotsun.SearchSolution1.ExactFullMatrix[PointtoVary1][PointtoVary1][IndextoVary1][IndextoVary1];
		}
		else
		{
			SecondDerivative11 = 2.0 * Hotsun.SearchSolution1.ExactDiagonalofMatrix[PointtoVary1][IndextoVary1][IndextoVary1];
		}
		MakeVectorGlobal(Hotsun.SearchSolution1.param, Hotsun.GlobalParameter);

		// Calculate Second Derivative
		SALSABLAS.zrword(Hotsun.UtilityLocalVector1);
		Hotsun.UtilityLocalVector1[PointtoVary1][IndextoVary1] = 1.0;
		MakeVectorGlobal(Hotsun.UtilityLocalVector1, Hotsun.UtilityGlobalVector1);
		GlobalMatrixVectorProduct.invoke(Hotsun.UtilityLocalVector2, Hotsun.SearchSolution1, true, Hotsun.GlobalParameter, Hotsun.UtilityGlobalVector1);
		double SecondDerivative21 = 2.0 * Hotsun.UtilityLocalVector2[PointtoVary2][IndextoVary2];
		double SecondDerivative11A = 2.0 * Hotsun.UtilityLocalVector2[PointtoVary1][IndextoVary1];
		SALSABLAS.zrword(Hotsun.UtilityLocalVector1);
		Hotsun.UtilityLocalVector1[PointtoVary2][IndextoVary2] = 1.0;
		MakeVectorGlobal(Hotsun.UtilityLocalVector1, Hotsun.UtilityGlobalVector1);
		GlobalMatrixVectorProduct.invoke(Hotsun.UtilityLocalVector2, Hotsun.SearchSolution1, true, Hotsun.GlobalParameter, Hotsun.UtilityGlobalVector1);
		double SecondDerivative12 = 2.0 * Hotsun.UtilityLocalVector2[PointtoVary1][IndextoVary1];

		if (Hotsun.npar < 10)
		{
			for (int i1 = 0; i1 < Hotsun.npar; i1++)
			{
				String linetoprint = "Row " + i1 + " ";
				for (int i2 = 0; i2 < Hotsun.npar; i2++)
				{
					double a;
					if (Hotsun.fullmatrixset)
					{
						a = 2.0 * Hotsun.SearchSolution1.ExactFullMatrix[i1][i2][IndextoVary1][IndextoVary1];
					}
					else
					{
						a = 2.0 * Hotsun.SearchSolution1.ExactDiagonalofMatrix[i1][IndextoVary1][IndextoVary1];
					}
					linetoprint += String.format("%.4E", a) + " * ";
				}
				SALSAUtility.SALSAPrint(1, linetoprint);
			}
		}

		double variablestep = 0.0;
		SALSAUtility.SALSAPrint(1, " Variable " + PointtoVary1 + " " + IndextoVary1 + " Value " + String.format("%.4E", Hotsun.SearchSolution1.param[PointtoVary1][IndextoVary1]) + " Increment " + String.format("%.4E", variableincrement));
		for (int step = 0; step < 2; step++)
		{
			Hotsun.SearchSolution1.param[PointtoVary1][IndextoVary1] += variableincrement;
			MakeVectorGlobal(Hotsun.SearchSolution1.param, Hotsun.GlobalParameter);
			Hotsun.zerocr = 0.0;
			ZeroSolution(Hotsun.SearchSolution1);
			localviolat = Calcfg.invoke(Hotsun.SearchSolution1);
			variablestep += variableincrement;
			double change = (Hotsun.zerocr - startingpoint);
			double newzerocr = Hotsun.zerocr;
			double NumericalFirstDerivative1 = change / variablestep;
			double NewFirstDerivative1 = 2.0 * Hotsun.SearchSolution1.first[PointtoVary1][IndextoVary1];
			double NumericalSecondDerivative11 = (NewFirstDerivative1 - FirstDerivative1) / variablestep; // Diagonal Derivative
			double NumericalSecondDerivative21 = (2.0 * Hotsun.SearchSolution1.first[PointtoVary2][IndextoVary2] - FirstDerivative2) / variablestep; // Off-Diagonal Derivative
			double FuncChangeestimate = FirstDerivative1 * variablestep + 0.5 * SecondDerivative11 * variablestep * variablestep;
			SALSAUtility.SALSAPrint(1, step + " step " + String.format("%.4E", variablestep) + " chi " + String.format("%.4E", newzerocr) + " Fact " + String.format("%.4E", change) + " Fest " + String.format("%.4E", FuncChangeestimate) + " D1est " + String.format("%.4E", NumericalFirstDerivative1) + " D1act " + String.format("%.4E", FirstDerivative1) + " S11est " + String.format("%.4E", NumericalSecondDerivative11) + " S11act " + String.format("%.4E", SecondDerivative11) + " S11actA " + String.format("%.4E", SecondDerivative11A) + " S21est " + String.format("%.4E", NumericalSecondDerivative21) + " S21act " + String.format("%.4E", SecondDerivative21) + " S12act " + String.format("%.4E", SecondDerivative12));
		}

		Hotsun.FullSecondDerivative = SaveFullSecondDerivative;
		Hotsun.UseDiagonalScaling = SaveUseDiagonalScaling;
		Hotsun.AddMarquardtQDynamically = SaveAddMarquardtQDynamically;
		Hotsun.doingderivtest = false;
		Hotsun.zerocr = savechisq;
		MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
		return;
	} // End DoDerivTest()


	//  Set up 2D Double MPI Packet
    public static MPI2DDoubleVectorPacket setup2DDoubleMPIPacket() {
        MPI2DDoubleVectorPacket togoVector = new MPI2DDoubleVectorPacket(SALSAUtility.PointCount_Largest, Hotsun.ParameterVectorDimension);
        togoVector.FirstPoint = SALSAUtility.PointStart_Process;
        togoVector.NumberofPoints = SALSAUtility.PointCount_Process;


        if (SALSAUtility.PointCount_Process != SALSAUtility.PointCount_Largest) {
            for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++) {
                togoVector.Marray[SALSAUtility.PointCount_Largest - 1][LocalVectorIndex] = 0.0;
            }
        }
        return togoVector;
    }

    //  Set up 1D String MPI Packet
    public static MPI1DStringVectorPacket setup1DStringMPIPacket() {
        MPI1DStringVectorPacket togoVector = new MPI1DStringVectorPacket(SALSAUtility.PointCount_Largest, Hotsun.ParameterVectorDimension);
        togoVector.FirstPoint = SALSAUtility.PointStart_Process;
        togoVector.NumberofPoints = SALSAUtility.PointCount_Process;

        if (SALSAUtility.PointCount_Process != SALSAUtility.PointCount_Largest) {
            togoVector.Marray[SALSAUtility.PointCount_Largest - 1] = "";
        }
        return togoVector;
    }


	// Globalize distributed local copies with 2D Double
	public static void MakeVectorGlobal(double[][] DistributedVector, double[][] GlobalVector)
	{
		if ((SALSAUtility.MPI_Size <= 1) || (!Hotsun.DecomposeParameters))
		{
			SALSABLAS.CopyVector(GlobalVector, DistributedVector, 0, Hotsun.Number_VectorParameters);
		}
		else
		{
			SALSABLAS.CopyVector(TogoDistributed2DDoubleVector.Marray, DistributedVector, 0, SALSAUtility.PointCount_Process);

			SALSAUtility.StartSubTimer(SALSAUtility.MPISENDRECEIVETiming);

			for (int MPICommunicationSteps = 0; MPICommunicationSteps < SALSAUtility.MPI_Size; MPICommunicationSteps++)
			{
				if (MPICommunicationSteps == SALSAUtility.MPI_Rank)
				{
					MPI2DDoubleVectorPacket.CopyMPI2DDoubleVectorPacket(FromAfar2DDoubleVector, TogoDistributed2DDoubleVector);
				}
				tangible.RefObject<MPI2DDoubleVectorPacket> tempRef_FromAfar2DDoubleVector = new tangible.RefObject<MPI2DDoubleVectorPacket>(FromAfar2DDoubleVector);
				SALSAUtility.MPI_communicator.<MPI2DDoubleVectorPacket>Broadcast(tempRef_FromAfar2DDoubleVector, MPICommunicationSteps);
				FromAfar2DDoubleVector = tempRef_FromAfar2DDoubleVector.argValue;
				SALSABLAS.CopyVector(GlobalVector, FromAfar2DDoubleVector.Marray, FromAfar2DDoubleVector.FirstPoint, FromAfar2DDoubleVector.NumberofPoints);
			} // End loop over MPIrankcount
			SALSAUtility.StopSubTimer(SALSAUtility.MPISENDRECEIVETiming);
		}
	} // End MakeVectorGlobal() 2D Double

	// Globalize distributed local copies with 1D String
	public static void MakeVectorGlobal(String[] DistributedVector, String[] GlobalVector)
	{
		if (SALSAUtility.MPI_Size <= 1)
		{
			SALSABLAS.CopyVector(GlobalVector, DistributedVector, 0, Hotsun.Number_VectorParameters);
		}
		else
		{
			SALSABLAS.CopyVector(TogoDistributed1DStringVector.Marray, DistributedVector, 0, SALSAUtility.PointCount_Process);

			SALSAUtility.StartSubTimer(SALSAUtility.MPISENDRECEIVETiming);

			for (int MPICommunicationSteps = 0; MPICommunicationSteps < SALSAUtility.MPI_Size; MPICommunicationSteps++)
			{
				if (MPICommunicationSteps == SALSAUtility.MPI_Rank)
				{
					MPI1DStringVectorPacket.CopyMPI1DStringVectorPacket(FromAfar1DStringVector, TogoDistributed1DStringVector);
				}
				tangible.RefObject<MPI1DStringVectorPacket> tempRef_FromAfar1DStringVector = new tangible.RefObject<MPI1DStringVectorPacket>(FromAfar1DStringVector);
				SALSAUtility.MPI_communicator.<MPI1DStringVectorPacket>Broadcast(tempRef_FromAfar1DStringVector, MPICommunicationSteps);
				FromAfar1DStringVector = tempRef_FromAfar1DStringVector.argValue;
				SALSABLAS.CopyVector(GlobalVector, FromAfar1DStringVector.Marray, FromAfar1DStringVector.FirstPoint, FromAfar1DStringVector.NumberofPoints);
			} // End loop over MPIrankcount
			SALSAUtility.StopSubTimer(SALSAUtility.MPISENDRECEIVETiming);
		}
	} // End MakeVectorGlobal() with 1D String

	//  Set up GLOBAL arrays Hotsun.diag and Hotsun.sqdginv
	public static void SetupDiagonalScaling(Desertwind Solution)
	{
		if (!Hotsun.DecomposeParameters)
		{
			for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					double tmp;

					if (Hotsun.fullmatrixset)
					{
						if (Hotsun.FullSecondDerivative)
						{
							tmp = Math.abs(Solution.ExactFullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex]);
						}
						else
						{
							tmp = Math.abs(Solution.FullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex]);
						}
					}
					else
					{
						if (Hotsun.FullSecondDerivative)
						{
							tmp = Math.abs(Solution.ExactDiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex]);
						}
						else
						{
							tmp = Math.abs(Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex]);
						}
					}

					Hotsun.diag[LongIndex][LocalVectorIndex] = tmp;
					if (Hotsun.FixedParameter[LongIndex][LocalVectorIndex])
					{
						Hotsun.sqdginv[LongIndex][LocalVectorIndex] = 0.0;
					}
					else
					{
						Hotsun.sqdginv[LongIndex][LocalVectorIndex] = 1.0 / Math.sqrt(tmp);
					}
				}
			}
			return;
		}

		// Parallel Local Calculation of Diagonal Scaling Elements
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LocalToProcessIndex = beginpoint; LocalToProcessIndex < indexlen + beginpoint; LocalToProcessIndex++)
			{
				int GlobalIndex = LocalToProcessIndex + SALSAUtility.PointStart_Process;
				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					double tmp;
					if (Hotsun.fullmatrixset)
					{
						if (Hotsun.FullSecondDerivative)
						{
							tmp = Math.abs(Solution.ExactFullMatrix[LocalToProcessIndex][LocalToProcessIndex][LocalVectorIndex][LocalVectorIndex]);
						}
						else
						{
							tmp = Math.abs(Solution.FullMatrix[LocalToProcessIndex][LocalToProcessIndex][LocalVectorIndex][LocalVectorIndex]);
						}
					}
					else
					{
						if (Hotsun.FullSecondDerivative)
						{
							tmp = Math.abs(Solution.ExactDiagonalofMatrix[LocalToProcessIndex][LocalVectorIndex][LocalVectorIndex]);
						}
						else
						{
							tmp = Math.abs(Solution.DiagonalofMatrix[LocalToProcessIndex][LocalVectorIndex][LocalVectorIndex]);
						}
					}
					TogoDiagVector.Marray[LocalToProcessIndex][LocalVectorIndex] = tmp;
					if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
					{
						TogoSqDgInvVector.Marray[LocalToProcessIndex][LocalVectorIndex] = 0.0;
					}
					else
					{
						TogoSqDgInvVector.Marray[LocalToProcessIndex][LocalVectorIndex] = 1.0 / Math.sqrt(tmp);
					}
				}
			}
		}
	   );

		// Convert into Global arrays
		if (SALSAUtility.MPI_Size <= 1)
		{ // No MPI
			SALSABLAS.CopyVector(Hotsun.diag, TogoDiagVector.Marray, TogoDiagVector.FirstPoint, TogoDiagVector.NumberofPoints);
			SALSABLAS.CopyVector(Hotsun.sqdginv, TogoSqDgInvVector.Marray, TogoSqDgInvVector.FirstPoint, TogoSqDgInvVector.NumberofPoints);
		}
		else
		{
			for (int MPICommunicationSteps = 0; MPICommunicationSteps < SALSAUtility.MPI_Size; MPICommunicationSteps++)
			{
				SALSAUtility.StartSubTimer(SALSAUtility.MPISENDRECEIVETiming);

				if (MPICommunicationSteps == SALSAUtility.MPI_Rank)
				{
					MPI2DDoubleVectorPacket.CopyMPI2DDoubleVectorPacket(FromAfar2DDoubleVector, TogoDiagVector);
				}

				tangible.RefObject<MPI2DDoubleVectorPacket> tempRef_FromAfar2DDoubleVector = new tangible.RefObject<MPI2DDoubleVectorPacket>(FromAfar2DDoubleVector);
				SALSAUtility.MPI_communicator.<MPI2DDoubleVectorPacket>Broadcast(tempRef_FromAfar2DDoubleVector, MPICommunicationSteps);
				FromAfar2DDoubleVector = tempRef_FromAfar2DDoubleVector.argValue;
				SALSABLAS.CopyVector(Hotsun.diag, FromAfar2DDoubleVector.Marray, FromAfar2DDoubleVector.FirstPoint, FromAfar2DDoubleVector.NumberofPoints);

				if (MPICommunicationSteps == SALSAUtility.MPI_Rank)
				{
					MPI2DDoubleVectorPacket.CopyMPI2DDoubleVectorPacket(FromAfar2DDoubleVector, TogoSqDgInvVector);
				}
				tangible.RefObject<MPI2DDoubleVectorPacket> tempRef_FromAfar2DDoubleVector2 = new tangible.RefObject<MPI2DDoubleVectorPacket>(FromAfar2DDoubleVector);
				SALSAUtility.MPI_communicator.<MPI2DDoubleVectorPacket>Broadcast(tempRef_FromAfar2DDoubleVector2, MPICommunicationSteps);
				FromAfar2DDoubleVector = tempRef_FromAfar2DDoubleVector2.argValue;
				SALSABLAS.CopyVector(Hotsun.sqdginv, FromAfar2DDoubleVector.Marray, FromAfar2DDoubleVector.FirstPoint, FromAfar2DDoubleVector.NumberofPoints);
				SALSAUtility.StopSubTimer(SALSAUtility.MPISENDRECEIVETiming);
			} // End loop over MPIrankcount
		}
	} // End SetupDiagonalScaling()

	//  Approximate error estimate in paramters
	public static void CalculateParameterErrors(Desertwind Solution)
	{
		// Parallel Calculation of Parameter Errors
		double NormalizeSQRTChisq = Math.sqrt(Hotsun.zeromn / (Hotsun.ndata - Hotsun.npar));

		if (!Hotsun.DecomposeParameters)
		{
			for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					double tmp;

					if (Hotsun.fullmatrixset)
					{
						tmp = Solution.FullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex];
					}
					else
					{
						tmp = Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex];
					}
					Hotsun.perr[LongIndex][LocalVectorIndex] = NormalizeSQRTChisq / Math.sqrt(tmp);
				}
			}
			return;
		}

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					double tmp;
					if (Hotsun.fullmatrixset)
					{
						tmp = Solution.FullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex];
					}
					else
					{
						tmp = Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex];
					}
					Hotsun.perr[LongIndex][LocalVectorIndex] = NormalizeSQRTChisq / Math.sqrt(Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex]);
				}
			}
		}
	   );

	} // End CalculateParameterErrors()

	public static void AddSubtractMarquardtFactor(Desertwind Solution, double factor)
	{
		if (Hotsun.AddMarquardtQDynamically)
		{
			return;
		}

		if (!Hotsun.DecomposeParameters)
		{
			for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
			{
				int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;

				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					if (Hotsun.fullmatrixset)
					{
						Solution.FullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex] += factor * Hotsun.Q * Hotsun.diag[GlobalIndex][LocalVectorIndex];
						Solution.ExactFullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex] += factor * Hotsun.Q * Hotsun.diag[GlobalIndex][LocalVectorIndex];
					}
					else
					{
						Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex] += factor * Hotsun.Q * Hotsun.diag[GlobalIndex][LocalVectorIndex];
						Solution.ExactDiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex] += factor * Hotsun.Q * Hotsun.diag[GlobalIndex][LocalVectorIndex];
					}
				}
			}
			return;
		}

		// Parallel Addition (factor = 1.0) or Subtraction (factor = -1.0) of Marquardt A
		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;
				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					if (Hotsun.fullmatrixset)
					{
						Solution.FullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex] += factor * Hotsun.Q * Hotsun.diag[GlobalIndex][LocalVectorIndex];
						Solution.ExactFullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex] += factor * Hotsun.Q * Hotsun.diag[GlobalIndex][LocalVectorIndex];
					}
					else
					{
						Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex] += factor * Hotsun.Q * Hotsun.diag[GlobalIndex][LocalVectorIndex];
						Solution.ExactDiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex] += factor * Hotsun.Q * Hotsun.diag[GlobalIndex][LocalVectorIndex];
					}
				}
			}
		}
	   );

	} // End AddinMarquardtFactor

	public static void SteepestDescentSolution(Desertwind Solution, Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
	{
		if (!Hotsun.DecomposeParameters)
		{ // Sequential Calculation of Steepest Descent
			for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					double tmp = Hotsun.sqdginv[LongIndex][LocalVectorIndex];
					Solution.xshift[LongIndex][LocalVectorIndex] = Solution.first[LongIndex][LocalVectorIndex] * tmp;
				}
			}
		}
		else
		{

			// Parallel Calculation of Steepest Descent
						Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
						{
							int indexlen = SALSAUtility.PointsperThread[ThreadNo];
							int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
							for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
							{
								int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;
								for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
								{
									double tmp = Hotsun.sqdginv[GlobalIndex][LocalVectorIndex];
									Solution.xshift[LongIndex][LocalVectorIndex] = Solution.first[LongIndex][LocalVectorIndex] * tmp;
								}
							}
						}
					   );

		}

		// Convert xshift into a Global Vector
		MakeVectorGlobal(Solution.xshift, Hotsun.UtilityGlobalVector1);
		Hotsun.UseDiagonalScaling = false;
		Hotsun.AddMarquardtQDynamically = false;

		//  Calculate Normalization
		// First Calculate ChisqMatrix dotted into xshift
		Hotsun.FullSecondDerivative = false;
		GlobalMatrixVectorProduct.invoke(Hotsun.UtilityLocalVector1, Solution, false, Hotsun.GlobalParameter, Hotsun.UtilityGlobalVector1);
		double Fdotxshift = SALSABLAS.VectorScalarProduct(Solution.xshift, Solution.first);
		double xshiftdotGdotxshift = SALSABLAS.VectorScalarProduct(Solution.xshift, Hotsun.UtilityLocalVector1);
		double factor = Fdotxshift / xshiftdotGdotxshift;
		SALSABLAS.LinearCombineVector(Solution.xshift, factor, Solution.xshift, 0.0, Solution.xshift);
	} // End SteepestDescentSolution()

	//  Find xnorm
	//  If method = 1 take existing xshift
	//  If method = 2 first set xshift from Best - Current Solution
	public static void Findxnorm(int method)
	{
		if (!Hotsun.DecomposeParameters)
		{
			double localnorm = 0.0;

			for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
			{
				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					if (Hotsun.FixedParameter[LongIndex][LocalVectorIndex])
					{
						Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] = 0;
						continue;
					}

					if (method == 2)
					{
						Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] = Hotsun.BestSolution.param[LongIndex][LocalVectorIndex] - Hotsun.CurrentSolution.param[LongIndex][LocalVectorIndex];
					}
					double tmp;

					if (Hotsun.fullmatrixset)
					{
						if (Hotsun.FullSecondDerivative)
						{
							tmp = Math.abs(Hotsun.CurrentSolution.ExactFullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex]);
						}
						else
						{
							tmp = Math.abs(Hotsun.CurrentSolution.FullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex]);
						}
					}
					else
					{
						if (Hotsun.FullSecondDerivative)
						{
							tmp = Math.abs(Hotsun.CurrentSolution.ExactDiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex]);
						}
						else
						{
							tmp = Math.abs(Hotsun.CurrentSolution.DiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex]);
						}
					}
					localnorm += Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] * Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] * tmp;
				}
			}
			Hotsun.xnorm = localnorm;
			return;
		}

		GlobalReductions.FindDoubleSum FindxNorm = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			double localnorm = 0.0;
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
			{
				int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;
				for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
				{
					if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
					{
						Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] = 0;
						continue;
					}
					if (method == 2)
					{
						Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] = Hotsun.BestSolution.param[LongIndex][LocalVectorIndex] - Hotsun.CurrentSolution.param[LongIndex][LocalVectorIndex];
					}
						double tmp;
						if (Hotsun.fullmatrixset)
						{
							if (Hotsun.FullSecondDerivative)
							{
								tmp = Math.abs(Hotsun.CurrentSolution.ExactFullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex]);
							}
							else
							{
								tmp = Math.abs(Hotsun.CurrentSolution.FullMatrix[LongIndex][LongIndex][LocalVectorIndex][LocalVectorIndex]);
							}
						}
						else
						{
							if (Hotsun.FullSecondDerivative)
							{
								tmp = Math.abs(Hotsun.CurrentSolution.ExactDiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex]);
							}
							else
							{
								tmp = Math.abs(Hotsun.CurrentSolution.DiagonalofMatrix[LongIndex][LocalVectorIndex][LocalVectorIndex]);
							}
						}
						localnorm += Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] * Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] * tmp;
				}
			}
			FindxNorm.addapoint(ThreadNo, localnorm);
		}
	   );

		FindxNorm.sumoverthreadsandmpi();
		Hotsun.xnorm = FindxNorm.Total;

	} // End Findxnorm

	// Even if diagonal scaling used in solver, this is not stored in DiagonalofMatrix or first
	public static void FindPredictedChisqChange(Desertwind Solution, Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
	{
		// Calculate xnorm pred1 and pred2
		Findxnorm(1);

		Hotsun.pred1 = 2.0 * SALSABLAS.VectorScalarProduct(Solution.xshift, Solution.first); // Linear Term pred1

		// Convert xshift into a Global Vector
		MakeVectorGlobal(Solution.xshift, Hotsun.UtilityGlobalVector1);
		Hotsun.UseDiagonalScaling = false;
		Hotsun.AddMarquardtQDynamically = false;
		boolean saveflag = Hotsun.FullSecondDerivative;
		Hotsun.FullSecondDerivative = false;
		GlobalMatrixVectorProduct.invoke(Hotsun.UtilityLocalVector1, Solution, false, Hotsun.GlobalParameter, Hotsun.UtilityGlobalVector1);
		Hotsun.pred2 = -SALSABLAS.VectorScalarProduct(Solution.xshift, Hotsun.UtilityLocalVector1);
		Hotsun.FullSecondDerivative = true;
		GlobalMatrixVectorProduct.invoke(Hotsun.UtilityLocalVector1, Solution, true, Hotsun.GlobalParameter, Hotsun.UtilityGlobalVector1);
		double FullPred = -SALSABLAS.VectorScalarProduct(Solution.xshift, Hotsun.UtilityLocalVector1);

		if (saveflag)
		{
			Hotsun.pred3 = Hotsun.pred2;
			Hotsun.pred2 = FullPred;
		}
		else
		{
			Hotsun.pred3 = FullPred;
		}
		Hotsun.FullSecondDerivative = saveflag;
	} // End FindPredictedChisqChange()

	//  CurrentSolution is current configuration wth value derivative
	//  OriginalSolution was previous configuration with value and  derivative
	//  LineFactor * xshift is returned solution that decreases Chisq by an additional Extradecrease in LineIteration iterations

	//  LineSearchMethod    1 No chisq increase. Method should not have been called. Shift is 1
	//  LineSearchMethod    40 + DerivativeSearchMethod Call DerivativeSearch for case where shift < 1; DerivativeSearch is method indicator returned by DerivativeSearch
	//  LineSearchMethod    3 Request to find a shift > 1 but Chisq decreases for shift = 1 so rejected. Return Shift = 1
	//  LineSearchMethod    80 + DerivativeSearchMethod Call DerivativeSearch for case where requested shift > 1 BUT derivative positive at shift of 1; So Shift < 1 returned; DerivativeSearch is method indicator returned by DerivativeSearch
	//  LineSearchMethod    120 + DerivativeSearchMethod Call DerivativeSearch for case where requested shift > 1 AND derivative negative at shift of 1; Shift > 1 possible; DerivativeSearch is method indicator returned by DerivativeSearch

	//  DerivativeSearchMethod 30 + GoldenSectionMethod if Right Derivative < 0 This is a value search between Middle and Right 
	//  DerivativeSearchMethod 20 + GoldenSectionMethod if Left Solution bad This is a value search between Left and Middle
	//  or DerivativeSearchMethod = SimpleDerivativeSearchMethod      search between Middle and Right        
	//  or DerivativeSearchMethod = SimpleDerivativeSearchMethod + 5     search between Right and Middle
	//  or DerivativeSearchMethod = SimpleDerivativeSearchMethod + 10     search between Right and Left

	//  SimpleDerivativeSearchMethod 1 Failure New point higher Take Middle  
	//  SimpleDerivativeSearchMethod 2 Solution is found point with lower chisq and Positive derivative 
	//  SimpleDerivativeSearchMethod 3 Solution is found point possibly with lower chisq -- Cut on position change
	//  SimpleDerivativeSearchMethod 4 Solution is found point possibly with lower chisq -- Cut on value change 
	//  SimpleDerivativeSearchMethod 5 Iteration Cut


	//  GoldenSectionMethod 1   Null Solution handed to routine
	//  GoldenSectionMethod 2   Inconsistent Data   Right Lowest and also Left < Middle
	//  GoldenSectionMethod 3   Inconsistent Data   Left Lowest
	//  GoldenSectionMethod 4   Inconsistent Data   Right Lowest
	//  GoldenSectionMethod 5   Inconsistent Data   Positions inconsistent
	//  GoldenSectionMethod 6   Consistent Data Cut on position values
	//  GoldenSectionMethod 7   Consistent Data Cut on Chisq values
	//  GoldenSectionMethod 8   Consistent Data Iteration Cut 

	public static boolean LineSearch(Desertwind CurrentSolution, Desertwind OriginalSolution, tangible.RefObject<Integer> LineSearchMethod, tangible.RefObject<Double> LineFactor, tangible.RefObject<Double> ExtraDecrease, tangible.RefObject<Integer> LineIterations, Hotsun.CalcfgSignature Calcfg, double EstimatedLineFactor)
	{
		SearchStuff Left = new SearchStuff();
		SearchStuff Middle = new SearchStuff();
		SearchStuff Right = new SearchStuff();
		SearchStuff BestFromLineSearch = new SearchStuff();
		boolean success;
		int LocalIterations = 0;

		LineFactor.argValue = 1.0;
		ExtraDecrease.argValue = 0.0;
		LineIterations.argValue = 0;
		LineSearchMethod.argValue = 0;
		int DerivSearchMethod = 0;

		// BeginningLinePositionSolution; Set to starting point of search (holds param of starting point and value of first derivative there)
		// EndingLinePositionSolution: Set to starting point of search (holds xshift from starting point and value and first derivative of end of line )
		Hotsun.BeginningLinePositionSolution = OriginalSolution;
		Hotsun.EndingLinePositionSolution = CurrentSolution;

		if (EstimatedLineFactor > 0.0)
		{ // Case of LineFactor < 1 from an initial step that failed
			boolean NoChisqIncrease = CurrentSolution.Chisquared <= OriginalSolution.Chisquared;
			tangible.RefObject<Boolean> tempRef_NoChisqIncrease = new tangible.RefObject<Boolean>(NoChisqIncrease);
			SALSAUtility.SynchronizeMPIvariable(tempRef_NoChisqIncrease);
			NoChisqIncrease = tempRef_NoChisqIncrease.argValue;

			if (NoChisqIncrease)
			{
				LineSearchMethod.argValue = 1;
				return false; // Method should not have been called
			}
			double alpha = Math.min(0.5, 2.0 * EstimatedLineFactor);
			CopySolution(Hotsun.SearchSolution1, OriginalSolution);
			CopySolution(Hotsun.SearchSolution3, CurrentSolution);

			tangible.RefObject<SearchStuff> tempRef_Left = new tangible.RefObject<SearchStuff>(Left);
			ExistingCalcDeriv_LineSearch(0.0, Hotsun.SearchSolution1, tempRef_Left);
			Left = tempRef_Left.argValue;
			tangible.RefObject<SearchStuff> tempRef_Middle = new tangible.RefObject<SearchStuff>(Middle);
			NewCalcDeriv_LineSearch(alpha, Hotsun.SearchSolution2, tempRef_Middle, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
		Middle = tempRef_Middle.argValue;
			tangible.RefObject<SearchStuff> tempRef_Right = new tangible.RefObject<SearchStuff>(Right);
			ExistingCalcDeriv_LineSearch(1.0, Hotsun.SearchSolution3, tempRef_Right);
			Right = tempRef_Right.argValue;

			tangible.RefObject<Integer> tempRef_DerivSearchMethod = new tangible.RefObject<Integer>(DerivSearchMethod);
			tangible.RefObject<SearchStuff> tempRef_BestFromLineSearch = new tangible.RefObject<SearchStuff>(BestFromLineSearch);
			tangible.RefObject<Integer> tempRef_LocalIterations = new tangible.RefObject<Integer>(LocalIterations);
			boolean DerivSearch3 = DerivativeSearch(Left.clone(), Middle.clone(), Right.clone(), tempRef_DerivSearchMethod, tempRef_BestFromLineSearch, tempRef_LocalIterations, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
		DerivSearchMethod = tempRef_DerivSearchMethod.argValue;
		BestFromLineSearch = tempRef_BestFromLineSearch.argValue;
		LocalIterations = tempRef_LocalIterations.argValue;
			++LocalIterations;
			LineSearchMethod.argValue = 40 + DerivSearchMethod;
			tangible.RefObject<Boolean> tempRef_DerivSearch3 = new tangible.RefObject<Boolean>(DerivSearch3);
			SALSAUtility.SynchronizeMPIvariable(tempRef_DerivSearch3);
			DerivSearch3 = tempRef_DerivSearch3.argValue;
			success = DerivSearch3;
		}
		else
		{ // Case of LineFactor > 1 from a good step that can be improved
			boolean NoChisqDecrease = CurrentSolution.Chisquared >= OriginalSolution.Chisquared;
			tangible.RefObject<Boolean> tempRef_NoChisqDecrease = new tangible.RefObject<Boolean>(NoChisqDecrease);
			SALSAUtility.SynchronizeMPIvariable(tempRef_NoChisqDecrease);
			NoChisqDecrease = tempRef_NoChisqDecrease.argValue;

			if (NoChisqDecrease)
			{
				LineSearchMethod.argValue = 3;
				return false; // Method should not have been called
			}
			Left.Solution = null;
			CopySolution(Hotsun.SearchSolution1, OriginalSolution);
			CopySolution(Hotsun.SearchSolution2, CurrentSolution);

			tangible.RefObject<SearchStuff> tempRef_Middle2 = new tangible.RefObject<SearchStuff>(Middle);
			ExistingCalcDeriv_LineSearch(0.0, Hotsun.SearchSolution1, tempRef_Middle2);
			Middle = tempRef_Middle2.argValue;
			tangible.RefObject<SearchStuff> tempRef_Right2 = new tangible.RefObject<SearchStuff>(Right);
			ExistingCalcDeriv_LineSearch(1.0, Hotsun.SearchSolution2, tempRef_Right2);
			Right = tempRef_Right2.argValue;

			boolean RightDerivativePositive = Right.deriv > 0;
			tangible.RefObject<Boolean> tempRef_RightDerivativePositive = new tangible.RefObject<Boolean>(RightDerivativePositive);
			SALSAUtility.SynchronizeMPIvariable(tempRef_RightDerivativePositive);
			RightDerivativePositive = tempRef_RightDerivativePositive.argValue;

			if (RightDerivativePositive)
			{ // Inconsistent derivatives

				tangible.RefObject<Integer> tempRef_DerivSearchMethod2 = new tangible.RefObject<Integer>(DerivSearchMethod);
				tangible.RefObject<SearchStuff> tempRef_BestFromLineSearch2 = new tangible.RefObject<SearchStuff>(BestFromLineSearch);
				tangible.RefObject<Integer> tempRef_LocalIterations2 = new tangible.RefObject<Integer>(LocalIterations);
				boolean DerivSearch1 = DerivativeSearch(Left.clone(), Middle.clone(), Right.clone(), tempRef_DerivSearchMethod2, tempRef_BestFromLineSearch2, tempRef_LocalIterations2, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
			DerivSearchMethod = tempRef_DerivSearchMethod2.argValue;
			BestFromLineSearch = tempRef_BestFromLineSearch2.argValue;
			LocalIterations = tempRef_LocalIterations2.argValue;
				LineSearchMethod.argValue = 80 + DerivSearchMethod;
				tangible.RefObject<Boolean> tempRef_DerivSearch1 = new tangible.RefObject<Boolean>(DerivSearch1);
				SALSAUtility.SynchronizeMPIvariable(tempRef_DerivSearch1);
				DerivSearch1 = tempRef_DerivSearch1.argValue;
				success = DerivSearch1;
			}
			else
			{ // Correct Sign of derivative -- explore larger values of alpha
				double alpha = 1.0;

				while (LineIterations.argValue <= 20)
				{ // Iterate over increasing alpha
					++LineIterations.argValue;
					alpha += 1.0;
					Desertwind SaveSolutionReference;

					if (Left.Solution == null)
					{
						SaveSolutionReference = Hotsun.SearchSolution3;
					}
					else
					{
						SaveSolutionReference = Left.Solution;
					}
					Left = Middle.clone();
					Middle = Right.clone();
					tangible.RefObject<SearchStuff> tempRef_Right3 = new tangible.RefObject<SearchStuff>(Right);
					NewCalcDeriv_LineSearch(alpha, SaveSolutionReference, tempRef_Right3, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
				Right = tempRef_Right3.argValue;

					//  Right.deriv > 0 is best case but also stop if value increases at Right even with unexpected negative derivative
					boolean TimetoBreak = (Right.deriv > 0) || (Right.value > Middle.value);
					tangible.RefObject<Boolean> tempRef_TimetoBreak = new tangible.RefObject<Boolean>(TimetoBreak);
					SALSAUtility.SynchronizeMPIvariable(tempRef_TimetoBreak);
					TimetoBreak = tempRef_TimetoBreak.argValue;

					if (TimetoBreak)
					{
						break;
					}
				}
				tangible.RefObject<Integer> tempRef_DerivSearchMethod3 = new tangible.RefObject<Integer>(DerivSearchMethod);
				tangible.RefObject<SearchStuff> tempRef_BestFromLineSearch3 = new tangible.RefObject<SearchStuff>(BestFromLineSearch);
				tangible.RefObject<Integer> tempRef_LocalIterations3 = new tangible.RefObject<Integer>(LocalIterations);
				boolean DerivSearch2 = DerivativeSearch(Left.clone(), Middle.clone(), Right.clone(), tempRef_DerivSearchMethod3, tempRef_BestFromLineSearch3, tempRef_LocalIterations3, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
			DerivSearchMethod = tempRef_DerivSearchMethod3.argValue;
			BestFromLineSearch = tempRef_BestFromLineSearch3.argValue;
			LocalIterations = tempRef_LocalIterations3.argValue;
				LineSearchMethod.argValue = 120 + DerivSearchMethod;
				tangible.RefObject<Boolean> tempRef_DerivSearch2 = new tangible.RefObject<Boolean>(DerivSearch2);
				SALSAUtility.SynchronizeMPIvariable(tempRef_DerivSearch2);
				DerivSearch2 = tempRef_DerivSearch2.argValue;
				success = DerivSearch2;
			}
		} // End case LineFactor > 1.0

		//  Copy Best Solution if needed
		LineIterations.argValue += LocalIterations;

		if (!success)
		{
			MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
			GlobalParameterSet = true;
			return false;
		}

		// Success -- Change Current Solution to Best from Line Search
		LineFactor.argValue = BestFromLineSearch.alpha;
		ExtraDecrease.argValue = CurrentSolution.Chisquared - BestFromLineSearch.value;
		SALSABLAS.LinearCombineVector(BestFromLineSearch.Solution.xshift, LineFactor.argValue, Hotsun.EndingLinePositionSolution.xshift, 0.0, Hotsun.EndingLinePositionSolution.xshift);
		CopySolution(CurrentSolution, BestFromLineSearch.Solution);
		MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
		GlobalParameterSet = true;
		Hotsun.pred1 *= LineFactor.argValue;
		Hotsun.pred2 *= LineFactor.argValue * LineFactor.argValue;
		Hotsun.pred3 *= LineFactor.argValue * LineFactor.argValue;

		return true;
	} // End Line Search

	public static boolean DerivativeSearch(SearchStuff Left, SearchStuff Middle, SearchStuff Right, tangible.RefObject<Integer> SearchMethod, tangible.RefObject<SearchStuff> Best, tangible.RefObject<Integer> Iterations, Hotsun.CalcfgSignature Calcfg)
	{ // Do a derivative Secant search but never stray outside interval by requiring End point derivatives opposite sign
		//  Revert to Value search if problems


		Iterations.argValue = 0;
		Best.argValue = Middle.clone();
		boolean WrongSignofRightDeriv = Right.deriv <= 0;
		tangible.RefObject<Boolean> tempRef_WrongSignofRightDeriv = new tangible.RefObject<Boolean>(WrongSignofRightDeriv);
		SALSAUtility.SynchronizeMPIvariable(tempRef_WrongSignofRightDeriv);
		WrongSignofRightDeriv = tempRef_WrongSignofRightDeriv.argValue;

		if (WrongSignofRightDeriv)
		{ // Inconsistent Derivatives -- Use value Search
			int GoldenSectionMethod = 0;
			tangible.RefObject<Integer> tempRef_GoldenSectionMethod = new tangible.RefObject<Integer>(GoldenSectionMethod);
			boolean GoldenResult = GoldenSectionValueSearch(Left.clone(), Middle.clone(), Right.clone(), tempRef_GoldenSectionMethod, Best.clone(), Iterations, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
		GoldenSectionMethod = tempRef_GoldenSectionMethod.argValue;
			SearchMethod.argValue = 20 + GoldenSectionMethod;
			tangible.RefObject<Boolean> tempRef_GoldenResult = new tangible.RefObject<Boolean>(GoldenResult);
			SALSAUtility.SynchronizeMPIvariable(tempRef_GoldenResult);
			GoldenResult = tempRef_GoldenResult.argValue;
			return GoldenResult;
		}

		int SimpleSearchMethod = 0;
		int LocalIterations = 0;

		boolean UnexpectedSignofMiddleDeriv = Middle.deriv > 0;
		tangible.RefObject<Boolean> tempRef_UnexpectedSignofMiddleDeriv = new tangible.RefObject<Boolean>(UnexpectedSignofMiddleDeriv);
		SALSAUtility.SynchronizeMPIvariable(tempRef_UnexpectedSignofMiddleDeriv);
		UnexpectedSignofMiddleDeriv = tempRef_UnexpectedSignofMiddleDeriv.argValue;

		if (UnexpectedSignofMiddleDeriv)
		{ // Best solution betwenn Left and Middle// First see if Deriv search possible. That needs Left to be defined with a negative derivative
			boolean noleftsolution = (Left.Solution == null);
			boolean derivsearchforbidden = noleftsolution;
			if (!derivsearchforbidden)
			{
				if (Left.deriv > 0)
				{
					derivsearchforbidden = true;
				}
			}
			tangible.RefObject<Boolean> tempRef_derivsearchforbidden = new tangible.RefObject<Boolean>(derivsearchforbidden);
			SALSAUtility.SynchronizeMPIvariable(tempRef_derivsearchforbidden);
			derivsearchforbidden = tempRef_derivsearchforbidden.argValue;

			if (derivsearchforbidden)
			{
				int GoldenSectionMethod = 0;
				tangible.RefObject<Integer> tempRef_GoldenSectionMethod2 = new tangible.RefObject<Integer>(GoldenSectionMethod);
				boolean GoldenResult = GoldenSectionValueSearch(Left.clone(), Middle.clone(), Right.clone(), tempRef_GoldenSectionMethod2, Best.clone(), Iterations, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
			GoldenSectionMethod = tempRef_GoldenSectionMethod2.argValue;
				SearchMethod.argValue = 30 + GoldenSectionMethod;
				tangible.RefObject<Boolean> tempRef_GoldenResult2 = new tangible.RefObject<Boolean>(GoldenResult);
				SALSAUtility.SynchronizeMPIvariable(tempRef_GoldenResult2);
				GoldenResult = tempRef_GoldenResult2.argValue;
				return GoldenResult;
			}
			if (Middle.value > Right.value)
			{
				tangible.RefObject<Integer> tempRef_SimpleSearchMethod = new tangible.RefObject<Integer>(SimpleSearchMethod);
				tangible.RefObject<Integer> tempRef_LocalIterations = new tangible.RefObject<Integer>(LocalIterations);
				boolean result2 = SimpleDerivativeSearch(Left.clone(), Right.clone(), tempRef_SimpleSearchMethod, Best.clone(), tempRef_LocalIterations, new Hotsun.Hotsun.CalcfgSignature(Calcfg), Middle.clone());
			SimpleSearchMethod = tempRef_SimpleSearchMethod.argValue;
			LocalIterations = tempRef_LocalIterations.argValue;
				Iterations.argValue += LocalIterations;
				SearchMethod.argValue = SimpleSearchMethod + 5;
				return result2;
			}
			tangible.RefObject<Integer> tempRef_SimpleSearchMethod2 = new tangible.RefObject<Integer>(SimpleSearchMethod);
			tangible.RefObject<Integer> tempRef_LocalIterations2 = new tangible.RefObject<Integer>(LocalIterations);
			boolean result1 = SimpleDerivativeSearch(Left.clone(), Middle.clone(), tempRef_SimpleSearchMethod2, Best.clone(), tempRef_LocalIterations2, new Hotsun.Hotsun.CalcfgSignature(Calcfg), Right.clone());
		SimpleSearchMethod = tempRef_SimpleSearchMethod2.argValue;
		LocalIterations = tempRef_LocalIterations2.argValue;
			Iterations.argValue += LocalIterations;
			SearchMethod.argValue = SimpleSearchMethod + 10;
			return result1;
		}
		tangible.RefObject<Integer> tempRef_SimpleSearchMethod3 = new tangible.RefObject<Integer>(SimpleSearchMethod);
		tangible.RefObject<Integer> tempRef_LocalIterations3 = new tangible.RefObject<Integer>(LocalIterations);
		boolean result = SimpleDerivativeSearch(Middle.clone(), Right.clone(), tempRef_SimpleSearchMethod3, Best.clone(), tempRef_LocalIterations3, new Hotsun.Hotsun.CalcfgSignature(Calcfg), Left.clone());
	SimpleSearchMethod = tempRef_SimpleSearchMethod3.argValue;
	LocalIterations = tempRef_LocalIterations3.argValue;
		Iterations.argValue += LocalIterations;
		SearchMethod.argValue = SimpleSearchMethod;
		return result;

	} // End DerivativeSearch

	public static boolean SimpleDerivativeSearch(SearchStuff OneSide, SearchStuff OtherSide, tangible.RefObject<Integer> SearchMethod, tangible.RefObject<SearchStuff> Best, tangible.RefObject<Integer> Iterations, Hotsun.CalcfgSignature Calcfg, SearchStuff NotUsed)
	{
		double x2 = OneSide.alpha;
		double x3 = OtherSide.alpha;
		double InitialRange = x3 - x2;
		double TargetDelta = 0.01 * (OtherSide.value - OneSide.value);
		double TargetValueChange = 0.01 * Math.abs(OtherSide.value - OneSide.value);
		SearchMethod.argValue = 5;
		Iterations.argValue = 0;

		while (Iterations.argValue < 10)
		{
			++Iterations.argValue;
			double x4 = x3 - (x3 - x2) * OtherSide.deriv / (OtherSide.deriv - OneSide.deriv);
			SearchStuff NewPoint = new SearchStuff();
			Desertwind NewSolution = matchingSolution(NotUsed.Solution, OneSide.Solution, OtherSide.Solution);
			tangible.RefObject<SearchStuff> tempRef_NewPoint = new tangible.RefObject<SearchStuff>(NewPoint);
			NewCalcDeriv_LineSearch(x4, NewSolution, tempRef_NewPoint, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
		NewPoint = tempRef_NewPoint.argValue;
			double oldtestvalue = Math.min(OneSide.value, OtherSide.value);

			boolean NewPointValueLarger = (NewPoint.value > OneSide.value);
			tangible.RefObject<Boolean> tempRef_NewPointValueLarger = new tangible.RefObject<Boolean>(NewPointValueLarger);
			SALSAUtility.SynchronizeMPIvariable(tempRef_NewPointValueLarger);
			NewPointValueLarger = tempRef_NewPointValueLarger.argValue;

			if (NewPointValueLarger)
			{
				if (NewPoint.deriv < 0)
				{
					Best.argValue = OneSide.clone();
					SearchMethod.argValue = 1;
					return true;
				}
				OtherSide = NewPoint.clone();
			}
			else
			{
				if (NewPoint.deriv > 0)
				{
					Best.argValue = NewPoint.clone();
					SearchMethod.argValue = 2;
					return true;
				}
				OneSide = NewPoint.clone();
			}
			x3 = OtherSide.alpha;
			x2 = OneSide.alpha;

			boolean casetobreak1 = ((OtherSide.value - OneSide.value) < TargetDelta) || ((x3 - x2) < 0.01 * InitialRange);
			tangible.RefObject<Boolean> tempRef_casetobreak1 = new tangible.RefObject<Boolean>(casetobreak1);
			SALSAUtility.SynchronizeMPIvariable(tempRef_casetobreak1);
			casetobreak1 = tempRef_casetobreak1.argValue;

			if (casetobreak1)
			{
				SearchMethod.argValue = 3;
				break;
			}

			double newtestvalue = Math.min(OneSide.value, OtherSide.value);
			boolean casetobreak2 = (newtestvalue > oldtestvalue - TargetValueChange);
			tangible.RefObject<Boolean> tempRef_casetobreak2 = new tangible.RefObject<Boolean>(casetobreak2);
			SALSAUtility.SynchronizeMPIvariable(tempRef_casetobreak2);
			casetobreak2 = tempRef_casetobreak2.argValue;

			if (casetobreak2)
			{
				SearchMethod.argValue = 4;
				break;
			}
		}
		Best.argValue = OneSide.clone();
		return true;
	}


	public static boolean GoldenSectionValueSearch(SearchStuff Left, SearchStuff Middle, SearchStuff Right, tangible.RefObject<Integer> GoldenSectionMethod, tangible.RefObject<SearchStuff> Best, tangible.RefObject<Integer> Iterations, Hotsun.CalcfgSignature Calcfg)
	{ // Do a value search to find minimum assuming three values are set
		//  Return End point if this is minimum (Derivative search will not do this if End Point has positive derivative

		double goldennumber = 1.618033989;
		GoldenSectionMethod.argValue = 0;
		Iterations.argValue = 0;
		Best.argValue = Middle.clone();

		if (Left.Solution == null)
		{
			GoldenSectionMethod.argValue = 1;
			return false;
		}

		//  First Check for Consistency of data
		boolean leftlow = (Left.value <= Middle.value);
		tangible.RefObject<Boolean> tempRef_leftlow = new tangible.RefObject<Boolean>(leftlow);
		SALSAUtility.SynchronizeMPIvariable(tempRef_leftlow);
		leftlow = tempRef_leftlow.argValue;

		if (leftlow)
		{
			boolean leftbiggerthanright = (Left.value > Right.value);
			tangible.RefObject<Boolean> tempRef_leftbiggerthanright = new tangible.RefObject<Boolean>(leftbiggerthanright);
			SALSAUtility.SynchronizeMPIvariable(tempRef_leftbiggerthanright);
			leftbiggerthanright = tempRef_leftbiggerthanright.argValue;

			if (leftbiggerthanright)
			{
				Best.argValue = Right.clone();
				GoldenSectionMethod.argValue = 2;
				return true;
			}
			else
			{
				Best.argValue = Left.clone();
				GoldenSectionMethod.argValue = 3;
				return false;
			}
		}

		boolean rightlow = (Right.value <= Middle.value);
		tangible.RefObject<Boolean> tempRef_rightlow = new tangible.RefObject<Boolean>(rightlow);
		SALSAUtility.SynchronizeMPIvariable(tempRef_rightlow);
		rightlow = tempRef_rightlow.argValue;

		if (rightlow)
		{
			GoldenSectionMethod.argValue = 4;
			Best.argValue = Right.clone();
			return true;
		}

		double x1 = Left.alpha;
		double x2 = Middle.alpha;
		double x3 = Right.alpha;
		boolean inconsistentpositions = (((x2 - x1) <= 0.0) || ((x3 - x2) <= 0.0));
		tangible.RefObject<Boolean> tempRef_inconsistentpositions = new tangible.RefObject<Boolean>(inconsistentpositions);
		SALSAUtility.SynchronizeMPIvariable(tempRef_inconsistentpositions);
		inconsistentpositions = tempRef_inconsistentpositions.argValue;

		if (inconsistentpositions)
		{ // Data Error in positions
			GoldenSectionMethod.argValue = 5;
			return false;
		}
		double InitialRange = x3 - x1;
		double TargetDelta = 0.01 * Math.max(Left.value - Middle.value, Right.value - Middle.value);
		double LastChange = Math.min(Left.value - Middle.value, Right.value - Middle.value);

		//  Data Consistent
		GoldenSectionMethod.argValue = 8;
		while (Iterations.argValue <= 10)
		{
			SearchStuff NewPoint = new SearchStuff();
			Desertwind NewSolution = matchingSolution(Left.Solution, Middle.Solution, Right.Solution);
			double x4;
			SearchStuff save = new SearchStuff();

			boolean leftintervallarger = ((x2 - x1) > (x3 - x2));
			tangible.RefObject<Boolean> tempRef_leftintervallarger = new tangible.RefObject<Boolean>(leftintervallarger);
			SALSAUtility.SynchronizeMPIvariable(tempRef_leftintervallarger);
			leftintervallarger = tempRef_leftintervallarger.argValue;

			if (leftintervallarger)
			{ // Left Interval Bigger -- place New Point to left
				x4 = x2 - (x3 - x2) / goldennumber;
				tangible.RefObject<SearchStuff> tempRef_NewPoint = new tangible.RefObject<SearchStuff>(NewPoint);
				NewCalcDeriv_LineSearch(x4, NewSolution, tempRef_NewPoint, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
			NewPoint = tempRef_NewPoint.argValue;

				boolean NewPointsmallerthanMiddle = (NewPoint.value < Middle.value);
				tangible.RefObject<Boolean> tempRef_NewPointsmallerthanMiddle = new tangible.RefObject<Boolean>(NewPointsmallerthanMiddle);
				SALSAUtility.SynchronizeMPIvariable(tempRef_NewPointsmallerthanMiddle);
				NewPointsmallerthanMiddle = tempRef_NewPointsmallerthanMiddle.argValue;

				if (NewPointsmallerthanMiddle)
				{ // Its x1 x4 x2
					save = Middle.clone();
					Middle = NewPoint.clone();
					Right = save.clone();
				}
				else
				{ // Its x4 x2 x3
					Left = NewPoint.clone();
				}
			}
			else
			{ // // Right Interval Bigger -- place New Point to right
				x4 = x2 + (x2 - x1) / goldennumber;
				tangible.RefObject<SearchStuff> tempRef_NewPoint2 = new tangible.RefObject<SearchStuff>(NewPoint);
				NewCalcDeriv_LineSearch(x4, NewSolution, tempRef_NewPoint2, new Hotsun.Hotsun.CalcfgSignature(Calcfg));
			NewPoint = tempRef_NewPoint2.argValue;

				boolean NewPointsmallerthanMiddle = (NewPoint.value < Middle.value);
				tangible.RefObject<Boolean> tempRef_NewPointsmallerthanMiddle2 = new tangible.RefObject<Boolean>(NewPointsmallerthanMiddle);
				SALSAUtility.SynchronizeMPIvariable(tempRef_NewPointsmallerthanMiddle2);
				NewPointsmallerthanMiddle = tempRef_NewPointsmallerthanMiddle2.argValue;

				if (NewPointsmallerthanMiddle)
				{ // Its x2 x4 x3
					save = Middle.clone();
					Middle = NewPoint.clone();
					Left = save.clone();
				}
				else
				{ // Its x1 x2 x4
					Right = NewPoint.clone();
				}
			}
			++Iterations.argValue;
			x1 = Left.alpha;
			x2 = Middle.alpha;
			x3 = Right.alpha;

			boolean casetobreakonsectionsize = ((x3 - x1) < 0.01 * InitialRange);
			tangible.RefObject<Boolean> tempRef_casetobreakonsectionsize = new tangible.RefObject<Boolean>(casetobreakonsectionsize);
			SALSAUtility.SynchronizeMPIvariable(tempRef_casetobreakonsectionsize);
			casetobreakonsectionsize = tempRef_casetobreakonsectionsize.argValue;

			if (casetobreakonsectionsize)
			{
				GoldenSectionMethod.argValue = 6;
				break;
			}


			double CurrentChange = Math.min(Left.value - Middle.value, Right.value - Middle.value);
			boolean casetobreakonvaluechange = ((CurrentChange <= TargetDelta) && (LastChange <= TargetDelta));
			tangible.RefObject<Boolean> tempRef_casetobreakonvaluechange = new tangible.RefObject<Boolean>(casetobreakonvaluechange);
			SALSAUtility.SynchronizeMPIvariable(tempRef_casetobreakonvaluechange);
			casetobreakonvaluechange = tempRef_casetobreakonvaluechange.argValue;

			if (casetobreakonvaluechange)
			{
				GoldenSectionMethod.argValue = 7;
				break;
			}
			LastChange = CurrentChange;
		}
		Best.argValue = Middle.clone();


		return true;
	} // End GoldenSectionValueSearch

	// Find available Solution from set of pre allocated ones
	public static Desertwind matchingSolution(Desertwind a, Desertwind b, Desertwind c)
	{
		Desertwind[] Possibles = new Desertwind[] {Hotsun.SearchSolution1, Hotsun.SearchSolution2, Hotsun.SearchSolution3, Hotsun.SearchSolution4};
		int[] used = new int[4];

		for (int i = 0; i < 4; i++)
		{
			used[i] = 0;
		}
		used[findsolutionindex(a, Possibles)] = 1;
		used[findsolutionindex(b, Possibles)] = 1;
		used[findsolutionindex(c, Possibles)] = 1;

		for (int i = 0; i < 4; i++)
		{
			if (used[i] == 0)
			{
				return Possibles[i];
			}
		}
		SALSAUtility.printAndThrowRuntimeException("Invalid Solution Set");


	} // End matchingSolution

	//  Find which of pre allocated Solutions correspondds to a particular one stored in trial
	public static int findsolutionindex(Desertwind trial, Desertwind[] Possibles)
	{
		int arraysize = Possibles.length;

		for (int i = 0; i < arraysize; i++)
		{
			if (trial.equals(Possibles[i]))
			{
				return i;
			}
		}

        SALSAUtility.printAndThrowRuntimeException("Invalid Solution Index");
	} // End findsolutionindex

	public static void ExistingCalcDeriv_LineSearch(double alpha, Desertwind SearchSolution, tangible.RefObject<SearchStuff> SearchPoint)
	{ // Calculate value and derivative at a point where Calcfg has been called

		SearchPoint.argValue.alpha = alpha;
		SearchPoint.argValue.value = SearchSolution.Chisquared;
		SearchPoint.argValue.Solution = SearchSolution;
		SearchPoint.argValue.deriv = -2.0 * SALSABLAS.VectorScalarProduct(SearchPoint.argValue.Solution.first, Hotsun.EndingLinePositionSolution.xshift);
	} // End CalcDeriv_LineSearch at a point where Calcfg has been called

	public static void NewCalcDeriv_LineSearch(double alpha, Desertwind NewSolution, tangible.RefObject<SearchStuff> SearchPoint, Hotsun.CalcfgSignature Calcfg)
	{ // Calculate value and derivative at a point where Calcfg has NOT been called

		SearchPoint.argValue.Solution = NewSolution;

		//  Set up parameters for next call to Calcfg
		//  This adds in estimated change to param and zeros first derivative, Second Derivative and Chisq (zerocr)
		SALSABLAS.LinearCombineVector(NewSolution.param, 1.0, Hotsun.BeginningLinePositionSolution.param, -alpha, Hotsun.EndingLinePositionSolution.xshift);
		MakeVectorGlobal(NewSolution.param, Hotsun.GlobalParameter);
		GlobalParameterSet = true;

		Hotsun.idata = 0;
		ZeroSolution(NewSolution);
		//  Call Calcfg to calculate Taylor expansion
		SALSAUtility.StartSubTimer(2);
		violat = Calcfg.invoke(NewSolution);
		SALSAUtility.StopSubTimer(2);

		SearchPoint.argValue.alpha = alpha;
		SearchPoint.argValue.value = NewSolution.Chisquared;
		SearchPoint.argValue.deriv = -2.0 * SALSABLAS.VectorScalarProduct(NewSolution.first, Hotsun.EndingLinePositionSolution.xshift);
		return;

	} // End CalcDeriv_LineSearch at a point where Calcfg has NOT been called

//C# TO JAVA CONVERTER WARNING: Java does not allow user-defined value types. The behavior of this class will differ from the original:
//ORIGINAL LINE: public struct SearchStuff
	public final static class SearchStuff
	{
		public double alpha; // Free Parameter in Line Search
		public double value; // Chisq Value
		public double deriv; // First Derivative wrt alpha
		public Desertwind Solution; // Full Solution Detail -- xshift NOT set typically

		public SearchStuff clone()
		{
			SearchStuff varCopy = new SearchStuff();

			varCopy.alpha = this.alpha;
			varCopy.value = this.value;
			varCopy.deriv = this.deriv;
			varCopy.Solution = this.Solution;

			return varCopy;
		}
	} // End Little struct SearchStuff

	public static void CopySolution(Desertwind Solution2, Desertwind Solution1)
	{
		int Numberparms = Hotsun.Number_VectorParameters;

		if (Hotsun.DecomposeParameters)
		{
			Numberparms = SALSAUtility.PointCount_Process;
		}
		SALSABLAS.CopyVector(Solution2.param, Solution1.param, 0, Numberparms);
		SALSABLAS.CopyVector(Solution2.first, Solution1.first, 0, Numberparms);
		SALSABLAS.CopyVector(Solution2.xshift, Solution1.xshift, 0, Numberparms);

		if (Hotsun.fullmatrixset)
		{
			SALSABLAS.CopyMatrix(Solution2.FullMatrix, Solution1.FullMatrix);
			SALSABLAS.CopyMatrix(Solution2.ExactFullMatrix, Solution1.ExactFullMatrix);
		}
		else
		{
			SALSABLAS.CopyVector(Solution2.DiagonalofMatrix, Solution1.DiagonalofMatrix, 0, Numberparms);
			SALSABLAS.CopyVector(Solution2.ExactDiagonalofMatrix, Solution1.ExactDiagonalofMatrix, 0, Numberparms);
		}

		Solution2.Chisquared = Solution1.Chisquared;
		Solution2.IterationCalculated = Solution1.IterationCalculated;
	}

	public static void SaveBestSolution(Desertwind CurrentSolution, Desertwind BestSolution)
	{
		CopySolution(BestSolution, CurrentSolution);
		Hotsun.zeromn = Hotsun.zerocr;
		Hotsun.Qbest = Hotsun.Q;
		return;

	} // End SaveBestSolution()

	public static void RestoreBestSolution(Desertwind CurrentSolution, Desertwind BestSolution)
	{
		CopySolution(CurrentSolution, BestSolution);
		Hotsun.zerocr = Hotsun.zeromn;

		MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
		GlobalParameterSet = true; // Set Indicator  that Global Parameters are set

	} // End RestoreBestSolution()

	public static void ZeroSolution(Desertwind CurrentSolution)
	{ //  xshift and param are NOT zeroed
		SALSABLAS.zrword(CurrentSolution.first);

		if (Hotsun.fullmatrixset)
		{
			SALSABLAS.zrword(CurrentSolution.FullMatrix);
			SALSABLAS.zrword(CurrentSolution.ExactFullMatrix);
		}
		else
		{
			SALSABLAS.zrword(CurrentSolution.DiagonalofMatrix);
			SALSABLAS.zrword(CurrentSolution.ExactDiagonalofMatrix);
		}
		CurrentSolution.Chisquared = 0.0;
		CurrentSolution.IterationCalculated = Hotsun.numit + 1;
		Hotsun.zerocr = 0.0;
	} // End ZeroSolution()

} // End Class ManxcatCentral // End namespace Manxcat