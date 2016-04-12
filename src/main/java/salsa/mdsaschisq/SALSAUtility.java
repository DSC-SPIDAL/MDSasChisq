package salsa.mdsaschisq;

import com.google.common.base.Stopwatch;
import com.google.common.base.Strings;
import mpi.Intracomm;
import mpi.MPIException;
import salsa.mpi.MpiOps;

import java.io.*;
import java.nio.ByteOrder;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

public class SALSAUtility {
    public static int checkerboard = 0; // If 0 store full matrix, If 1 lower triangular, If 2 use load balanced checkerboard pattern
    public static Intracomm MPI_communicator = null; //MPI communicator
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
    public static int PointCount_Transforming_File = 0;
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

    // TODO - distance type - short
    public static short[][] PointDistances; // Point Distances
    public static int MatrixBreakFactor = 1;
    public static int[] ArrayDivision;
    public static int[] Indexsubtraction;
    public static int DivisionSize;

    public static Date startTime;
    public static Date endTime;
    public static Stopwatch mainTimer;
    public static Stopwatch PreciseTimer; //    Hold Precise Timing
    public static int NumberofSubTimings = 0; // Number of subtimings
    public static Stopwatch[] SubTimers; // Timing Objects
    public static long mainDuration = 0; // duration measured by mainTimer in milliseconds
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
    public static ArrayList<String> CosmicOutput = new ArrayList<>(1000); // Monitoring Output
    public static boolean ConsoleDebugOutput = true; // If true send Monitoring output to console
    public static int DebugPrintOption = 2; // Control Printing (= 0 None, ==1 Summary, = 2 Full)

    public static MpiOps mpiOps;
    public static int dataTypeSize;
    public static ByteOrder endianness;


    // TODO - Debugs
    public static void debugPrintCameHere(String args){
//        debugPrintCameHere(args, true);
    }
    public static void debugPrintCameHere(String args, int onRank) {
//        if (onRank != -2 && MPI_Rank != onRank) return;
//        System.out.println("***Rank=" + MPI_Rank + " Came here " + (!Strings.isNullOrEmpty(args) ? args : ""));
    }
    public static void debugPrintCameHere(String args, boolean onlyRankZero){
//        debugPrintCameHere(args, onlyRankZero ? 0 : -2);
    }

    public static void SetupDistanceWeights() {
        // Set up Distance Weights
        if (!Strings.isNullOrEmpty(ManxcatCentral.config.DistanceWeightsCuts)) {
            Pattern pattern = Pattern.compile("[,]");
            String[] weightCuts = pattern.split(ManxcatCentral.config.DistanceWeightsCuts.trim());
            NumberDistanceWeightCuts = weightCuts.length;
            DistanceWeightCuts = new double[NumberDistanceWeightCuts];
            ActualWeightCuts = new double[NumberDistanceWeightCuts + 1];
            for (int i = 0; i < NumberDistanceWeightCuts; i++) {
                DistanceWeightCuts[i] = Double.parseDouble(weightCuts[i]);
            }
        }
    }

    public static void printException(Exception e) {
        System.out.println("SALSA Error " + e.getMessage());
    }

    public static void printAndThrowRuntimeException(RuntimeException e) {
        System.out.println("SALSA Error " + e.getMessage());
        throw e;
    }

    public static void printAndThrowRuntimeException(String message) {
        System.out.println("SALSA Error " + message);
        throw new RuntimeException(message);
    } // end printAndThrowRuntimeException

    // printOption = 0 Essential Printout
    // printOption = 1 Summary Printout
    // printOption = 2 Only if full print out requested
    public static void SALSAPrint(int printOption, String stuffToPrint) {
        if (MPI_Rank != 0 || DebugPrintOption < printOption) {
            return;
        }
        CosmicOutput.add(stuffToPrint);

        if (ConsoleDebugOutput) {
            System.out.println(stuffToPrint);
        }
    }

    public static int SALSAGracefulend(String directory) {
        int readInstruction = -1;
        Path file = Paths.get(directory, ManxcatConstants.GRACEFULEND_FILE_NAME);

        if (SALSAUtility.MPI_Rank == 0 && Files.isRegularFile(file)) {
            try (BufferedReader br = Files.newBufferedReader(file, Charset.defaultCharset())) {
                Pattern pattern = Pattern.compile("[\t ]");
                // Read first line of file
                String inputLineStr = br.readLine();
                if (inputLineStr != null && inputLineStr.length() > 0) {
                    String[] split = pattern.split(inputLineStr.trim());
                    if ((split.length > 0) && (split[0].length() > 0)) {
                        readInstruction = Integer.parseInt(split[0]);
                        SALSAPrint(0, " Termination Instruction " + readInstruction);
                    }
                }
            } catch (IOException e) {
                SALSAUtility.printAndThrowRuntimeException(" Failure reading " + file + " " + e.toString());
            }
        }
        return readInstruction;
    }

    public static void SALSAStatus(String directory, String Message) {
        if (SALSAUtility.MPI_Rank != 0) {
            return;
        }
        Path file = Paths.get(directory, ManxcatConstants.STATUS_FILE_NAME);

        try (PrintWriter writer = new PrintWriter(file.toString())) {
            writer.println(Message);
        } catch (FileNotFoundException e) {
            SALSAUtility.printAndThrowRuntimeException("Failed writing status data" + e);
        }
    }

    public static void InitializeTiming(int InputNumberofTimers) {
        NumberofSubTimings = InputNumberofTimers;
        SubTimers = new Stopwatch[NumberofSubTimings]; // Timing Objects

        SubDurations = new double[NumberofSubTimings]; // Hold partial timing
        SubTimingEnable = new boolean[NumberofSubTimings];
        SubTimingNames = new String[NumberofSubTimings];

        for (int itimer = 0; itimer < NumberofSubTimings; itimer++) {
            SubTimers[itimer] = Stopwatch.createUnstarted();
            SubDurations[itimer] = 0.0;
            SubTimingEnable[itimer] = true;
        }

        PreciseTimer = Stopwatch.createStarted();
        mainTimer = Stopwatch.createStarted();
        startTime = new Date();
    } // End InitializeTiming

    public static void SetUpSubTimer(int TimingIndex, String TimingLabel) {
        if (TimingIndex >= NumberofSubTimings) {
            printAndThrowRuntimeException("Error in Timing Index " + TimingIndex + " Max " + NumberofSubTimings);
            return;
        }
        SubTimingNames[TimingIndex] = TimingLabel;
        SubTimingEnable[TimingIndex] = true;
        SubDurations[TimingIndex] = 0.0;
    } // End SetUpSubTimer

    public static void SetUpMPISubTimers(int StartTimingIndex, String MPISetLabel) {
        if ((StartTimingIndex + 4) >= NumberofSubTimings) {
            printAndThrowRuntimeException(
                    "Error in  MPI Timing Index " + (StartTimingIndex + 2) + " Max " + NumberofSubTimings);
            return;
        }
        SetUpSubTimer(StartTimingIndex, MPISetLabel + "MPI Reduce");
        SetUpSubTimer(StartTimingIndex + 1, MPISetLabel + "MPI SR");
        SetUpSubTimer(StartTimingIndex + 2, MPISetLabel + "MPI Bcast");
        SetUpSubTimer(StartTimingIndex + 3, MPISetLabel + "MPI Global Reductions");
        SetUpSubTimer(StartTimingIndex + 4, MPISetLabel + "Thread Global Reductions");

        if ((!MPISetLabel.equals("")) && (!MPISetLabel.equals("Lib "))) {
            return;
        }
        MPIREDUCETiming = StartTimingIndex;
        MPISENDRECEIVETiming = StartTimingIndex + 1;
        MPIBROADCASTTiming = StartTimingIndex + 2;
        MPIREDUCETiming1 = StartTimingIndex + 3;
        ThreadTiming = StartTimingIndex + 4;

    } // End SetUpMPISubTimers

    public static void InterimTiming() {
        PreciseTimer.stop();
        HPDuration += PreciseTimer.elapsed(TimeUnit.MICROSECONDS) * 0.001;
        PreciseTimer.reset();
        PreciseTimer.start();
    } // end EndTiming

    public static void EndTiming(boolean fullStop) {
        PreciseTimer.stop();
        HPDuration += PreciseTimer.elapsed(TimeUnit.MICROSECONDS) * 0.001;
        PreciseTimer.reset();
        if (fullStop){
            mainTimer.stop();
            mainDuration = mainTimer.elapsed(TimeUnit.MILLISECONDS);
            mainTimer.reset();
        } else {
            PreciseTimer.start();
        }
        endTime = new Date();
    } // end EndTiming

    public static void StartSubTimer(int TimingIndex) {
        if (TimingIndex < 0) {
            return;
        }

        if (SubTimingEnable[TimingIndex]) {
            SubTimers[TimingIndex].start();
        }
    } // End StartSubTimer

    public static double StopSubTimer(int TimingIndex) {
        if (TimingIndex < 0) {
            return 0.0;
        }

        if (SubTimingEnable[TimingIndex]) {
            SubTimers[TimingIndex].stop();
            double elapsed = SubTimers[TimingIndex].elapsed(TimeUnit.MICROSECONDS) * 0.001;
            SubDurations[TimingIndex] += elapsed;
            SubTimers[TimingIndex].reset();
            return elapsed;
        }
        return 0.0;

    } // End StopSubTimer


    public static int synchronizeMPIVariable(int sync) throws MPIException {
        if (SALSAUtility.MPI_Size > 1) {
            SALSAUtility.StartSubTimer(SALSAUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - int
            sync = SALSAUtility.mpiOps.broadcast(sync, 0);
            SALSAUtility.StopSubTimer(SALSAUtility.MPIBROADCASTTiming);
        }
        return sync;
    }

    public static boolean synchronizeMPIVariable(boolean sync) throws MPIException {
        if (SALSAUtility.MPI_Size > 1) {
            SALSAUtility.StartSubTimer(SALSAUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - boolean
            sync = SALSAUtility.mpiOps.broadcast(sync, 0);
            SALSAUtility.StopSubTimer(SALSAUtility.MPIBROADCASTTiming);
        }
        return sync;
    }

    // Write performance results of process 0 into a file
    public static void WriteTiming_Cluster(String fileName, String RunLabel, int RunNumber, String DataFileName,
                                           String ProcessorName) {
        if (Strings.isNullOrEmpty(fileName)) return;
        File file = new File(fileName);
        boolean isNewFile = !file.isFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            if (isNewFile) {
                // Formats must match those in following data write
                writer.println(String.format(
                        "%1$-15s%2$-15s%3$-8s%4$-12s%5$-8s%6$-12s%7$-12s%8$-16s%9$-12s%10$-25s%11$-25s%12$-25s",
                        "Duration(ms)", "HPDuration(ms)", "Thread#", "MPIperNode", "Node#", "Pt#/Process", "Pt#/Global",
                        "Run Label", "Run Number", "DataFileName", "CurrentTime", "ProcessorName"));

                for (int subtimer = 0; subtimer < SALSAUtility.NumberofSubTimings; subtimer++) {
                    if (!SALSAUtility.SubTimingEnable[subtimer]) {
                        continue;
                    }

                    writer.println(String.format("%1$-15s%2$-8s", SALSAUtility.SubTimingNames[subtimer],
                            "Ratio")); // Subset Time    Name
                }
                writer.println(" "); // Blank line
                writer.println(" "); // Blank line
            }

            // Format syntax (Variable Number , Spaces with - left justified and no dash right aligned
            // Processor Used -  Current Time -  Name of input data file -  Run Number -  Label for Run -
            // Global pointsn -  Local points -  Node# aquired -  Process# per node -  Thread# per process -
            // High performance timer -  Total time
            writer.println(String.format(
                    "%1$-15s%2$-15s%3$-8s%4$-12s%5$-8s%6$-12s%7$-12s%8$-16s%9$-12s%10$-25s%11$-25s%12$-25s",
                    Math.round(mainDuration), Math.round(SALSAUtility.HPDuration * Math.pow(10, 0)) / Math.pow(10, 0),
                    SALSAUtility.ThreadCount, SALSAUtility.MPIperNodeCount, SALSAUtility.NodeCount,
                    SALSAUtility.PointCount_Process, SALSAUtility.PointCount_Global, RunLabel, RunNumber, DataFileName,
                    startTime, ProcessorName));

            for (int subtimer = 0; subtimer < SALSAUtility.NumberofSubTimings; subtimer++) {
                if (!SALSAUtility.SubTimingEnable[subtimer]) {
                    continue;
                }

                double SubTime = SALSAUtility.SubDurations[subtimer];

                double SubRatio = SubTime / (HPDuration);
                // Subset Time
                writer.println(String.format("%1$-15s%2$-8s", Math.round(SubTime), String.format("%.4f", SubRatio)));
            }
            writer.println(" "); // Blank line
        } catch (FileNotFoundException e) {
            printAndThrowRuntimeException("Failed writing data" + e);
        }
    }

    public static String formatElapsedMillis(long elapsed) {
        String format = "%dd:%02dH:%02dM:%02dS:%03dmS";
        short millis = (short) (elapsed % (1000.0));
        elapsed = (elapsed - millis) / 1000; // remaining elapsed in seconds
        byte seconds = (byte) (elapsed % 60.0);
        elapsed = (elapsed - seconds) / 60; // remaining elapsed in minutes
        byte minutes = (byte) (elapsed % 60.0);
        elapsed = (elapsed - minutes) / 60; // remaining elapsed in hours
        byte hours = (byte) (elapsed % 24.0);
        long days = (elapsed - hours) / 24; // remaining elapsed in days
        return String.format(format, days, hours, minutes, seconds, millis);
    }
} // End class SalsaUtility