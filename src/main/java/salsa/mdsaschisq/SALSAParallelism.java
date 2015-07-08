package salsa.mdsaschisq;

import mpi.MPI;
import mpi.MPIException;
import salsa.common.DistanceReader;
import salsa.mpi.MpiOps;

import java.io.IOException;

import static edu.rice.hj.Module0.finalizeHabanero;
import static edu.rice.hj.Module0.initializeHabanero;

/**
 * Set Parallelism related Parameters
 */
public class SALSAParallelism {

    //  Distance related variables
    public static short[][] buffer;

    public static void SetupParallelism(String[] args) throws MPIException {
        //  Set up MPI
        MPI.Init(args);
        SALSAUtility.MPI_communicator = MPI.COMM_WORLD; //initializing MPI world communicator
        SALSAUtility.MPI_Rank = SALSAUtility.MPI_communicator.getRank(); // Rank of this process
        SALSAUtility.MPI_Size = SALSAUtility.MPI_communicator.getSize(); // Number of MPI Processes
        SALSAUtility.mpiOps = new MpiOps(SALSAUtility.MPI_communicator);
        // Set up MPI
        SALSAUtility.MPIperNodeCount =
            SALSAUtility.MPI_Size / SALSAUtility.NodeCount;
        ManxcatCentral.config.MPIperNodeCount = SALSAUtility.MPIperNodeCount;

        if ((SALSAUtility.MPIperNodeCount * SALSAUtility.NodeCount) !=
            SALSAUtility.MPI_Size) {
            SALSAUtility.printAndThrowRuntimeException(
                "Inconsistent MPI counts Nodes " + SALSAUtility.NodeCount +
                " Size " + SALSAUtility.MPI_Size);
        }

        SALSAUtility.ParallelPattern = "Machine:" + MPI.getProcessorName() +
                                       " " +
                                       SALSAUtility.ThreadCount + "x" + SALSAUtility.MPIperNodeCount +
                                       "x" +
                                       SALSAUtility.NodeCount;
        if (SALSAUtility.MPI_Rank == 0) {
            // TODO - distance type - short
            SALSAUtility.SALSAPrint(
                0, " Distance Data Type: " +
                   (ManxcatCentral.config.dataTypeSize == 2 ? "short"
                                                            : ManxcatCentral.config.dataTypeSize));
            SALSAUtility.SALSAPrint(0, SALSAUtility.ParallelPattern);
        }
    }

    public static void TearDownParallelism() throws MPIException {
        // End MPI
        MPI.Finalize();
    }

    public static void SetParallelDecomposition() {
        //	First divide points among processes
        Range[] processRanges = RangePartitioner.Partition(SALSAUtility.PointCount_Global, SALSAUtility.MPI_Size);
        Range processRange = processRanges[SALSAUtility.MPI_Rank]; // The answer for this process

        SALSAUtility.PointStart_Process = processRange.getStartIndex();
        SALSAUtility.PointCount_Process = processRange.getLength();
        SALSAUtility.VariedPointStart_Process = SALSAUtility.PointStart_Process;
        SALSAUtility.VariedPointCount_Process = SALSAUtility.PointCount_Process;
        SALSAUtility.PointCount_Largest = Integer.MIN_VALUE;

        for (Range r : processRanges) {
            SALSAUtility.PointCount_Largest = Math.max(r.getLength(), SALSAUtility.PointCount_Largest);
        }

        // We need points per process for all processes as used by load balancing algorithm and then to reorder
        // distance matrix consistently
        SALSAUtility.PointsperProcess = new int[SALSAUtility.MPI_Size];
        SALSAUtility.PointsperThreadperProcess = new int[SALSAUtility.MPI_Size][];

        for (int i = 0; i < SALSAUtility.MPI_Size; i++) {
            SALSAUtility.PointsperProcess[i] = processRanges[i].getLength();
            Range[] threadprocessRanges = RangePartitioner.Partition(processRanges[i], SALSAUtility.ThreadCount);
            SALSAUtility.PointsperThreadperProcess[i] = new int[SALSAUtility.ThreadCount];
            for (int j = 0; j < SALSAUtility.ThreadCount; j++) {
                SALSAUtility.PointsperThreadperProcess[i][j] = threadprocessRanges[j].getLength();
            }
        }

        //	Now divide points among threads for this process
        Range[] threadRanges = RangePartitioner
                .Partition(processRanges[SALSAUtility.MPI_Rank], SALSAUtility.ThreadCount);
        SALSAUtility.PointsperThread = new int[SALSAUtility.ThreadCount];
        SALSAUtility.StartPointperThread = new int[SALSAUtility.ThreadCount];

        for (int j = 0; j < SALSAUtility.ThreadCount; j++) {
            SALSAUtility.PointsperThread[j] = threadRanges[j].getLength();
            SALSAUtility.StartPointperThread[j] = threadRanges[j].getStartIndex();
        }
    }

    //  Read Distance Data
    public static void ReadDataFromFile(String fname) throws IOException, MPIException {
        if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0)) {
            SALSAUtility.SALSAPrint(1, "Starting to read data: " + " Distance Cut " + String
                    .format("%.3f", SALSAUtility.DistanceCut));
        }
        double countremoveddistances = 0.0;
        double counttotaldistances = 0.0;

        SALSAUtility.StoredDistanceOption = Math.max(2, SALSAUtility.StoredDistanceOption);
        if (SALSAUtility.PointCount_Global == SALSAUtility.NumberOriginalPoints) {
            SALSAUtility.DiskDistanceOption = Math.max(2, SALSAUtility.DiskDistanceOption);
        }

        // Remove unsupported options
        if (SALSAUtility.DiskDistanceOption == 3) {
            SALSAUtility.StoredDistanceOption = 3;
        }
        if (SALSAUtility.StoredDistanceOption == 3) {
            SALSAUtility.CalcFixedCrossFixed = false;
        }

        // Set sizes of matrix on disk
        int colcount = SALSAUtility.PointCount_Global;
        if (SALSAUtility.DiskDistanceOption == 1) {
            colcount = SALSAUtility.NumberOriginalPoints;
        }

        boolean Oneread = true;
        if (SALSAUtility.StoredDistanceOption != SALSAUtility.DiskDistanceOption) {
            Oneread = false;
        }
        if (SALSAUtility.Usedreordered) {
            Oneread = false;
        }

        if (Oneread) {
            SALSAUtility.PointDistances = DistanceReader.readRowRangeAsShort(fname, SALSAUtility.PointStart_Process,
                                                                             SALSAUtility.PointStart_Process +
                                                                                     SALSAUtility.PointCount_Process
                                                                                     - 1,
                                                                             colcount,
                                                                             SALSAUtility.endianness);
            int numberofcolumns = SALSAUtility.PointCount_Global;
            for (int GlobalPointIndex = SALSAUtility.PointStart_Process; GlobalPointIndex < SALSAUtility
                    .PointStart_Process + SALSAUtility.PointCount_Process; GlobalPointIndex++) {
                int rowindex = GlobalPointIndex;
                if (SALSAUtility.StoredDistanceOption == 2) {
                    rowindex = rowindex - SALSAUtility.PointStart_Process;
                }
                if (SALSAUtility.StoredDistanceOption == 3) {
                    int originalpoint = SALSAUtility.UsedPointtoOriginalPointMap[rowindex];
                    int variedpoint = SALSAUtility.OriginalPointDisposition[originalpoint] - SALSAUtility.SALSASHIFT;
                    if (variedpoint < 0) {
                        SALSAUtility.printAndThrowRuntimeException(
                                " Illegal Distance Request Used Point " + rowindex + " Original " + originalpoint);

                    }
                    rowindex = variedpoint - SALSAUtility.VariedPointStart_Process;
                }
                for (int columnindex = 0; columnindex < numberofcolumns; columnindex++) {
                    short temp = SALSAUtility.PointDistances[rowindex][columnindex];
                    counttotaldistances = counttotaldistances + 1.0;

                    if (SALSAUtility.DistanceCut > 0.0) {
                        double distancevalue = (temp / (Short.MAX_VALUE * 1.0));
                        if (distancevalue > SALSAUtility.DistanceCut) {
                            SALSAUtility.PointDistances[rowindex][columnindex] = Short.MAX_VALUE;
                            countremoveddistances = countremoveddistances + 1.0;
                        }
                    }
                }
            }
        } else {
            int colsread = SALSAUtility.PointCount_Global;
            if (SALSAUtility.DiskDistanceOption == 1) {
                colsread = SALSAUtility.NumberOriginalPoints;
            }

            if (SALSAUtility.StoredDistanceOption == 2) {
                SALSAUtility.PointDistances = new short[SALSAUtility.PointCount_Process][];
            }
            if (SALSAUtility.StoredDistanceOption == 3) {
                SALSAUtility.PointDistances = new short[SALSAUtility.VariedPointCount_Process][];
            }

            for (int GlobalPointIndex = SALSAUtility.PointStart_Process; GlobalPointIndex < SALSAUtility
                    .PointStart_Process + SALSAUtility.PointCount_Process; GlobalPointIndex++) {
                int rowtostore = GlobalPointIndex;
                if (SALSAUtility.StoredDistanceOption == 2) {
                    rowtostore = GlobalPointIndex - SALSAUtility.PointStart_Process;
                }
                if (SALSAUtility.StoredDistanceOption == 3) {
                    int OriginalIndex = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
                    rowtostore = SALSAUtility.OriginalPointDisposition[OriginalIndex] - SALSAUtility.SALSASHIFT;
                    if (rowtostore < 0) {
                        continue;
                    }
                    rowtostore = rowtostore - SALSAUtility.VariedPointStart_Process;
                }


                int rowtoread = -1;
                if (SALSAUtility.DiskDistanceOption == 1) {
                    rowtoread = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
                }
                if (SALSAUtility.DiskDistanceOption == 2) {
                    rowtoread = SALSAUtility.ActualtoNaiveUsedOrder[GlobalPointIndex];
                }
                if (SALSAUtility.DiskDistanceOption == 3) {
                    int OriginalIndex = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
                    rowtoread = SALSAUtility.OriginalPointDisposition[OriginalIndex] - SALSAUtility.SALSASHIFT;
                    if (rowtoread < 0) {
                        continue;
                    }
                    rowtoread = SALSAUtility.ActualtoNaiveUsedOrder[rowtoread];
                }
                SALSAUtility.PointDistances[rowtostore] = new short[colcount];
                buffer = DistanceReader.readRowRangeAsShort(fname, rowtoread, rowtoread, colcount,
                                                            SALSAUtility.endianness);
                // Store buffer in PointDistances
                for (int colIndex = 0; colIndex < colsread; colIndex++) {
                    int coltostore = colIndex;
                    if (SALSAUtility.DiskDistanceOption == 1) {
                        coltostore = SALSAUtility.OriginalPointtoUsedPointMap[colIndex];
                        if (coltostore < 0) {
                            continue;
                        }
                    } else if (SALSAUtility.DiskDistanceOption > 1) {
                        coltostore = SALSAUtility.NaivetoActualUsedOrder[colIndex];
                    }
                    SALSAUtility.PointDistances[rowtostore][coltostore] = buffer[0][colIndex];
                    counttotaldistances = counttotaldistances + 1.0;

                    if (SALSAUtility.DistanceCut > 0.0) {
                        short temp = SALSAUtility.PointDistances[rowtostore][coltostore];
                        double distancevalue = (temp / (Double.MAX_VALUE * 1.0));
                        if (distancevalue > SALSAUtility.DistanceCut) {
                            SALSAUtility.PointDistances[rowtostore][coltostore] = Short.MAX_VALUE;
                            countremoveddistances = countremoveddistances + 1.0;
                        }
                    }
                }
            }
        }
        SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming);
        counttotaldistances = SALSAUtility.mpiOps.allReduce(counttotaldistances, MPI.SUM);
        countremoveddistances = SALSAUtility.mpiOps.allReduce(countremoveddistances, MPI.SUM);
        SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming);
        double fractionLeft = 1.0 - countremoveddistances / counttotaldistances;
        if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0)) {
            SALSAUtility.SALSAPrint(1, "Total Distances " + String
                    .format("%.0f", counttotaldistances) + " Distances Removed on Input " + String
                    .format("%.0f", countremoveddistances) + " Fraction Left " + String.format("%.5f", fractionLeft));
        }
    }

    /**
     * Gets the distance value corresponding to row and column of Used points
     *
     * @param row row number from Used points
     * @param col col number from Used points
     * @return -1 if distance undefined, else the distance value corresponding to row and col
     */
    public static double getDistanceValue(int row, int col) {
        if (SALSAUtility.StoredDistanceOption == 2) {
            row = row - SALSAUtility.PointStart_Process;
        }
        if (SALSAUtility.StoredDistanceOption == 3) {
            int originalPoint = SALSAUtility.UsedPointtoOriginalPointMap[row];
            int variedPoint = SALSAUtility.OriginalPointDisposition[originalPoint] - SALSAUtility.SALSASHIFT;
            if (variedPoint < 0) {
                SALSAUtility.printAndThrowRuntimeException(
                        " Illegal Distance Request Used Point " + row + " Original " + originalPoint);

            }
            row = variedPoint - SALSAUtility.VariedPointStart_Process;
        }

        short Temp = SALSAUtility.PointDistances[row][col];

        // Return -1.0 if input distance < 0 or if equals TDistance.MaxValue
        if (Temp == Short.MAX_VALUE) {
            if (SALSAUtility.DistanceCut > 0.0) {
                if (SALSAUtility.UndefinedDistanceValue < 0.0) {
                    return -1.0;
                }
                if (SALSAUtility.DistanceProcessingOption == 2) {
                    return SALSAUtility.UndefinedDistanceValue * SALSAUtility.UndefinedDistanceValue;
                } else {
                    return SALSAUtility.UndefinedDistanceValue;
                }
            } else {
                return 1.0;
            }
        }


        if (Temp < 0) {
            if (SALSAUtility.UndefinedDistanceValue < 0.0) {
                return -1.0;
            }
            if (SALSAUtility.DistanceProcessingOption == 2) {
                return SALSAUtility.UndefinedDistanceValue * SALSAUtility.UndefinedDistanceValue;
            } else {
                return SALSAUtility.UndefinedDistanceValue;
            }
        }

        double distanceValue = (Temp / (Short.MAX_VALUE * 1.0));

        if (SALSAUtility.DistanceProcessingOption == 2) {
            distanceValue = distanceValue * distanceValue;
        }
        return distanceValue;
    }

    /**
     * Puts a given distance value to the corresponding row and col in Used points
     *
     * @param row   row number from Used points
     * @param col   col number from Used points
     * @param value the distance value to be put
     */
    public static void putDistanceValue(int row, int col, double value) {

        if (SALSAUtility.StoredDistanceOption == 2) {
            row = row - SALSAUtility.PointStart_Process;
        }
        if (SALSAUtility.StoredDistanceOption == 3) {
            int originalpoint = SALSAUtility.UsedPointtoOriginalPointMap[row];
            int variedpoint = SALSAUtility.OriginalPointDisposition[originalpoint] - SALSAUtility.SALSASHIFT;
            if (variedpoint < 0) {
                SALSAUtility.printAndThrowRuntimeException(
                        " Illegal Distance Put Request Used Point " + row + " Original " + originalpoint);

            }
            row = variedpoint - SALSAUtility.VariedPointStart_Process;
        }

        Short Temp;

        if (value > 1.0) {
            SALSAUtility.printAndThrowRuntimeException(" Illegal Distance value Put Request Used Point " + String
                    .format("%.4f", value) + " Coordinates " + row + " " + col);

        }
        if ((value == 1.0) && (SALSAUtility.DistanceCut > 0.0) && (SALSAUtility.UndefinedDistanceValue < 0.0)) { //
            // Inconsistent value
            SALSAUtility.printAndThrowRuntimeException(" Illegal Distance value 1.0 Put Request Used Point " + String
                    .format("%.4f", value) + " Coordinates " + row + " " + col);

        }
        if (value < 0.0) {
            if (SALSAUtility.DistanceCut < 0.0) {
                SALSAUtility.printAndThrowRuntimeException(
                        " Illegal Distance value < 0.0 Put Request Used Point " + String
                                .format("%.4f", value) + " Coordinates " + row + " " + col);

            }

            Temp = Short.MAX_VALUE;
        } else {
            Temp = (short) (value * Short.MAX_VALUE);
        }
        SALSAUtility.PointDistances[row][col] = Temp;
    }
}