package salsa.mdsaschisq;

import com.google.common.base.Strings;
import edu.rice.hj.api.SuspendableException;
import mpi.MPIException;
import tangible.RefObject;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Hashtable;
import java.util.regex.Pattern;

import static edu.rice.hj.Module1.forallChunked;

public class ManxcatMDS {
    //  Input Variables
    public static int Chisqnorm = 0; // Option for Chisq Norm
    public static int DistanceFormula = 1; // =2 Square of Euclidean Distance, = 1 EuclideanDistance
    public static int WeightingOption = 0; // Weight Option
    public static int NumLinksHistogramBins = 100; // Bins in Links Histogram
    public static int LinkCutforCenter = 20; // Ensure Center has at least this many links

    //  Calculated Variables
    public static int[] PointStatus; // Set status of points = -1 Deleted and Fixed = 0 normal, = 1 Center, =2 x axis, =3 in x,y plane NOT DECOMPOSED
    public static double[] PointWeights; // Intrinsic weight of each point
    public static double SystemRadius = 0.0; //  Radius of System
    public static double SystemMax = 0.0; // Maximumum value of distance
    public static double SystemAverage = 0.0; //  Average Distance in System
    public static double SystemSigma = 0.0; // Standard Deviation in System
    public static double SQRTSystemAverage = 0.0; //  SQRT Average Distance in System
    public static double MinimumDistance; // Minimum Distance for weight calculation
    public static double SQRTMinimumDistance; // Minimum Distance for weight calculation
    public static double SQUAREMinimumDistance; // Minimum Distance for weight calculation
    public static boolean FindPointstoFix = true; // If true find key points to fix

    public static void SetupHotsunforMDS() {

        int NumberofPoints = SALSAUtility.PointCount_Global;
        Hotsun.ParameterVectorDimension = ManxcatCentral.config.LocalVectorDimension;
        Hotsun.ndata = NumberofPoints * (NumberofPoints - 1L) / 2L;
        Hotsun.Number_VectorParameters = NumberofPoints;
        Hotsun.npar = Hotsun.Number_VectorParameters * Hotsun.ParameterVectorDimension;
        SALSAUtility.sequentialBLAS = false;
        Hotsun.DecomposeParameters = true;
        Hotsun.fullmatrixset = false;
    }

    //  Set Chisq version of MDS
    public static void SetupMDSasChisq() throws MPIException {
        int NumberofPoints = SALSAUtility.PointCount_Global;
        Hotsun.ndata = NumberofPoints * (NumberofPoints - 1L) / 2L;

        // Initialize Linear Algebra
        MDSLinearAlgebra.Initialize();

        // Read in distance data
        SALSAUtility.DistanceProcessingOption = ManxcatCentral.config.DistanceProcessingOption;
        try {
            SALSAParallelism.ReadDataFromFile(ManxcatCentral.config.DistanceMatrixFile);
        } catch (IOException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        } catch (MPIException e) {
            e.printStackTrace();
        }

        // If FindPointstoFix true, all fixed parameters are set to 0 UNLESS initialized differently

        FindPointstoFix = SALSAUtility.NumberFixedPoints <= 0;

        if (ManxcatCentral.config.InitializationOption > 0) {
            ManxcatCentral.config.InitializationLoops = 1; // Can't loop if fixed initial positions
        }


        // set up Fixed and Deleted Parameters
        PointStatus = new int[SALSAUtility.PointCount_Global];
        for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++) {
            PointStatus[GlobalPointIndex] = 0;
        }

        SALSAUtility.SetupDistanceWeights();

        if (SALSAUtility.TransformMethod == 10 || (SALSAUtility.NumberDistanceWeightCuts > 0)) {
            for (int cutloop = 0; cutloop < SALSAUtility.NumberDistanceWeightCuts; cutloop++) {
                SALSAUtility.DistanceWeightCuts[cutloop] = Math.pow(SALSAUtility.DistanceWeightCuts[cutloop],
                        SALSAUtility.TransformParameter);
            }
        } else if (SALSAUtility.TransformMethod == 8 || SALSAUtility.TransformMethod == 9) {
            SALSAUtility.NumberDistanceWeightCuts = 0;
        }


        // Loop over initial analysis until no more deleted points
        SALSAUtility.SALSAPrint(1, "\nInitial Processing Parameters\nDistance Cut " + String
                .format("%.3f",
                        SALSAUtility.DistanceCut) + " Link Cut " + SALSAUtility.LinkCut + " Allowed Deleted Fraction " + String
                .format("%.3f", SALSAUtility.AllowedDeletedFraction) + " Undefined Distance Value " + String
                .format("%.3f",
                        SALSAUtility.UndefinedDistanceValue) + "\nDistance Transformation Method " + SALSAUtility.TransformMethod + " with Parameter " + String


                .format("%.4f", SALSAUtility.TransformParameter));
        double MissingDistances = 0.0;
        int DisconnectedPoints = 0;
        int DisconnectedLoopCount = 1;
        while (true) {
            int oldDisconnectedPoints = DisconnectedPoints;
            SALSAUtility.SALSAPrint(1,
                    "\n ******* Loop over Identification of Disconnected Points " + DisconnectedLoopCount);
            RefObject<Double> tempRef_SystemAverage = new RefObject<Double>(SystemAverage);
            RefObject<Double> tempRef_SystemMax = new RefObject<Double>(SystemMax);
            RefObject<Double> tempRef_SystemSigma = new RefObject<Double>(SystemSigma);
            RefObject<Integer> tempRef_DisconnectedPoints = new RefObject<Integer>(DisconnectedPoints);
            RefObject<Double> tempRef_MissingDistances = new RefObject<Double>(MissingDistances);
            ManxcatMDSBasicDataProcessing.IntialDistanceAnalysis(tempRef_SystemAverage, tempRef_SystemMax,
                    tempRef_SystemSigma, tempRef_DisconnectedPoints, tempRef_MissingDistances);
            SystemAverage = tempRef_SystemAverage.argValue;
            SystemMax = tempRef_SystemMax.argValue;
            SystemSigma = tempRef_SystemSigma.argValue;
            DisconnectedPoints = tempRef_DisconnectedPoints.argValue;
            MissingDistances = tempRef_MissingDistances.argValue;
            SQRTSystemAverage = Math.sqrt(SystemAverage);
            if (oldDisconnectedPoints >= DisconnectedPoints) {
                break;
            }
            int AllowedDisconnectedNumber = (int) SALSAUtility.AllowedDeletedFraction + SALSAUtility.PointCount_Global;
            if (DisconnectedPoints > AllowedDisconnectedNumber) {
                SALSAUtility.printAndThrowRuntimeException(
                        "Must stop as Number of Disconnected Points " + DisconnectedPoints + " Exceeds Cut " +
                                AllowedDisconnectedNumber
                );

            }
            ++DisconnectedLoopCount;
        }
        int deletedpoints = 0;
        for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++) {
            if (PointStatus[GlobalPointIndex] == -1) {
                for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++) {
                    Hotsun.FixedParameter[GlobalPointIndex][LocalVectorIndex] = true;
                }
                ++deletedpoints;
            }
        }
        if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0)) {
            SALSAUtility.SALSAPrint(1, "Deleted Points due to bad link count " + deletedpoints);
        }

        // If necessary Clean up and Transform Distances overriding previous averages
        RefObject<Double> tempRef_SystemAverage2 = new RefObject<Double>(SystemAverage);
        RefObject<Double> tempRef_SystemMax2 = new RefObject<Double>(SystemMax);
        RefObject<Double> tempRef_SystemSigma2 = new RefObject<Double>(SystemSigma);
        RefObject<Integer> tempRef_DisconnectedPoints2 = new RefObject<Integer>(DisconnectedPoints);
        RefObject<Double> tempRef_MissingDistances2 = new RefObject<Double>(MissingDistances);
        ManxcatMDSBasicDataProcessing.CleanandTransformDistances(tempRef_SystemAverage2, tempRef_SystemMax2,
                tempRef_SystemSigma2, tempRef_DisconnectedPoints2, tempRef_MissingDistances2);
        SystemAverage = tempRef_SystemAverage2.argValue;
        SystemMax = tempRef_SystemMax2.argValue;
        SystemSigma = tempRef_SystemSigma2.argValue;
        DisconnectedPoints = tempRef_DisconnectedPoints2.argValue;
        MissingDistances = tempRef_MissingDistances2.argValue;

        // If ManxcatCentral.MetadataforRun.MinimumDistance  is positive, it is absolute Minimum Distance
        // If ManxcatCentral.MetadataforRun.MinimumDistance  is negative, it is - multiplier of System Average
        if (ManxcatCentral.config.MinimumDistance < 0) {
            MinimumDistance = -ManxcatCentral.config.MinimumDistance * SystemAverage;
        } else {
            MinimumDistance = ManxcatCentral.config.MinimumDistance;
        }
        SQRTMinimumDistance = Math.sqrt(MinimumDistance);
        SQUAREMinimumDistance = MinimumDistance * MinimumDistance;

        //  Now find centers and histogram distances
        int PointsinDistanceHistogram = ManxcatCentral.config.HistogramBinCount;
        double Histmin = 0.0;
        double Histmax = SystemMax;
        RefObject<Double> tempRef_Histmin = new RefObject<Double>(Histmin);
        RefObject<Double> tempRef_Histmax = new RefObject<Double>(Histmax);
        ManxcatMDSBasicDataProcessing.SetUpHistogramRange(PointsinDistanceHistogram, tempRef_Histmin, tempRef_Histmax);
        Histmin = tempRef_Histmin.argValue;
        Histmax = tempRef_Histmax.argValue;
        double[] Bincounts = new double[2 + PointsinDistanceHistogram];
        int Center = 0;
        RefObject<Integer> tempRef_Center = new RefObject<Integer>(Center);
        RefObject<Double> tempRef_SystemRadius = new RefObject<Double>(SystemRadius);
        RefObject<double[]> tempRef_Bincounts = new RefObject<double[]>(Bincounts);
        ManxcatMDSBasicDataProcessing.FindCenter(tempRef_Center, tempRef_SystemRadius, Histmin, Histmax,
                PointsinDistanceHistogram, tempRef_Bincounts);
        Center = tempRef_Center.argValue;
        SystemRadius = tempRef_SystemRadius.argValue;
        Bincounts = tempRef_Bincounts.argValue;

        PointStatus[Center] = 1;
        if (FindPointstoFix) {
            for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++) {
                Hotsun.FixedParameter[Center][LocalVectorIndex] = true;
            }
        }

        //  Set up Weights
        PointWeights = new double[SALSAUtility.PointCount_Global];
        WeightingOption = ManxcatCentral.config.WeightingOption;
        SetupWeightings(PointWeights);

        double DistancesNearEachOther = 0.0;
        int NotLonelyPoints = 0;
        int xAxis = 0;
        double xAxisExtent = 0.0;
        RefObject<Integer> tempRef_xAxis = new RefObject<Integer>(xAxis);
        RefObject<Double> tempRef_xAxisExtent = new RefObject<Double>(xAxisExtent);
        RefObject<Double> tempRef_DistancesNearEachOther = new RefObject<Double>(DistancesNearEachOther);
        RefObject<Integer> tempRef_NotLonelyPoints = new RefObject<Integer>(NotLonelyPoints);
        ManxcatMDSBasicDataProcessing.FindxAxis(Center, tempRef_xAxis, tempRef_xAxisExtent,
                tempRef_DistancesNearEachOther, tempRef_NotLonelyPoints);
        xAxis = tempRef_xAxis.argValue;
        xAxisExtent = tempRef_xAxisExtent.argValue;
        DistancesNearEachOther = tempRef_DistancesNearEachOther.argValue;
        NotLonelyPoints = tempRef_NotLonelyPoints.argValue;
        PointStatus[xAxis] = 2;


        if (FindPointstoFix) {
            for (int LocalVectorIndex = 1; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++) {
                Hotsun.FixedParameter[xAxis][LocalVectorIndex] = true;
            }
        }

        int xyPlane = 0;
        double xyPlaneExtent = 0.0;
        RefObject<Integer> tempRef_xyPlane = new RefObject<Integer>(xyPlane);
        RefObject<Double> tempRef_xyPlaneExtent = new RefObject<Double>(xyPlaneExtent);
        ManxcatMDSBasicDataProcessing.FindxyPlane(Center, xAxis, tempRef_xyPlane, tempRef_xyPlaneExtent);
        xyPlane = tempRef_xyPlane.argValue;
        xyPlaneExtent = tempRef_xyPlaneExtent.argValue;
        PointStatus[xyPlane] = 3;

        if (FindPointstoFix) {
            for (int LocalVectorIndex = 2; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++) {
                Hotsun.FixedParameter[xyPlane][LocalVectorIndex] = true;
            }
        }

        if (!FindPointstoFix) {
            if (SALSAUtility.NumberFixedPoints < (Hotsun.ParameterVectorDimension - 1)) {
                SALSAUtility.printAndThrowRuntimeException("Not enough fixed points " + SALSAUtility.NumberFixedPoints);

            }
        }

        if (SALSAUtility.NumberFixedPoints > 0) {
            for (int FixedIndex = 0; FixedIndex < SALSAUtility.NumberFixedPoints; FixedIndex++) {
                int ParameterIndex = SALSAUtility.FixedPointOriginal[FixedIndex];
                ParameterIndex = SALSAUtility.OriginalPointtoUsedPointMap[ParameterIndex];
                // This is Used point number
                if (ParameterIndex < 0) {
                    SALSAUtility.printAndThrowRuntimeException("Illegal Fixed Point " + FixedIndex);

                }
                for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++) {
                    Hotsun.FixedParameter[ParameterIndex][LocalVectorIndex] = true;
                }
            }
        }

        // Normalize Chisq for MDS
        ManxcatCentral.ChisqPrintConstant = 0.5 * ManxcatCentral.ChisqPrintConstant / (Hotsun.ndata - 0.5 * MissingDistances - Hotsun.npar); // Print Chisq per point (factor of .5 as each point doubled)
        Chisqnorm = ManxcatCentral.config.Chisqnorm;
        DistanceFormula = ManxcatCentral.config.DistanceFormula;

        if (DistanceFormula == 0) {
            DistanceFormula = 2;
        }

        if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0)) {
            double tmp = SALSAUtility.NumberVariedPoints;
            double fraction = MissingDistances / (tmp * (tmp - 1.0));
            double EstimatedDimension = 2.0 * SystemAverage * SystemAverage / (SystemSigma * SystemSigma);
            SALSAUtility.SALSAPrint(1,
                    "\nAFTER CLEAN Disconnected Points " + DisconnectedPoints + " Missing Distances " + String

                            .format("%.0f", MissingDistances) + " Fraction " + String
                            .format("%.4f", fraction)
            );
            SALSAUtility.SALSAPrint(1, "AFTER TRANSFORM Max " + String.format("%.4E", SystemMax) + " Average " + String
                    .format("%.4E", SystemAverage) + " Sigma " + String
                    .format("%.4E", SystemSigma) + " Estimated Dimension " + String
                    .format("%.2f",
                            EstimatedDimension) + "\n Center " + Center + " With Link Cut " + LinkCutforCenter + " " +
                    "Radius " + String
                    .format("%.4E", SystemRadius) + " xAxis " + xAxis + " " + String
                    .format("%.4E", xAxisExtent) + " xyPlane " + xyPlane + " " + String
                    .format("%.4E", xyPlaneExtent) + "\n Minimum Distance " + String
                    .format("%.4E", MinimumDistance) + " Distances Less than this " + String
                    .format("%.0f", DistancesNearEachOther) + " Points Affected " + NotLonelyPoints);
            String histogramcounts = "";
            for (int binloop = 0; binloop < (2 + PointsinDistanceHistogram); binloop++) {
                histogramcounts += String.format("%.0f", Bincounts[binloop]) + ", ";
            }
            double BinSize = (Histmax - Histmin) / PointsinDistanceHistogram;
            SALSAUtility.SALSAPrint(1,
                    "\nDistance Histogram Min " + String.format("%.4E", Histmin) + " Max " + String.format("%.4E",
                            Histmax) + " Binsize " + String.format("%.4E",
                            BinSize) + " #Counts with under/overflow\n" + histogramcounts);
        }
        return;

    } //  End SetupChisq

    public static void Sequel() throws MPIException {
        /* Sequel of ManxcatMDS work */

        SALSAUtility.SALSAPrint(1, "\nStarting density and histogram generation");

        java.util.Hashtable orignalPnumToCnumTable = new java.util.Hashtable();
        if (SALSAUtility.IsClustersSelected) {
            ReadClusterFile(orignalPnumToCnumTable);
        }

        double xminWhole = Double.MAX_VALUE;
        double xmaxWhole = -Double.MAX_VALUE;
        double yminWhole = Double.MAX_VALUE;
        double ymaxWhole = -Double.MAX_VALUE;

        double xminSelected = Double.MAX_VALUE;
        double xmaxSelected = -Double.MAX_VALUE;
        double yminSelected = Double.MAX_VALUE;
        double ymaxSelected = -Double.MAX_VALUE;

        double xminSelectedInter = Double.MAX_VALUE;
        double xmaxSelectedInter = -Double.MAX_VALUE;
        double yminSelectedInter = Double.MAX_VALUE;
        double ymaxSelectedInter = -Double.MAX_VALUE;

        SALSAUtility.SALSAPrint(1, "\n\tFinding min/max values");

        RefObject<Double> tempRef_xminWhole = new RefObject<Double>(xminWhole);
        RefObject<Double> tempRef_xmaxWhole = new RefObject<Double>(xmaxWhole);
        RefObject<Double> tempRef_yminWhole = new RefObject<Double>(yminWhole);
        RefObject<Double> tempRef_ymaxWhole = new RefObject<Double>(ymaxWhole);
        RefObject<Double> tempRef_xminSelected = new RefObject<Double>(xminSelected);
        RefObject<Double> tempRef_xmaxSelected = new RefObject<Double>(xmaxSelected);
        RefObject<Double> tempRef_yminSelected = new RefObject<Double>(yminSelected);
        RefObject<Double> tempRef_ymaxSelected = new RefObject<Double>(ymaxSelected);
        RefObject<Double> tempRef_xminSelectedInter = new RefObject<Double>(xminSelectedInter);
        RefObject<Double> tempRef_xmaxSelectedInter = new RefObject<Double>(xmaxSelectedInter);
        RefObject<Double> tempRef_yminSelectedInter = new RefObject<Double>(yminSelectedInter);
        RefObject<Double> tempRef_ymaxSelectedInter = new RefObject<Double>(ymaxSelectedInter);
        ManxcatMDSBasicDataProcessing.FindMinMaxForDensity(tempRef_xminWhole, tempRef_xmaxWhole, tempRef_yminWhole,
                tempRef_ymaxWhole, tempRef_xminSelected, tempRef_xmaxSelected, tempRef_yminSelected,
                tempRef_ymaxSelected, tempRef_xminSelectedInter, tempRef_xmaxSelectedInter, tempRef_yminSelectedInter,
                tempRef_ymaxSelectedInter, orignalPnumToCnumTable);
        xminWhole = tempRef_xminWhole.argValue;
        xmaxWhole = tempRef_xmaxWhole.argValue;
        yminWhole = tempRef_yminWhole.argValue;
        ymaxWhole = tempRef_ymaxWhole.argValue;
        xminSelected = tempRef_xminSelected.argValue;
        xmaxSelected = tempRef_xmaxSelected.argValue;
        yminSelected = tempRef_yminSelected.argValue;
        ymaxSelected = tempRef_ymaxSelected.argValue;
        xminSelectedInter = tempRef_xminSelectedInter.argValue;
        xmaxSelectedInter = tempRef_xmaxSelectedInter.argValue;
        yminSelectedInter = tempRef_yminSelectedInter.argValue;
        ymaxSelectedInter = tempRef_ymaxSelectedInter.argValue;
        SALSAUtility.SALSAPrint(1, String.format("\n\t\txminWhole: %1$s\txmaxWhole: %2$s", xminWhole, xmaxWhole));
        SALSAUtility.SALSAPrint(1, String.format("\t\tyminWhole: %1$s\tymaxWhole: %2$s", yminWhole, ymaxWhole));
        SALSAUtility.SALSAPrint(1,
                String.format("\n\t\txminSelected: %1$s\txmaxSelected: %2$s", xminSelected, xmaxSelected));
        SALSAUtility.SALSAPrint(1,
                String.format("\t\tyminSelected: %1$s\tymaxSelected: %2$s", yminSelected, ymaxSelected));
        SALSAUtility.SALSAPrint(1,
                String.format("\n\t\txminSelectedInter: %1$s\txmaxSelectedInter: %2$s", xminSelectedInter,
                        xmaxSelectedInter));
        SALSAUtility.SALSAPrint(1,
                String.format("\t\tyminSelectedInter: %1$s\tymaxSelectedInter: %2$s", yminSelectedInter,
                        ymaxSelectedInter));

        SALSAUtility.SALSAPrint(1, "\n\tDone");


        double deltaxWhole = (xmaxWhole - xminWhole) / SALSAUtility.Xres;
        double deltayWhole = (ymaxWhole - yminWhole) / SALSAUtility.Yres;
        double deltasWhole = deltaxWhole * deltayWhole;

        double deltaxSelected = (xmaxSelected - xminSelected) / SALSAUtility.Xres;
        double deltaySelected = (ymaxSelected - yminSelected) / SALSAUtility.Yres;
        double deltasSelected = deltaxSelected * deltaySelected;

        double deltaxSelectedInter = (xmaxSelectedInter - xminSelectedInter) / SALSAUtility.Xres;
        double deltaySelectedInter = (ymaxSelectedInter - yminSelectedInter) / SALSAUtility.Yres;
        double deltasSelectedInter = deltaxSelectedInter * deltaySelectedInter;

        double[][] densityMatrixWhole = null, densityMatrixSelected = null, densityMatrixSelectedInter = null;
        double[] xHistogramWhole = null, yHistogramWhole = null;
        double[] xHistogramSelected = null, yHistogramSelected = null;
        double[] xHistogramSelectedInter = null, yHistogramSelectedInter = null;

        double countWhole = 0;
        double countSelected = 0;
        double countSelectedInter = 0;

        SALSAUtility.SALSAPrint(1, "\n\tGenerating density matrices and histograms");
        RefObject<double[][]> tempRef_densityMatrixWhole = new RefObject<double[][]>(densityMatrixWhole);
        RefObject<double[]> tempRef_xHistogramWhole = new RefObject<double[]>(xHistogramWhole);
        RefObject<double[]> tempRef_yHistogramWhole = new RefObject<double[]>(yHistogramWhole);
        RefObject<double[][]> tempRef_densityMatrixSelected = new RefObject<double[][]>(densityMatrixSelected);
        RefObject<double[]> tempRef_xHistogramSelected = new RefObject<double[]>(xHistogramSelected);
        RefObject<double[]> tempRef_yHistogramSelected = new RefObject<double[]>(yHistogramSelected);
        RefObject<double[][]> tempRef_densityMatrixSelectedInter = new RefObject<double[][]>(
                densityMatrixSelectedInter);
        RefObject<double[]> tempRef_xHistogramSelectedInter = new RefObject<double[]>(xHistogramSelectedInter);
        RefObject<double[]> tempRef_yHistogramSelectedInter = new RefObject<double[]>(yHistogramSelectedInter);
        RefObject<Double> tempRef_countWhole = new RefObject<Double>(countWhole);
        RefObject<Double> tempRef_countSelected = new RefObject<Double>(countSelected);
        RefObject<Double> tempRef_countSelectedInter = new RefObject<Double>(countSelectedInter);
        ManxcatMDSBasicDataProcessing.GenerateDensityMatrix(tempRef_densityMatrixWhole, tempRef_xHistogramWhole,
                tempRef_yHistogramWhole, tempRef_densityMatrixSelected, tempRef_xHistogramSelected,
                tempRef_yHistogramSelected, tempRef_densityMatrixSelectedInter, tempRef_xHistogramSelectedInter,
                tempRef_yHistogramSelectedInter, xminWhole, xmaxWhole, yminWhole, ymaxWhole, xminSelected, xmaxSelected,
                yminSelected, ymaxSelected, xminSelectedInter, xmaxSelectedInter, yminSelectedInter, ymaxSelectedInter,
                tempRef_countWhole, tempRef_countSelected, tempRef_countSelectedInter, orignalPnumToCnumTable);
        densityMatrixWhole = tempRef_densityMatrixWhole.argValue;
        xHistogramWhole = tempRef_xHistogramWhole.argValue;
        yHistogramWhole = tempRef_yHistogramWhole.argValue;
        densityMatrixSelected = tempRef_densityMatrixSelected.argValue;
        xHistogramSelected = tempRef_xHistogramSelected.argValue;
        yHistogramSelected = tempRef_yHistogramSelected.argValue;
        densityMatrixSelectedInter = tempRef_densityMatrixSelectedInter.argValue;
        xHistogramSelectedInter = tempRef_xHistogramSelectedInter.argValue;
        yHistogramSelectedInter = tempRef_yHistogramSelectedInter.argValue;
        countWhole = tempRef_countWhole.argValue;
        countSelected = tempRef_countSelected.argValue;
        countSelectedInter = tempRef_countSelectedInter.argValue;

        SALSAUtility.SALSAPrint(1, "\tDone");
        if (SALSAUtility.MPI_Rank == 0) {
			/* Density matrices and histograms should be good here */
            String outDir = (new java.io.File(ManxcatCentral.config.ReducedVectorOutputFileName)).getParent();
            SALSAUtility.SALSAPrint(1, "\n\tGenerating files for whole sample");
            GenerateDensityDataFile(densityMatrixWhole, xHistogramWhole, yHistogramWhole, xminWhole, xmaxWhole,
                    yminWhole, ymaxWhole, deltaxWhole, deltayWhole, deltasWhole, countWhole, outDir, "whole");
            SALSAUtility.SALSAPrint(1, "\tDone");
            if (SALSAUtility.IsClustersSelected) {
                SALSAUtility.SALSAPrint(1, "\n\tGenerating files for selected cluster (intra)");
                GenerateDensityDataFile(densityMatrixSelected, xHistogramSelected, yHistogramSelected, xminSelected,
                        xmaxSelected, yminSelected, ymaxSelected, deltaxSelected, deltaySelected, deltasSelected,
                        countSelected, outDir, "selected");
                SALSAUtility.SALSAPrint(1, "\tDone");

                SALSAUtility.SALSAPrint(1, "\n\tGenerating files for selected cluster (inter)");
                GenerateDensityDataFile(densityMatrixSelectedInter, xHistogramSelectedInter, yHistogramSelectedInter,
                        xminSelectedInter, xmaxSelectedInter, yminSelectedInter, ymaxSelectedInter, deltaxSelectedInter,
                        deltaySelectedInter, deltasSelectedInter, countSelectedInter, outDir, "selected-inter");
                SALSAUtility.SALSAPrint(1, "\tDone");
            }
        }
    }

    private static void GenerateDensityDataFile(double[][] densityMatrx, double[] xHist, double[] yHist, double xmin,
                                                double xmax, double ymin, double ymax, double deltax, double deltay,
                                                double deltas, double count, String outDir, String prefix) {
        double cellmax = 0.0;
        for (int i = 0; i < SALSAUtility.Yres; i++) {
            for (int j = 0; j < SALSAUtility.Xres; j++) {
                if (densityMatrx[i][j] > cellmax) {
                    cellmax = densityMatrx[i][j];
                }
            }
        }
        double cellmean = count / (SALSAUtility.Xres * SALSAUtility.Yres);
        double power = cellmax < (SALSAUtility.Alpha * cellmean) ? 1.0 : (Math.log(SALSAUtility.Alpha) / Math.log(
                cellmax / cellmean));
        // Constant value by which the number of points in a 2D square is multiplied.
        // The resulting value is independent of the total number of points as well as
        // the x,y resolution. The mult value is a factor changing the z value scale.
        double c = 1.0 / cellmax;

        // Output density values
        System.out.println("****************************************");
        System.out.println("DataSet\t" + prefix);
        System.out.println("Count\t" + count);
        System.out.println("Deltas\t" + deltas);
        System.out.println("CellMean\t" + cellmean);
        System.out.println("CellMax\t" + cellmax);
        System.out.println("Power\t" + power);
        System.out.println("Const\t" + c);
        for (int i = 0; i < 10; i++) {
            double density = i / 10.0;
            double densityToCount = Math.pow(density, (1 / power)) / c;
            System.out.println(density + "\t" + densityToCount);
        }
        System.out.println("****************************************");


        int xpointcount = 2 * SALSAUtility.Xres;
        int ypointcount = 2 * SALSAUtility.Yres;

        Path densityFile = Paths.get(outDir, prefix + "-density.txt");

        Path xHistFile = Paths.get(outDir, prefix + "-xHist.txt");

        Path yHistFile = Paths.get(outDir, prefix + "-yHist.txt");

        Path scriptFile = Paths.get(outDir, prefix + "-plot.txt");


        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(densityFile, Charset.defaultCharset()));
             PrintWriter xHistWriter = new PrintWriter(Files.newBufferedWriter(xHistFile, Charset.defaultCharset()));
             PrintWriter yHistWriter = new PrintWriter(Files.newBufferedWriter(yHistFile, Charset.defaultCharset()));
             PrintWriter scriptWriter = new PrintWriter(Files.newBufferedWriter(scriptFile, Charset.defaultCharset()))) {

            writer.println("#xcoord\tycoord\thistogramValue");
            xHistWriter.println("#xval\thistogramvalue");
            yHistWriter.println("#yval\thistogramvalue");

            // Generating x histogram
            double xoffset = xmin + 0.5 * deltax;
            for (int i = 0; i < SALSAUtility.Xres; ++i) {
                double xcoord = xoffset + i * deltax;
                if (SALSAUtility.Normalize) {
                    xcoord /= xmax;
                }
                xHistWriter.println(xcoord + "\t" + xHist[i]);
            }

            // Generating y histogram
            double yoffset = ymin + 0.5 * deltay;
            for (int i = 0; i < SALSAUtility.Yres; ++i) {
                double ycoord = yoffset + i * deltay;
                if (SALSAUtility.Normalize) {
                    ycoord /= ymax;
                }
                yHistWriter.println(ycoord + "\t" + yHist[i]);
            }

            for (int i = 0; i < xpointcount; i++) {
                double x = xmin + ((IsOdd(i) ? (i + 1) / 2 : i / 2) * deltax);
                int cellx = IsOdd(i) ? (i - 1) / 2 : i / 2;

                for (int j = 0; j < ypointcount; j++) {
                    double y = ymin + ((IsOdd(j) ? (j + 1) / 2 : j / 2) * deltay);
                    int celly = IsOdd(j) ? (j - 1) / 2 : j / 2;

                    double cellvalue = Math.pow((densityMatrx[celly][cellx] * c), power);

                    cellvalue = cellvalue > SALSAUtility.Pcutf ? SALSAUtility.Pcutf : cellvalue;

                    if (SALSAUtility.Normalize) {
                        writer.println(x / xmax + "\t" + y / ymax + "\t" + cellvalue);
                    } else {
                        writer.println(x + "\t" + y + "\t" + cellvalue);
                    }
                }
                writer.println();
            }

            // Fill up the remaining region from beyond x=xmax and y=ymax as zero
            writer.println();
            if (SALSAUtility.Normalize) {
                writer.println(xmin / xmax + "\t" + ymax / ymax + "\t" + 0.0);
                writer.println(xmin / xmax + "\t" + SALSAUtility.Xmaxbound + "\t" + 0.0);
                writer.println();
                writer.println(xmax / xmax + "\t" + ymax / ymax + "\t" + 0.0);
                writer.println(xmax / xmax + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                writer.println();
                writer.println(xmax / xmax + "\t" + ymin / ymax + "\t" + 0.0);
                writer.println(xmax / xmax + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                writer.println();
                writer.println(SALSAUtility.Xmaxbound + "\t" + ymin / ymax + "\t" + 0.0);
                writer.println(SALSAUtility.Xmaxbound + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
            } else {
                writer.println(xmin + "\t" + ymax + "\t" + 0.0);
                writer.println(xmin + "\t" + SALSAUtility.Xmaxbound + "\t" + 0.0);
                writer.println();
                writer.println(xmax + "\t" + ymax + "\t" + 0.0);
                writer.println(xmax + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                writer.println();
                writer.println(xmax + "\t" + ymin + "\t" + 0.0);
                writer.println(xmax + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                writer.println();
                writer.println(SALSAUtility.Xmaxbound + "\t" + ymin + "\t" + 0.0);
                writer.println(SALSAUtility.Xmaxbound + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
            }


            WriteGnuplotScript(prefix, densityFile, xHistFile, yHistFile, scriptWriter);
        } catch (IOException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }

    private static void WriteGnuplotScript(String prefix, Path densityFile, Path xHistFile, Path yHistFile,
                                           PrintWriter scriptWriter) {

        scriptWriter.println("set terminal png truecolor nocrop font arial 14 size 1200,1200 xffffff");
        scriptWriter.println();

        String pngfile = prefix + (SALSAUtility.Normalize ? "-normalized" : "") + "-plot.png";
        scriptWriter.println("set output '" + pngfile + "'");

        scriptWriter.println("set size 1.0, 1.0");
        scriptWriter.println("set multiplot");

        scriptWriter.println();

        // Title box
        scriptWriter.println("set origin 0.0, 0.85");
        scriptWriter.println("set size 0.95, 0.1");
        scriptWriter.println("set border linecolor rgbcolor \"white\"");
        scriptWriter.println("unset key");
        String title = SALSAUtility.ManxcatRunName + " (" + prefix + (SALSAUtility.Normalize ? "-normalized" : "") + ")";
        scriptWriter.println("set title \"" + title + "\" textcolor rgbcolor \"black\"");
        scriptWriter.println("plot [0:1] [0:1] 0.0 lt rgb \"white\"");

        scriptWriter.println("set border linecolor rgbcolor \"black\"");
        scriptWriter.println("set dummy u,v");
        scriptWriter.println("unset key");
        scriptWriter.println("set size ratio 1.0");
        scriptWriter.println("set style fill  solid 0.85 noborder");
        scriptWriter.println("set style line 1 lt 1 lw 4");
        scriptWriter.println("set pm3d map");
        scriptWriter.println("set palette rgbformulae 30,31,32 model RGB negative");

        scriptWriter.println();

        // Y histogram (rotated)
        scriptWriter.println("set origin 0.0, 0.45");
        scriptWriter.println("set size 0.45, 0.45");
        scriptWriter.println("set xtics rotate by -90");
        String xlabel = "Count";
        scriptWriter.println("set xlabel \"" + xlabel + "\" textcolor rgbcolor \"black\"");
        String ylabel = "Euclidean Distance";
        scriptWriter.println("set ylabel \"" + ylabel + "\" textcolor rgbcolor \"black\"");
        title = "Histogram (rotated) of " + ylabel;
        scriptWriter.println("set title \"" + title + "\" textcolor rgbcolor \"black\"");
        scriptWriter.println(
                "plot [][:" + SALSAUtility.Ymaxbound + "] '" + yHistFile.getFileName() + "' using 2:1 with filledcurves y1 lt rgb \"black\"");

        scriptWriter.println("set xtics rotate by 0");
        scriptWriter.println();

        // Density plot
        scriptWriter.println("set origin 0.45, 0.45");
        scriptWriter.println("set size 0.5, 0.5");
        xlabel = "Original Distance";
        scriptWriter.println("set xlabel \"" + xlabel + "\" textcolor rgbcolor \"black\"");
        ylabel = "Euclidena Distance";
        scriptWriter.println("set ylabel \"" + ylabel + "\" textcolor rgbcolor \"black\"");
        title = "Heat Map of " + ylabel + " vs " + xlabel;
        scriptWriter.println("set title \"" + title + "\" textcolor rgbcolor \"black\"");
        scriptWriter.println(
                "splot [:" + SALSAUtility.Xmaxbound + "] [:" + SALSAUtility.Ymaxbound + "] '" + densityFile.getFileName() + "'");


        scriptWriter.println();

        // Y histogram (unrotated)
        scriptWriter.println("set origin 0.0, 0.0");
        scriptWriter.println("set size 0.45, 0.45");
        xlabel = "Euclidean Distance";
        scriptWriter.println("set xlabel \"" + xlabel + "\" textcolor rgbcolor \"black\"");
        ylabel = "Count";
        scriptWriter.println("set ylabel \"" + ylabel + "\" textcolor rgbcolor \"black\"");
        title = "Histogram of " + xlabel;
        scriptWriter.println("set title \"" + title + "\" textcolor rgbcolor \"black\"");
        scriptWriter.println(
                "plot [:" + SALSAUtility.Ymaxbound + "] []'" + yHistFile.getFileName() + "' with filledcurves x1 lt rgb \"black\"");


        scriptWriter.println();

        // X histogram
        scriptWriter.println("set origin 0.45, 0.0");
        scriptWriter.println("set size 0.45, 0.45");
        xlabel = "Original Distance";
        scriptWriter.println("set xlabel \"" + xlabel + "\" textcolor rgbcolor \"black\"");
        ylabel = "Count";
        scriptWriter.println("set ylabel \"" + ylabel + "\" textcolor rgbcolor \"black\"");
        title = "Histogram of " + xlabel;
        scriptWriter.println("set title \"" + title + "\" textcolor rgbcolor \"black\"");
        scriptWriter.println(
                "plot [:" + SALSAUtility.Xmaxbound + "] []'" + xHistFile.getFileName() + "' with filledcurves x1 lt rgb \"black\"");

        scriptWriter.println();

        scriptWriter.println("unset multiplot");
    }

    private static boolean IsOdd(int value) {
        return (value & 1) == 1;
    }

    private static void ReadClusterFile(Hashtable orignalPnumToCnumTable) {
        try (BufferedReader reader = Files.newBufferedReader(Paths.get(SALSAUtility.ClusterFile), Charset.defaultCharset())) {
            Pattern pattern = Pattern.compile("[\t ]");
            String line;
            while ((line = reader.readLine()) != null) {
                if (!Strings.isNullOrEmpty(line)) {
                    String[] splits = pattern.split(line.trim());
                    int pnum = Integer.parseInt(splits[0]);
                    int cnum = Integer.parseInt(splits[1]);
                    if (orignalPnumToCnumTable.contains(pnum)) {
                        throw new RuntimeException("Point numbers in the cluster file should be unique");
                    }
                    orignalPnumToCnumTable.put(pnum, cnum);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static boolean Calcfg(Desertwind Solution) throws MPIException {
        // Assume zerocr and first/DiagonalofMatrix are zeroed before call
        // Here we only calculate diagonal elements of Chisq matrix

        boolean violat = false;
        Hotsun.succ = true;
        Hotsun.zerocr = 0.0;
        Hotsun.idata = 0; // Not used in MDS

        GlobalReductions.FindDoubleSum Findzerocr = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);

        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
                    {
                        double WeightFunction1 = 0, WeightFunction2 = 0;
                        double localzerocr = 0.0;
                        int indexlen = SALSAUtility.PointsperThread[threadIndex];
                        int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;

                        for (int LocalToProcessIndex1 = beginpoint; LocalToProcessIndex1 < indexlen + beginpoint; LocalToProcessIndex1++) {
                            int GlobalIndex1 = LocalToProcessIndex1 + SALSAUtility.PointStart_Process;
                            int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalIndex1];
                            boolean vary1 = true;

                            if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] < SALSAUtility.SALSASHIFT) {
                                vary1 = false;
                            }

                            // Skip Fixed rows if StoredDistanceOption =3
                            if ((SALSAUtility.StoredDistanceOption == 3) && !vary1) {
                                continue;
                            }

                            for (int GlobalIndex2 = 0; GlobalIndex2 < SALSAUtility.PointCount_Global; GlobalIndex2++) {
                                if (GlobalIndex1 == GlobalIndex2) {
                                    continue;
                                }

                                int OriginalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalIndex2];
                                boolean vary2 = true;

                                if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] < SALSAUtility.SALSASHIFT) {
                                    vary2 = false;
                                }

                                // Varied (row) - Fixed (column) Contribution to Chisq Doubled if StoredDistanceOption 3 and
                                // so no row for fixed parameters
                                // vary1 forced to be true in this scenario
                                boolean symmetrize = false;

                                if ((SALSAUtility.StoredDistanceOption == 3) && !vary2) {
                                    symmetrize = true;
                                }

                                // If requested skip chisq for fixed cross fixed points
                                if (!SALSAUtility.CalcFixedCrossFixed) {
                                    if (!vary1 && !vary2) {
                                        continue;
                                    }
                                }

                                RefObject<Double> tempRef_WeightFunction1 = new RefObject<>(WeightFunction1);
                                RefObject<Double> tempRef_WeightFunction2 = new RefObject<>(WeightFunction2);
                                boolean tempVar = !SetChisqWeights(GlobalIndex1, GlobalIndex2, tempRef_WeightFunction1,
                                        tempRef_WeightFunction2);
                                WeightFunction1 = tempRef_WeightFunction1.argValue;
                                WeightFunction2 = tempRef_WeightFunction2.argValue;

                                if (tempVar) {
                                    continue;
                                }

                                // Calculate contribution to Chisq
                                double SquaredDistance = 0.0;

                                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < Hotsun.ParameterVectorDimension;
                                     LocalVectorIndex1++) {
                                    double tmp = Hotsun.GlobalParameter[GlobalIndex1][LocalVectorIndex1] - Hotsun
                                            .GlobalParameter[GlobalIndex2][LocalVectorIndex1];
                                    SquaredDistance += tmp * tmp;
                                }

                                double DistanceFudge1 = 1.0;
                                double DistanceFudge2 = 1.0;
                                double ActualDistance = SquaredDistance;

                                if (DistanceFormula == 1) {
                                    SquaredDistance += SQUAREMinimumDistance;
                                    ActualDistance = Math.sqrt(SquaredDistance);
                                    DistanceFudge1 = 0.5 / ActualDistance;
                                    DistanceFudge2 = DistanceFudge1 / SquaredDistance;
                                }

                                double funcl = WeightFunction1 - ActualDistance * WeightFunction2;
                                double symmetryweight = 1.0;

                                if (symmetrize) {
                                    symmetryweight = 2.0;
                                }

                                double increment = funcl * funcl * symmetryweight;
                                if (SALSAUtility.chisqcomponent > 0) {
                                    if ((SALSAUtility.chisqcomponent == 1) && (!vary1 || !vary2)) {
                                        continue;
                                    }
                                    if ((SALSAUtility.chisqcomponent == 2) && (!vary1 || vary2)) {
                                        continue;
                                    }
                                    if ((SALSAUtility.chisqcomponent == 3) && (vary1 || !vary2)) {
                                        continue;
                                    }
                                    if ((SALSAUtility.chisqcomponent == 4) && (vary1 || vary2)) {
                                        continue;
                                    }
                                }
                                localzerocr += increment;

                                if (SALSAUtility.chisqcomponent > 0) {
                                    continue;
                                }

                                if (!vary1) {
                                    continue;
                                }

                                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < Hotsun.ParameterVectorDimension; LocalVectorIndex1++) {
                                    double tmp1 = (Hotsun.GlobalParameter[GlobalIndex1][LocalVectorIndex1] - Hotsun
                                            .GlobalParameter[GlobalIndex2][LocalVectorIndex1]);

                                    //  Calculate First Derivative
                                    if (!Hotsun.FixedParameter[GlobalIndex1][LocalVectorIndex1]) {
                                        Solution.first[LocalToProcessIndex1][LocalVectorIndex1] += -4.0 * tmp1 * funcl * DistanceFudge1 * WeightFunction2;
                                    }

                                    // Calculate Diagonal Matrix Contributions to Second Derivative
                                    for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < Hotsun.ParameterVectorDimension; LocalVectorIndex2++) {
                                        double tmp2 = (Hotsun.GlobalParameter[GlobalIndex1][LocalVectorIndex2] - Hotsun
                                                .GlobalParameter[GlobalIndex2][LocalVectorIndex2]);
                                        double ApproximateCalc = 8.0 * tmp1 * tmp2 * DistanceFudge1 * DistanceFudge1 *
                                                WeightFunction2 * WeightFunction2;
                                        Solution.DiagonalofMatrix[LocalToProcessIndex1][LocalVectorIndex1][LocalVectorIndex2] += ApproximateCalc;
                                        double correction = 0.0;
                                        if ((DistanceFormula == 2) && (LocalVectorIndex1 == LocalVectorIndex2)) {
                                            correction = -4.0 * funcl * WeightFunction2;
                                        }
                                        if ((DistanceFormula == 1) && (LocalVectorIndex1 == LocalVectorIndex2)) {
                                            correction = -4.0 * funcl * WeightFunction2 * DistanceFudge1;
                                        }
                                        if (DistanceFormula == 1) {
                                            correction += 4.0 * funcl * tmp1 * tmp2 * DistanceFudge2 * WeightFunction2;
                                        }
                                        Solution.ExactDiagonalofMatrix[LocalToProcessIndex1][LocalVectorIndex1][LocalVectorIndex2] += ApproximateCalc + correction;
                                    }
                                }
                            }
                        }
                        Findzerocr.addAPoint(threadIndex, localzerocr);
                    }
            );
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

        Findzerocr.sumOverThreadsAndMPI();
        Hotsun.zerocr = Findzerocr.Total;
        Solution.Chisquared = Hotsun.zerocr;
        return violat;

    }


    //  DistanceProcessingOption = 1 and DistanceFormula = 1
    //  Chisqnorm = 0 ** (Euclidean Distance - Observed Distance) / Observed Distance
    //  Chisqnorm = 1 ** (Euclidean Distance - Observed Distance) / (Observed Distance)^0.75
    //  Chisqnorm = 2 ** (Euclidean Distance - Observed Distance) / (Observed Distance)^0.5     Sammon
    //  Chisqnorm = 3 ** (Euclidean Distance - Observed Distance)                               SMACOF
    //  Chisqnorm = 4 ** (Euclidean Distance - Observed Distance) * (Observed Distance)^0.5
    //
    //  DistanceProcessingOption = 2 and DistanceFormula = 2
    //  Chisqnorm = 0 ** (Euclidean Distance - Observed Distance) / Observed Distance
    //  Chisqnorm = 1 ** (Euclidean Distance - Observed Distance) / (Observed Distance)^0.5     Sammon
    //  Chisqnorm = 2 ** (Euclidean Distance - Observed Distance)                               SMACOF
    //  Chisqnorm = 3 ** (Euclidean Distance - Observed Distance) * Observed Distance
    //  Chisqnorm = 4 ** (Euclidean Distance - Observed Distance) * (Observed Distance)^2
    //
    //  DistanceProcessingOption = 1 and DistanceFormula = 2
    //  Chisqnorm = 0 ** (Euclidean Distance - SQRT(Observed Distance)) /  Observed Distance^0.5
    //  Chisqnorm = 1 ** (Euclidean Distance - SQRT(Observed Distance)) / (Observed Distance)^0.25
    //  Chisqnorm = 2 ** (Euclidean Distance - SQRT(Observed Distance))
    //  Chisqnorm = 3 ** (Euclidean Distance - SQRT(Observed Distance)) * (Observed Distance)^0.5
    //  Chisqnorm = 4 ** (Euclidean Distance - SQRT(Observed Distance)) *  Observed Distance

    //  Return false if interpoint distance not set
    public static boolean SetChisqWeights(int GlobalPointIndex1, int GlobalPointIndex2,
                                          RefObject<Double> WeightFunction1, RefObject<Double> WeightFunction2) {
        double InterPointDistance = SALSAParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
        if (InterPointDistance < -0.5) {
            WeightFunction1.argValue = 0.0;
            WeightFunction2.argValue = 0.0;
            return false;
        }

        double Weighting = 1.0;
        if (SALSAUtility.NumberDistanceWeightCuts > 0) {
            int distpos = SALSAUtility.NumberDistanceWeightCuts;
            for (int weightbin = 0; weightbin < SALSAUtility.NumberDistanceWeightCuts; weightbin++) {
                if (SALSAUtility.DistanceWeightCuts[weightbin] > InterPointDistance) {
                    distpos = weightbin;
                    break;
                }
            }
            Weighting = SALSAUtility.ActualWeightCuts[distpos];
        }

        if (WeightingOption > 0) {
            Weighting = Math.sqrt(PointWeights[GlobalPointIndex1] * PointWeights[GlobalPointIndex2]);
        }

        if (Chisqnorm == 0) {
            if (InterPointDistance > MinimumDistance) {
                WeightFunction1.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier;
                WeightFunction2.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier / InterPointDistance;
                return true;
            }
            WeightFunction2.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier / MinimumDistance;
            WeightFunction1.argValue = WeightFunction2.argValue * InterPointDistance;
            return true;
        }

        if (Chisqnorm == 1) { // Sammon if use Distance Squared with DistanceFormula == 2
            if (InterPointDistance > MinimumDistance) {
                double SQRTDistance = Math.sqrt(InterPointDistance);
                WeightFunction2.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier / (SQRTDistance * Math.sqrt(
                        SQRTSystemAverage * SQRTDistance));
                WeightFunction1.argValue = WeightFunction2.argValue * InterPointDistance;
                return true;
            }
            WeightFunction2.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier / (SQRTMinimumDistance * Math.sqrt(
                    SQRTSystemAverage * SQRTMinimumDistance));
            WeightFunction1.argValue = WeightFunction2.argValue * InterPointDistance;
            return true;
        }

        if (Chisqnorm == 2) { // SMACOF if use Distance Squared  with DistanceFormula == 2
            // Sammon if use NonSquared Distance with DistanceFormula == 1

            if (InterPointDistance > MinimumDistance) {
                double SQRTDistance = Math.sqrt(InterPointDistance);
                WeightFunction2.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier / (SQRTDistance * SQRTSystemAverage);
                WeightFunction1.argValue = WeightFunction2.argValue * InterPointDistance;
                return true;
            }
            WeightFunction2.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier / (SQRTMinimumDistance * SQRTSystemAverage);
            WeightFunction1.argValue = WeightFunction2.argValue * InterPointDistance;
            return true;
        }

        if (Chisqnorm == 3) { // SMACOF if use Nonsquared Distance   with DistanceFormula == 1
            WeightFunction2.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier / SystemAverage;
            WeightFunction1.argValue = WeightFunction2.argValue * InterPointDistance;
            return true;
        }

        if (Chisqnorm == 4) { // Even more large distance weight than SMACOF if DistanceFormula == 1
            WeightFunction2.argValue = Weighting * ManxcatCentral.ChisqFunctionCalcMultiplier * Math.sqrt(
                    InterPointDistance) / (SQRTSystemAverage * SystemAverage);
            WeightFunction1.argValue = WeightFunction2.argValue * InterPointDistance;
            return true;
        }

        //  Default
        WeightFunction1.argValue = Weighting;
        WeightFunction2.argValue = Weighting;
        return true;

    }

    //  Set initial value of param in Solution
    //  InitializationOption =0 None (except fixed points) =1 Old File Type Model =2 New File Type
    //  InitializationOption 1 must have a file with length identical to number of used points
    //  For fixed points, any values read in with InitializationOption 1 or 2 overwrite those read in SALSAProcessVariedandFixed
    //  Such overwriting must happen for InitializationOption 1 but only happens for InitializationOption 2 if fixed points read
    public static void InitializeParameters(Desertwind Solution, int CountStartingPoints) {

        double RadiusUsed = Math.sqrt(SystemRadius);
        java.util.Random Randobject = new java.util.Random();
        String[] InitialMDSString = new String[SALSAUtility.PointCount_Global];
        String OriginalMDSFileName = ManxcatCentral.config.InitializationFileName;

        if (!OriginalMDSFileName.contains(":") && !OriginalMDSFileName.contains("$")) {
            OriginalMDSFileName = ManxcatCentral.config.ControlDirectoryName + File.separatorChar + OriginalMDSFileName;
        }
        int InitializationNumberofPoints = 0;

        if (ManxcatCentral.config.InitializationOption == 1) { // This file must be same length as number of used points

            int[] ColorValue = new int[SALSAUtility.PointCount_Global];
            RefObject<Integer> tempRef_InitializationNumberofPoints = new RefObject<Integer>(
                    InitializationNumberofPoints);
            boolean tempVar = ManxcatMDSDataProcessing.ReadMDSCluster_File(OriginalMDSFileName, InitialMDSString,
                    ColorValue, tempRef_InitializationNumberofPoints);
            InitializationNumberofPoints = tempRef_InitializationNumberofPoints.argValue;
            if (tempVar) {
                System.out.println(OriginalMDSFileName + " Read Successfully Points " + InitializationNumberofPoints);
            }
            if (SALSAUtility.PointCount_Global != InitializationNumberofPoints) {
                SALSAUtility.printAndThrowRuntimeException(
                        " Inconsistent Initialization File Point Counts " + InitializationNumberofPoints + " Expected" +
                                " is " + SALSAUtility.PointCount_Global
                );

            }
        }

        if (ManxcatCentral.config.InitializationOption == 2) {
            int InitializationFileType = -1;
            SALSAFileProperties InitializationFileProperties = new SALSAFileProperties();
            SALSADataPointProperties[] InitializationPointProperties = new SALSADataPointProperties[SALSAUtility.NumberOriginalPoints];
            RefObject<Integer> tempRef_InitializationFileType = new RefObject<Integer>(InitializationFileType);
            RefObject<SALSADataPointProperties[]> tempRef_InitializationPointProperties = new RefObject<SALSADataPointProperties[]>(
                    InitializationPointProperties);
            RefObject<Integer> tempRef_InitializationNumberofPoints2 = new RefObject<Integer>(
                    InitializationNumberofPoints);
            SALSA_Properties.ReadDataPointFile(OriginalMDSFileName, tempRef_InitializationFileType,
                    InitializationFileProperties, tempRef_InitializationPointProperties,
                    tempRef_InitializationNumberofPoints2);
            InitializationFileType = tempRef_InitializationFileType.argValue;
            InitializationPointProperties = tempRef_InitializationPointProperties.argValue;
            InitializationNumberofPoints = tempRef_InitializationNumberofPoints2.argValue;

            if ((SALSAUtility.NumberOriginalPoints < InitializationNumberofPoints) || (InitializationFileProperties.NumberOriginalPoints != SALSAUtility.NumberOriginalPoints)) {
                SALSAUtility.printAndThrowRuntimeException(
                        " Inconsistent Initialization File Point Counts " + InitializationNumberofPoints + " or " +
                                InitializationFileProperties.NumberOriginalPoints + " Expected is " + SALSAUtility
                                .NumberOriginalPoints
                );

            }

            for (int InitialIndex = 0; InitialIndex < InitializationNumberofPoints; InitialIndex++) {
                int OriginalIndex = InitializationPointProperties[InitialIndex].OriginalPointNumber;
                if (!InitializationPointProperties[InitialIndex].valuesset) {
                    continue;
                }

                // Consider all used points
                int UsedIndex = SALSAUtility.OriginalPointtoUsedPointMap[OriginalIndex];
                if (UsedIndex < 0) {
                    continue;
                }

                SALSAUtility.GlobalPointProperties[UsedIndex].x = InitializationPointProperties[InitialIndex].x;
                SALSAUtility.GlobalPointProperties[UsedIndex].y = InitializationPointProperties[InitialIndex].y;
                SALSAUtility.GlobalPointProperties[UsedIndex].z = InitializationPointProperties[InitialIndex].z;
                SALSAUtility.GlobalPointProperties[UsedIndex].valuesset = true;
                SALSAUtility.GlobalPointProperties[UsedIndex].source = InitializationPointProperties[InitialIndex].source;

                if (InitializationPointProperties[InitialIndex].errorsset) {
                    SALSAUtility.GlobalPointProperties[UsedIndex].xerr = InitializationPointProperties[InitialIndex].xerr;
                    SALSAUtility.GlobalPointProperties[UsedIndex].yerr = InitializationPointProperties[InitialIndex].yerr;
                    SALSAUtility.GlobalPointProperties[UsedIndex].zerr = InitializationPointProperties[InitialIndex].zerr;
                    SALSAUtility.GlobalPointProperties[UsedIndex].errorsset = true;
                }
                SALSAUtility.GlobalPointProperties[UsedIndex].family1 = InitializationPointProperties[InitialIndex].family1;
                SALSAUtility.GlobalPointProperties[UsedIndex].family2 = InitializationPointProperties[InitialIndex].family2;
                SALSAUtility.GlobalPointProperties[UsedIndex].cluster = InitializationPointProperties[InitialIndex].cluster;
                SALSAUtility.GlobalPointProperties[UsedIndex].familylabel1 = InitializationPointProperties[InitialIndex].familylabel1;
                SALSAUtility.GlobalPointProperties[UsedIndex].familylabel2 = InitializationPointProperties[InitialIndex].familylabel2;
                SALSAUtility.GlobalPointProperties[UsedIndex].clusterlabel = InitializationPointProperties[InitialIndex].clusterlabel;
                SALSAUtility.GlobalPointProperties[UsedIndex].pointlabel = InitializationPointProperties[InitialIndex].pointlabel;
                SALSAUtility.GlobalPointProperties[UsedIndex].PointType = InitializationPointProperties[InitialIndex].PointType;
            }
        }

        final Pattern pattern = Pattern.compile("[\t ]");
        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
                    {
                        int indexlen = SALSAUtility.PointsperThread[threadIndex];
                        int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                        for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++) {
                            int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;
                            if (ManxcatCentral.config.InitializationOption == 1) {
                                for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++) {
                                    String[] split = pattern.split(
                                            InitialMDSString[SALSAUtility.ActualtoNaiveUsedOrder[GlobalIndex]]);
                                    Solution.param[LongIndex][LocalVectorIndex] = Double.parseDouble(
                                            split[LocalVectorIndex + 1]);
                                }
                            }
                            if (ManxcatCentral.config.InitializationOption == 2) {
                                if (!SALSAUtility.GlobalPointProperties[GlobalIndex].valuesset) {
                                    SALSAUtility.printAndThrowRuntimeException(
                                            "Error in Initialized Points -- Unset point " + GlobalIndex);
                                }
                                Solution.param[LongIndex][0] = SALSAUtility.GlobalPointProperties[GlobalIndex].x;
                                Solution.param[LongIndex][1] = SALSAUtility.GlobalPointProperties[GlobalIndex].y;
                                if (Hotsun.ParameterVectorDimension > 2) {
                                    Solution.param[LongIndex][2] = SALSAUtility.GlobalPointProperties[GlobalIndex].z;
                                }
                            }
                            if (ManxcatCentral.config.InitializationOption == 0) {
                                double bigorsmall = 2.0 * Randobject.nextDouble();
                                for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++) {
                                    if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex]) {
                                        if (FindPointstoFix) {
                                            Solution.param[LongIndex][LocalVectorIndex] = 0.0;
                                        } else {
                                            if (!SALSAUtility.GlobalPointProperties[GlobalIndex].valuesset) {
                                                SALSAUtility.printAndThrowRuntimeException(
                                                        "Error in Fixed Points -- Unset point " + GlobalIndex);
                                            }
                                            if (LocalVectorIndex == 0) {
                                                Solution.param[LongIndex][LocalVectorIndex] = SALSAUtility.GlobalPointProperties[GlobalIndex].x;
                                            }
                                            if (LocalVectorIndex == 1) {
                                                Solution.param[LongIndex][LocalVectorIndex] = SALSAUtility.GlobalPointProperties[GlobalIndex].y;
                                            }
                                            if (LocalVectorIndex == 2) {
                                                Solution.param[LongIndex][LocalVectorIndex] = SALSAUtility.GlobalPointProperties[GlobalIndex].z;
                                            }
                                        }
                                    } else {
                                        Solution.param[LongIndex][LocalVectorIndex] = bigorsmall * Randobject.nextDouble() * RadiusUsed;
                                    }
                                }
                            }
                        }
                    }
            );
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }

    }

    public static void SetupWeightings(double[] WeightsasRead) {
        if (WeightingOption == 0) {
            for (int Globalindex = 0; Globalindex < SALSAUtility.PointCount_Global; Globalindex++) {
                WeightsasRead[Globalindex] = 1.0;
            }
            return;
        }

        String WeightFileName = ManxcatCentral.config.WeightingFileName;

        double sumofweights = 0.0;
        int NumberofLines = 0;

        Path weightFilePath = Paths.get(WeightFileName);
        if (!Files.exists(weightFilePath)) {
            SALSAUtility.printAndThrowRuntimeException("File " + WeightFileName + " does not exists.");
        } else {
            try (BufferedReader reader = Files.newBufferedReader(weightFilePath, Charset.defaultCharset())) {

                // Read contents of a file, line by line, into a string
                String inputLineStr;
                while ((inputLineStr = reader.readLine()) != null) {
                    inputLineStr = inputLineStr.trim();

                    if (inputLineStr.length() < 1) {
                        continue; //replace empty line
                    }

                    WeightsasRead[NumberofLines] = Double.parseDouble(inputLineStr);
                    sumofweights += WeightsasRead[NumberofLines];
                    ++NumberofLines;
                }
                reader.close();
            } catch (IOException e) {
                SALSAUtility.printAndThrowRuntimeException(
                        "Failed to read data from " + WeightFileName + " " + e.toString());
            }
        }


        if (NumberofLines != SALSAUtility.PointCount_Global) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Incorrect Weight count read " + NumberofLines + " Expected " + SALSAUtility.PointCount_Global);
        }
        double AveragetoOne = NumberofLines / sumofweights;
        double minweight = sumofweights;
        double maxweight = 0.0;

        for (int Globalindex = 0; Globalindex < SALSAUtility.PointCount_Global; Globalindex++) {
            WeightsasRead[Globalindex] *= AveragetoOne;
            minweight = Math.min(minweight, WeightsasRead[Globalindex]);
            maxweight = Math.max(maxweight, WeightsasRead[Globalindex]);
        }
        SALSAUtility.SALSAPrint(1,
                "File " + WeightFileName + " Non trivial Point Weights Maximum " + String.format("%.3f",
                        maxweight) + " Minimum " + String.format("%.4E", minweight)
        );
    }

    //  Calculate Matrix Global Vector product storing as a distributed vector
    public static void GlobalMatrixVectorProduct(double[][] DistributedVector, Desertwind Solution, boolean useexact,
                                                 double[][] GlobalxVector, double[][] GlobalVectoronRight) {
        double[][][] MatrixDiagonals;

        if (useexact) {
            MatrixDiagonals = Solution.ExactDiagonalofMatrix;
        } else {
            MatrixDiagonals = Solution.DiagonalofMatrix;
        }

        try {
            forallChunked(0, SALSAUtility.ThreadCount - 1, (threadIndex) ->
                    {
                        double WeightFunction1 = 0, WeightFunction2 = 0;
                        int indexlen = SALSAUtility.PointsperThread[threadIndex];
                        int beginpoint = SALSAUtility.StartPointperThread[threadIndex] - SALSAUtility.PointStart_Process;
                        for (int LocalProcessIndex1 = beginpoint; LocalProcessIndex1 < indexlen + beginpoint; LocalProcessIndex1++) {
                            int GlobalIndex1 = LocalProcessIndex1 + SALSAUtility.PointStart_Process;
                            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < Hotsun.ParameterVectorDimension; LocalVectorIndex1++) {
                                if (Hotsun.FixedParameter[GlobalIndex1][LocalVectorIndex1]) {
                                    DistributedVector[LocalProcessIndex1][LocalVectorIndex1] = 0.0;
                                } else {
                                    double tmp = 0.0;
                                    for (int GlobalIndex2 = 0; GlobalIndex2 < SALSAUtility.PointCount_Global; GlobalIndex2++) {
                                        int OriginalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalIndex2];
                                        if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] < SALSAUtility.SALSASHIFT) {
                                            continue;
                                        }
                                        RefObject<Double> tempRef_WeightFunction1 = new RefObject<>(WeightFunction1);
                                        RefObject<Double> tempRef_WeightFunction2 = new RefObject<>(WeightFunction2);
                                        boolean tempVar = !SetChisqWeights(GlobalIndex1, GlobalIndex2,
                                                tempRef_WeightFunction1, tempRef_WeightFunction2);
                                        WeightFunction1 = tempRef_WeightFunction1.argValue;
                                        WeightFunction2 = tempRef_WeightFunction2.argValue;
                                        if (tempVar) {
                                            continue;
                                        }
                                        double DistanceFudge1 = 1.0;
                                        double DistanceFudge2 = 1.0;
                                        double SquaredDistance = 0.0;
                                        double funcl = 0.0;
                                        if (GlobalIndex1 != GlobalIndex2) {
                                            for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < Hotsun.ParameterVectorDimension; LocalVectorIndex2++) {
                                                double tmp1 = GlobalxVector[GlobalIndex1][LocalVectorIndex2] - GlobalxVector[GlobalIndex2][LocalVectorIndex2];
                                                SquaredDistance += tmp1 * tmp1;
                                            }
                                            double ActualDistance = SquaredDistance;
                                            if (DistanceFormula == 1) {
                                                SquaredDistance += SQUAREMinimumDistance;
                                                ActualDistance = Math.sqrt(SquaredDistance);
                                                DistanceFudge1 = 0.5 / ActualDistance;
                                                DistanceFudge2 = DistanceFudge1 / SquaredDistance;
                                            }
                                            funcl = WeightFunction1 - WeightFunction2 * ActualDistance;
                                        }
                                        for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < Hotsun.ParameterVectorDimension; LocalVectorIndex2++) {
                                            if (Hotsun.FixedParameter[GlobalIndex2][LocalVectorIndex2]) {
                                                continue;
                                            }
                                            double MatrixElement;
                                            if (GlobalIndex1 == GlobalIndex2) {
                                                MatrixElement = MatrixDiagonals[LocalProcessIndex1][LocalVectorIndex1][LocalVectorIndex2];
                                                if (Hotsun.UseDiagonalScaling) {
                                                    MatrixElement *= Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex1] * Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex2];
                                                }
                                                if (Hotsun.AddMarquardtQDynamically && (LocalVectorIndex1 == LocalVectorIndex2)) {
                                                    MatrixElement += Hotsun.Q;
                                                }
                                            } else {
                                                double correction = 0.0;
                                                double VectorCrossProduct = (GlobalxVector[GlobalIndex1][LocalVectorIndex1] - GlobalxVector[GlobalIndex2][LocalVectorIndex1]) * (GlobalxVector[GlobalIndex1][LocalVectorIndex2] - GlobalxVector[GlobalIndex2][LocalVectorIndex2]);
                                                MatrixElement = -8.0 * VectorCrossProduct * DistanceFudge1 * DistanceFudge1 * WeightFunction2 * WeightFunction2;
                                                if (Hotsun.FullSecondDerivative) {
                                                    if ((DistanceFormula == 2) && (LocalVectorIndex1 == LocalVectorIndex2)) {
                                                        correction = 4.0 * funcl * WeightFunction2;
                                                    }
                                                    if ((DistanceFormula == 1) && (LocalVectorIndex1 == LocalVectorIndex2)) {
                                                        correction = 4.0 * funcl * WeightFunction2 * DistanceFudge1;
                                                    }
                                                    if (DistanceFormula == 1) {
                                                        correction += -4.0 * funcl * VectorCrossProduct * WeightFunction2 * DistanceFudge2;
                                                    }
                                                    MatrixElement += correction;
                                                }
                                                if (Hotsun.UseDiagonalScaling) {
                                                    MatrixElement *= Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex1] * Hotsun.sqdginv[GlobalIndex2][LocalVectorIndex2];
                                                }
                                            }
                                            tmp += MatrixElement * GlobalVectoronRight[GlobalIndex2][LocalVectorIndex2];
                                        }
                                    }
                                    DistributedVector[LocalProcessIndex1][LocalVectorIndex1] = tmp;
                                }
                            }
                        }
                    }
            );
        } catch (SuspendableException e) {
            SALSAUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }


    // Solve Matrix Equations
    //  The solution is rescaled so correct for native matrix
    //  However Diagonally scaled matrix used in solver as long as Hotsun.UseDiagonalScalinginSolvers = true
    public static boolean SolveMatrix(double[][] Answer, Desertwind Solution) throws MPIException {
        // Scale RHS Vector
        MDSLinearAlgebra.DiagScaleVector(Hotsun.UtilityLocalVector1, Solution.first, Hotsun.sqdginv);

        //  Solve Scaled Matrix Equations
        int RealMatrixSize = Hotsun.npar - ((Hotsun.ParameterVectorDimension + 1) * Hotsun.ParameterVectorDimension) / 2;
        RealMatrixSize = Math.min(RealMatrixSize, Hotsun.CGIterationLimit);
        boolean matrixsuccess;
        RefObject<Integer> tempRef_NumberofCGIterations = new RefObject<Integer>(Hotsun.NumberofCGIterations);
        matrixsuccess = MDSLinearAlgebra.ConjugateGradientSolver(Answer, Solution, Hotsun.FullSecondDerivative,
                Hotsun.GlobalParameter, Hotsun.UtilityLocalVector1, tempRef_NumberofCGIterations, RealMatrixSize,
                Hotsun.CGResidualLimit);
        Hotsun.NumberofCGIterations = tempRef_NumberofCGIterations.argValue;

        if (!matrixsuccess) {
            ++Hotsun.TotalCGFailures;
        }
        Hotsun.TotalCGIterations += Hotsun.NumberofCGIterations;


        //  Correct answer for diagonal scaling
        MDSLinearAlgebra.DiagScaleVector(Answer, Answer, Hotsun.sqdginv);
        return matrixsuccess;

    }

    public static void FindQlimits(Desertwind Solution, RefObject<Double> Qhigh, RefObject<Double> Qlow,
                                   RefObject<Integer> ReasontoStop1, RefObject<Integer> ReasontoStop2) throws MPIException {
        if (Hotsun.FullSecondDerivative) {
            RefObject<Double> tempRef_ChisqMatrixTrace = new RefObject<>(Hotsun.ChisqMatrixTrace);
            RefObject<Double> tempRef_ChisqMatrixNorm = new RefObject<>(Hotsun.ChisqMatrixNorm);
            MDSLinearAlgebra.FindTraceandNorm(Solution.ExactDiagonalofMatrix, Hotsun.GlobalParameter,
                    tempRef_ChisqMatrixTrace, tempRef_ChisqMatrixNorm);
            Hotsun.ChisqMatrixTrace = tempRef_ChisqMatrixTrace.argValue;
            Hotsun.ChisqMatrixNorm = tempRef_ChisqMatrixNorm.argValue;
        } else {
            RefObject<Double> tempRef_ChisqMatrixTrace2 = new RefObject<Double>(Hotsun.ChisqMatrixTrace);
            RefObject<Double> tempRef_ChisqMatrixNorm2 = new RefObject<Double>(Hotsun.ChisqMatrixNorm);
            MDSLinearAlgebra.FindTraceandNorm(Solution.DiagonalofMatrix, Hotsun.GlobalParameter,
                    tempRef_ChisqMatrixTrace2, tempRef_ChisqMatrixNorm2);
            Hotsun.ChisqMatrixTrace = tempRef_ChisqMatrixTrace2.argValue;
            Hotsun.ChisqMatrixNorm = tempRef_ChisqMatrixNorm2.argValue;
        }

        Qhigh.argValue = 0.0;
        Qlow.argValue = 0.0;
        double PowerEigenvalue1 = 0.0;
        ReasontoStop2.argValue = -3;

        RefObject<Double> tempRef_PowerEigenvalue1 = new RefObject<Double>(PowerEigenvalue1);
        ReasontoStop1.argValue = MDSLinearAlgebra.PowerIterate(Solution, 0, 0.0, tempRef_PowerEigenvalue1);
        PowerEigenvalue1 = tempRef_PowerEigenvalue1.argValue;

        if (ReasontoStop1.argValue > 0) {
            Hotsun.TotalPowerIterations += ReasontoStop1.argValue;
        } else {
            ++Hotsun.TotalPowerFailures;
        }

        if (ReasontoStop1.argValue > 0) {
            Qhigh.argValue = PowerEigenvalue1 - Hotsun.addonforQcomputation;
        } else {
            Qhigh.argValue = Math.min(Hotsun.ChisqMatrixNorm, Hotsun.ChisqMatrixTrace);
            Qlow.argValue = Hotsun.Qrange * Qhigh.argValue;
            return;
        }

        double PowerEigenvalue2 = 0.0;
        RefObject<Double> tempRef_PowerEigenvalue2 = new RefObject<Double>(PowerEigenvalue2);
        ReasontoStop2.argValue = MDSLinearAlgebra.PowerIterate(Solution, 1, Qhigh.argValue, tempRef_PowerEigenvalue2);
        PowerEigenvalue2 = tempRef_PowerEigenvalue2.argValue;

        if (ReasontoStop2.argValue > 0) {
            Hotsun.TotalPowerIterations += ReasontoStop2.argValue;
        } else {
            Hotsun.TotalPowerFailures += 1;
        }

        if (ReasontoStop2.argValue > 0) {
            Qlow.argValue = Qhigh.argValue - PowerEigenvalue2;
        } else {
            Qlow.argValue = Hotsun.Qrange * Qhigh.argValue;
        }
        boolean Qlowtoosmall = Qlow.argValue < (Hotsun.Qrange * Qhigh.argValue);
        Qlowtoosmall = SALSAUtility.synchronizeMPIVariable(Qlowtoosmall);

        if (Qlowtoosmall) {
            Qlow.argValue = Hotsun.Qrange * Qhigh.argValue;
        }
    }

    public static void FillupHotsun() {
        if (ManxcatCentral.config.InitializationOption != 1 || !(new java.io.File(
                ManxcatCentral.config.InitializationFileName)).isFile()) {
            throw new RuntimeException("Simple initialization file necessary with the specified processing option");
        }
        try (BufferedReader reader = Files.newBufferedReader(Paths.get(ManxcatCentral.config.InitializationFileName), Charset.defaultCharset())) {
            Pattern pattern = Pattern.compile("[\t ]");
            Hotsun.GlobalParameter = new double[SALSAUtility.PointCount_Global][];
            String line;
            while ((line = reader.readLine()) != null) {
                if (!Strings.isNullOrEmpty(line)) {
                    String[] splits = pattern.split(line.trim());
                    if (splits.length == 5) {
                        int originalIndex = Integer.parseInt(splits[0]);
                        double[] vector = new double[]{Double.parseDouble(splits[1]), Double.parseDouble(
                                splits[2]), Double.parseDouble(splits[3])};
                        if (originalIndex >= 0 && originalIndex < SALSAUtility.OriginalPointtoUsedPointMap.length) {
                            int globalIndex = SALSAUtility.OriginalPointtoUsedPointMap[originalIndex];
                            int usedIndex = SALSAUtility.NaivetoActualUsedOrder[globalIndex];
                            Hotsun.GlobalParameter[usedIndex] = vector;
                        }
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}