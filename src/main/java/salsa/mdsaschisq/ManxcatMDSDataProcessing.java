package salsa.mdsaschisq;

import com.google.common.base.Strings;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;
import tangible.RefObject;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Date;
import java.util.regex.Pattern;

public class ManxcatMDSDataProcessing {
    public static void UpdateManxcatMDS_Option12() {
        if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("traditionaldirectory")) {
            UpdateManxcatMDS_Option12_TraditionalDirectory();
            return;
        }

        if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("cluster")) {
            UpdateManxcatMDS_Option12_Cluster();
            return;
        }

        if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("cutfile")) {
            UpdateManxcatMDS_Option12_FamilybyCuts();
            return;
        }

        if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("function")) {
            if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("statistics")) {
                UpdateManxcatMDS_Option12_Functions(true);
            } else {
                UpdateManxcatMDS_Option12_Functions(false);
            }
        }
    } // End  UpdateManxcatMDS_Option12()

    public static void UpdateManxcatMDS_Option12_Functions(boolean statistics) {
        //  Setup weight: Only used in statistics
        double[] weights = new double[SALSAUtility.PointCount_Global];
        ManxcatMDS.WeightingOption = ManxcatCentral.config.WeightingOption;
        ManxcatMDS.SetupWeightings(weights);
        String TypicalMDSFileName = "";

        //  Setup MDS: Only used in statistics
        double[][] MDSPoints = new double[SALSAUtility.PointCount_Global][Hotsun.ParameterVectorDimension];

        if (statistics) { // Set up statistics by reading MDS file
            TypicalMDSFileName = ManxcatCentral.ResultDirectoryName + "\\SIMPLE" + ManxcatCentral.config
                    .ReducedVectorOutputFileName;

            if (!(new java.io.File(TypicalMDSFileName)).isFile()) {
                SALSAUtility.printAndThrowRuntimeException(" File " + TypicalMDSFileName + " Does Not Exist");


            }

            int NumberMDSPoints = 0;
            String[] InitialMDSString = new String[SALSAUtility.PointCount_Global];
            int[] ColorValue = new int[SALSAUtility.PointCount_Global];

            RefObject<Integer> tempRef_NumberMDSPoints = new RefObject<Integer>(NumberMDSPoints);
            boolean tempVar = ReadMDSCluster_File(TypicalMDSFileName, InitialMDSString, ColorValue,
                    tempRef_NumberMDSPoints);
            NumberMDSPoints = tempRef_NumberMDSPoints.argValue;
            if (tempVar) {
                SALSAUtility.SALSAPrint(0, "MDS File " + TypicalMDSFileName + " Read Successfully");
            }

            Pattern pattern = Pattern.compile("[\t ]");
            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) { // Extract MDS POints
                String inputLineStr = InitialMDSString[PointIndex].trim();
                String[] split = pattern.split(inputLineStr);

                if (split.length < Hotsun.ParameterVectorDimension) {
                    SALSAUtility.printAndThrowRuntimeException(
                            " Point " + PointIndex + " has wrong number of Entries " + inputLineStr + " Entries" +
                                    split.length);


                }

                for (int VectorIndex = 0; VectorIndex < Hotsun.ParameterVectorDimension; VectorIndex++) {
                    MDSPoints[PointIndex][VectorIndex] = Double.parseDouble(split[VectorIndex].trim());
                }
            }
        } // End set up of Statistics

        //  Set Mapping
        int[] ClusterLabel = new int[SALSAUtility.PointCount_Global];
        int NumDivisions = 10;
        int BasicSize = SALSAUtility.PointCount_Global / NumDivisions;
        int LeftOver = SALSAUtility.PointCount_Global - BasicSize * NumDivisions;
        int DivisionCount = 0;
        int Limit = 0;

        for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
            if (PointIndex >= Limit) {
                ++DivisionCount;
                Limit += BasicSize;

                if (LeftOver >= DivisionCount) {
                    ++Limit;
                }
            }
            ClusterLabel[PointIndex] = DivisionCount;
        }

        String FunctionFileName = ManxcatCentral.config.ConversionInformation;

        if (!FunctionFileName.contains(":")) {
            FunctionFileName = ManxcatCentral.config.ControlDirectoryName + File.separatorChar + FunctionFileName;
        }

        FunctionKeyLabelPair[] functionKeysAndLabels = new FunctionKeyLabelPair[SALSAUtility.PointCount_Global];

        int pickout = -1;

        while (true) {
            ++pickout;
            boolean endofdata = true;
            int NumberofLines = 0;
            double sumabs = 0.0;
            double totalweight = 0.0;
            double meangamma = 0.0;
            double vargamma = 0.0;
            double[] meanxyz = new double[Hotsun.ParameterVectorDimension];
            double[][] varxyz = new double[Hotsun.ParameterVectorDimension][Hotsun.ParameterVectorDimension];
            double[] correlxyzgamma = new double[Hotsun.ParameterVectorDimension];

            if (statistics) {
                for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++) {
                    meanxyz[VectorIndex1] = 0.0;

                    for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++) {
                        varxyz[VectorIndex1][VectorIndex2] = 0.0;
                    }
                    correlxyzgamma[VectorIndex1] = 0.0;
                }
            }

            Path functionFilePath = Paths.get(FunctionFileName);
            if (!Files.exists(functionFilePath)) {
                SALSAUtility.printAndThrowRuntimeException("File " + FunctionFileName + " does not exists.");
            } else {
                try (BufferedReader reader = Files.newBufferedReader(functionFilePath)) {
                    Pattern pattern = Pattern.compile("[ ]");
                    // Read contents of a file, line by line, into a string
                    String inputLineStr;

                    while ((inputLineStr = reader.readLine()) != null) {
                        inputLineStr = inputLineStr.trim();
                        String[] split = pattern.split(inputLineStr);

                        if (split.length <= pickout) {
                            break;
                        }
                        endofdata = false;
                        String usethis = split[pickout].trim();

                        if (usethis.length() < 1) {
                            continue; //replace empty line
                        }
                        double gamma = Double.parseDouble(usethis);
                        functionKeysAndLabels[NumberofLines].setKey(gamma);
                        functionKeysAndLabels[NumberofLines].setLabel(NumberofLines);
                        sumabs += Math.abs(gamma);

                        if (statistics) {
                            double wgt = weights[NumberofLines];
                            meangamma += gamma * wgt;
                            vargamma += wgt * gamma * gamma;
                            totalweight += wgt;

                            for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++) {
                                meanxyz[VectorIndex1] += MDSPoints[NumberofLines][VectorIndex1] * wgt;
                                correlxyzgamma[VectorIndex1] += MDSPoints[NumberofLines][VectorIndex1] * gamma * wgt;

                                for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension;
                                     VectorIndex2++) {
                                    varxyz[VectorIndex1][VectorIndex2] += MDSPoints[NumberofLines][VectorIndex1] *
                                            MDSPoints[NumberofLines][VectorIndex2] * wgt;
                                }
                            }
                        }
                        ++NumberofLines;
                    }
                    reader.close();
                } catch (IOException e) {
                    SALSAUtility.printAndThrowRuntimeException(
                            "Failed to read data from " + FunctionFileName + " " + e.toString());
                }
            }


            if (endofdata) {
                break;
            }

            if (NumberofLines != SALSAUtility.PointCount_Global) {
                SALSAUtility.printAndThrowRuntimeException(
                        "Incorrect Function count read " + NumberofLines + " Expected " + SALSAUtility
                                .PointCount_Global);
            }

            // Set Statistics
            if (statistics) {
                double[] alpha = new double[Hotsun.ParameterVectorDimension];
                double alphanorm = 0.0;
                meangamma = meangamma / totalweight;
                vargamma = (vargamma / totalweight) - meangamma * meangamma;

                for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++) {
                    meanxyz[VectorIndex1] = meanxyz[VectorIndex1] / totalweight;
                }

                for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++) {
                    for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++) {
                        varxyz[VectorIndex1][VectorIndex2] = (varxyz[VectorIndex1][VectorIndex2] / totalweight) -
                                meanxyz[VectorIndex1] * meanxyz[VectorIndex2];
                    }
                    correlxyzgamma[VectorIndex1] = correlxyzgamma[VectorIndex1] / totalweight - meangamma *
                            meanxyz[VectorIndex1];
                    alpha[VectorIndex1] = correlxyzgamma[VectorIndex1] / varxyz[VectorIndex1][VectorIndex1];
                    alphanorm += alpha[VectorIndex1] * alpha[VectorIndex1];
                }

                // invert Matrix to find best alpha
                double[][] ConventionalFirst = new double[Hotsun.ParameterVectorDimension][1];
                double[][] ConventionalAnswer = new double[Hotsun.ParameterVectorDimension][];
                double[][] ConventionalMatrix = new double[Hotsun.ParameterVectorDimension][Hotsun
                        .ParameterVectorDimension];

                for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++) {
                    ConventionalAnswer[VectorIndex1] = new double[1];
                    ConventionalFirst[VectorIndex1][0] = correlxyzgamma[VectorIndex1] / Math.sqrt(
                            varxyz[VectorIndex1][VectorIndex1]);

                    for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++) {
                        ConventionalMatrix[VectorIndex1][VectorIndex2] = varxyz[VectorIndex1][VectorIndex2] / Math.sqrt(
                                varxyz[VectorIndex1][VectorIndex1] * varxyz[VectorIndex2][VectorIndex2]);
                    }
                }
                DoubleMatrix cMatrix = new DoubleMatrix(ConventionalMatrix);
                DoubleMatrix RightHandSide = new DoubleMatrix(ConventionalFirst);
                ConventionalAnswer = Solve.solve(cMatrix, RightHandSide).toArray2();

                alphanorm = 0.0;

                for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++) {
                    alpha[VectorIndex1] = ConventionalAnswer[VectorIndex1][0] / Math.sqrt(
                            varxyz[VectorIndex1][VectorIndex1]);
                    alphanorm += alpha[VectorIndex1] * alpha[VectorIndex1];
                }

                double Fullcorrelation = 0.0;
                double varalphaxyz = 0.0;

                for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++) {
                    alpha[VectorIndex1] = alpha[VectorIndex1] / Math.sqrt(alphanorm);
                    Fullcorrelation += alpha[VectorIndex1] * correlxyzgamma[VectorIndex1];
                }

                for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++) {
                    for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++) {
                        varalphaxyz += alpha[VectorIndex1] * alpha[VectorIndex2] * varxyz[VectorIndex1][VectorIndex2];
                    }
                }
                System.out.println(varalphaxyz + " " + Fullcorrelation + " " + vargamma);
                Fullcorrelation = Fullcorrelation / (Math.sqrt(vargamma * varalphaxyz));
                SALSAUtility.SALSAPrint(0,
                        "Column " + pickout + " File " + FunctionFileName + " MDS File " + TypicalMDSFileName + " " +
                                "Total Weight " + totalweight + " Correlation " + String
                                .format("%.4f", Fullcorrelation) + " Direction " + String
                                .format("%.4f", alpha[0]) + " " + String
                                .format("%.4f", alpha[1]) + " " + String.format("%.4f", alpha[2])
                );
            }

            // Set Divisions
            Arrays.sort(functionKeysAndLabels);
            int[] PointColors = new int[SALSAUtility.PointCount_Global];

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
                int UnsortedPosition = functionKeysAndLabels[PointIndex].getLabel();
                PointColors[UnsortedPosition] = ClusterLabel[PointIndex];
            }
            String OutputFileExtension = ManxcatCentral.config.ConversionInformation;
            OutputFileExtension = OutputFileExtension.replace(".dat", "");
            OutputFileExtension = OutputFileExtension.replace(".txt", "");
            OutputFileExtension = OutputFileExtension + "-" + pickout + ".txt";
            String labelFileDirectory = ManxcatCentral.config.ClusterDirectory; // Directory for Colors attached to
            // Points
            String OutputFileName = labelFileDirectory + File.separatorChar + OutputFileExtension;
            WritePointCluster(OutputFileName, PointColors, SALSAUtility.PointCount_Global);

            double test = 0.1 * sumabs / NumberofLines;

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
                int UnsortedPosition = functionKeysAndLabels[PointIndex].getLabel();

                if (Math.abs(functionKeysAndLabels[UnsortedPosition].getKey()) < test) {
                    PointColors[UnsortedPosition] = 0;
                }
            }

            OutputFileName = OutputFileName.replace(".txt", "NZ.txt");
            WritePointCluster(OutputFileName, PointColors, SALSAUtility.PointCount_Global);
        }
    }

    public static void UpdateManxcatMDS_Option12_FamilybyCuts() {
        String CutFileName = ManxcatCentral.config.ConversionInformation;

        if (!CutFileName.contains(":")) {
            CutFileName = ManxcatCentral.config.ControlDirectoryName + File.separatorChar + CutFileName;
        }

        String[] inputdata = new String[SALSAUtility.PointCount_Global];
        int NumberofCuts = 0;

        String FamilyLabel = "Family";
        int FamilyNumber = 1;

        if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("family2")) {
            FamilyNumber = 2;
        }

        RefObject<Integer> tempRef_NumberofCuts = new RefObject<Integer>(NumberofCuts);
        ReadString_File(CutFileName, inputdata, tempRef_NumberofCuts);
        NumberofCuts = tempRef_NumberofCuts.argValue;

        if (NumberofCuts <= 0) {
            SALSAUtility.printAndThrowRuntimeException(
                    "Too few points " + NumberofCuts + " in Cut File " + CutFileName);


        }

        int[] CutStart = new int[NumberofCuts];
        String[] CutLabel = new String[NumberofCuts];
        int CountFamilies = 0;
        int itmpold = -2;

        Pattern pattern = Pattern.compile("[\t ]");
        for (int CutPos = 0; CutPos < NumberofCuts; CutPos++) {
            String[] split = pattern.split(inputdata[CutPos]);
            int itmp = Integer.parseInt(split[0]);

            if (itmp < 0) {
                FamilyLabel = split[1];
                continue;
            }

            if (itmp <= itmpold) {
                continue;
            }
            CutStart[CountFamilies] = itmp;
            CutLabel[CountFamilies] = FamilyLabel + "-" + split[1];
            ++CountFamilies;
            itmpold = itmp;
        }

        if (FamilyNumber == 1) {
            SALSAUtility.GlobalFileProperties.FamilyName1 = FamilyLabel + " from " + CutFileName;
        }

        if (FamilyNumber == 2) {
            SALSAUtility.GlobalFileProperties.FamilyName2 = FamilyLabel + " from " + CutFileName;
        }

        int minlabel = 1;
        int CurrentFamily = 0;
        int[] CategoryIndex = new int[SALSAUtility.PointCount_Global];

        for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
            if (CurrentFamily < CountFamilies - 1) {
                if ((PointIndex + 1) >= CutStart[CurrentFamily + 1]) {
                    ++CurrentFamily;
                }
            }
            CategoryIndex[PointIndex] = minlabel + CurrentFamily;
            continue;
        }

        int[] FamilyCounts = new int[CountFamilies];

        for (int Icount = 0; Icount < CountFamilies; Icount++) {
            FamilyCounts[Icount] = 0;
        }

        for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
            ++FamilyCounts[CategoryIndex[PointIndex] - minlabel];
        }

        for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
            int Icount = CategoryIndex[PointIndex];

            if (FamilyNumber == 1) {
                SALSAUtility.GlobalPointProperties[PointIndex].family1 = Icount;
                SALSAUtility.GlobalPointProperties[PointIndex].familylabel1 = CutLabel[Icount - minlabel] + " " +
                        Icount + " out of " + CountFamilies + " Count " + FamilyCounts[Icount - minlabel];
            }

            if (FamilyNumber == 2) {
                SALSAUtility.GlobalPointProperties[PointIndex].family2 = Icount;
                SALSAUtility.GlobalPointProperties[PointIndex].familylabel2 = CutLabel[Icount - minlabel] + " " +
                        Icount + " out of " + CountFamilies + " Count " + FamilyCounts[Icount - minlabel];
            }
        }

        String cleandate = SALSAUtility.startTime.toString();
        cleandate = cleandate.replace(":", ".");
        cleandate = cleandate.replace(" ", "-");
        String OldComment = SALSAUtility.GlobalFileProperties.Comment;

        if (!OldComment.equals("")) {
            OldComment += "\n";
        }
        OldComment += "Family" + FamilyNumber + " " + FamilyLabel + " Information added from " + CutFileName + " Time" +
                " " +
                "" + cleandate;
        SALSAUtility.GlobalFileProperties.Comment = OldComment;

        //  Write new label file
        String labelfilename = ManxcatCentral.config.DataLabelsFileName;

        if (!labelfilename.contains(":")) {
            labelfilename = ManxcatCentral.config.ControlDirectoryName + File.separatorChar + labelfilename;
        }
        SALSA_Properties.WriteDataPointFile(labelfilename, ManxcatCentral.config.Write2Das3D, "colon,SameFileName ",
                SALSAUtility.GlobalFileProperties, SALSAUtility.GlobalPointProperties,
                SALSAUtility.NumberOriginalPoints);
        SALSAUtility.SALSAPrint(0,
                "Family" + FamilyNumber + " " + FamilyLabel + " Info Added to " + labelfilename + " from " +
                        CutFileName);
        ManxcatCentral.config.Comment += "\nFamily" + FamilyNumber + " " + FamilyLabel + " Info Added to " +
                labelfilename + " from " + CutFileName;
    } // End UpdateManxcatMDS_Option12_FamilybyCuts

    public static void UpdateManxcatMDS_Option12_Cluster() {
        //  Read Cluster Information with File specified in ConversionInformation
        String ClusterFileName = ManxcatCentral.config.ConversionInformation;
        int NumberOfClusterLines = 0;
        String[] CategoryLabel = new String[SALSAUtility.PointCount_Global];
        RefObject<Integer> tempRef_NumberOfClusterLines = new RefObject<Integer>(NumberOfClusterLines);
        ReadClusterLabel_File(ClusterFileName, CategoryLabel, tempRef_NumberOfClusterLines);
        NumberOfClusterLines = tempRef_NumberOfClusterLines.argValue;

        if (NumberOfClusterLines != SALSAUtility.PointCount_Global) {
            SALSAUtility.printAndThrowRuntimeException(
                    " Illegal Count " + NumberOfClusterLines + " in File " + ClusterFileName);


        }

        int[] CategoryIndex = new int[SALSAUtility.PointCount_Global];
        int minlabel = SALSAUtility.PointCount_Global;
        int maxlabel = 0;

        for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
            CategoryIndex[PointIndex] = Integer.parseInt(CategoryLabel[PointIndex]);

            if (minlabel > CategoryIndex[PointIndex]) {
                minlabel = CategoryIndex[PointIndex];
            }

            if (maxlabel < CategoryIndex[PointIndex]) {
                maxlabel = CategoryIndex[PointIndex];
            }
        }

        int TotalNumberClusters = maxlabel - minlabel + 1;
        int[] ClusterCounts = new int[TotalNumberClusters];

        for (int Icount = 0; Icount < TotalNumberClusters; Icount++) {
            ClusterCounts[Icount] = 0;
        }

        for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
            ++ClusterCounts[CategoryIndex[PointIndex] - minlabel];
        }

        for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
            int Icount = CategoryIndex[PointIndex];
            SALSAUtility.GlobalPointProperties[PointIndex].cluster = Icount;
            SALSAUtility.GlobalPointProperties[PointIndex].clusterlabel = "Cluster " + Icount + " out of " +
                    TotalNumberClusters + " Count " + ClusterCounts[Icount - minlabel];
        }
        SALSAUtility.GlobalFileProperties.ClusterName = "Clusters from " + ClusterFileName;
        SALSAUtility.GlobalFileProperties.ClusterStartIndex = minlabel;

        String cleandate = SALSAUtility.startTime.toString();
        cleandate = cleandate.replace(":", ".");
        cleandate = cleandate.replace(" ", "-");
        String OldComment = SALSAUtility.GlobalFileProperties.Comment;

        if (!OldComment.equals("")) {
            OldComment += "\n";
        }
        OldComment += "Cluster Information added from " + ClusterFileName + " Time " + cleandate;
        SALSAUtility.GlobalFileProperties.Comment = OldComment;

        //  Write new label file
        String labelfilename = ManxcatCentral.config.DataLabelsFileName;

        if (!labelfilename.contains(":") && !labelfilename.contains("$")) {
            labelfilename = ManxcatCentral.config.ControlDirectoryName + File.separatorChar + labelfilename;
        }
        SALSA_Properties.WriteDataPointFile(labelfilename, ManxcatCentral.config.Write2Das3D, "colon,SameFileName ",
                SALSAUtility.GlobalFileProperties, SALSAUtility.GlobalPointProperties,
                SALSAUtility.NumberOriginalPoints);
        SALSAUtility.SALSAPrint(0, "Cluster Info Added to " + labelfilename + " from " + ClusterFileName);
        ManxcatCentral.config.Comment += "\nCluster Info Added to " + labelfilename + " from " + ClusterFileName;
    } // End UpdateManxcatMDS_Option12_Cluster()

    public static void UpdateManxcatMDS_Option12_TraditionalDirectory() {
        String labelFileDirectory = ManxcatCentral.config.ClusterDirectory; // Input Directory
        String MDSandClusterDirectory = labelFileDirectory + File.separatorChar + ManxcatCentral.config.RunSetLabel +
                "-R" + ManxcatCentral.config.RunNumber + "-ManxcatMDS"; // Output Directory

        File MDSandClusterDirectoryFile = new File(MDSandClusterDirectory);
        if (MDSandClusterDirectoryFile.isDirectory()) {
            SALSAUtility.SALSAPrint(0, "The directory " + MDSandClusterDirectory + " exists");
        } else {
            boolean success = MDSandClusterDirectoryFile.mkdir();
            SALSAUtility.SALSAPrint(0,
                    "The directory " + MDSandClusterDirectory + " was " + (success ? ("successful at " + new Date(
                            MDSandClusterDirectoryFile.lastModified()).toString()) : "unsuccessful")
            );
        }

        String TypicalMDSFileName = ManxcatCentral.config.ReducedVectorOutputFileName;

        if (!TypicalMDSFileName.contains(":")) {
            TypicalMDSFileName = ManxcatCentral.ResultDirectoryName + "\\SIMPLE" + TypicalMDSFileName;
        }

        if (!(new java.io.File(TypicalMDSFileName)).isFile()) {
            SALSAUtility.printAndThrowRuntimeException(" File " + TypicalMDSFileName + " Does Not Exist");


        }

        int NumberMDSPoints = 0;
        String[] InitialMDSString = new String[SALSAUtility.PointCount_Global];
        int[] ColorValue = new int[SALSAUtility.PointCount_Global];

        RefObject<Integer> tempRef_NumberMDSPoints = new RefObject<Integer>(NumberMDSPoints);
        boolean tempVar = ReadMDSCluster_File(TypicalMDSFileName, InitialMDSString, ColorValue,
                tempRef_NumberMDSPoints);
        NumberMDSPoints = tempRef_NumberMDSPoints.argValue;
        if (tempVar) {
            System.out.println(TypicalMDSFileName + " Read Successfully");
        }

        File[] fileEntries = new File(labelFileDirectory).listFiles();

        for (File file : fileEntries) {
            String LabelFileName = file.toString();
            if (LabelFileName.contains(".dot")) {
                continue;
            }

            if (LabelFileName.contains("Summary")) {
                continue;
            }

            if (LabelFileName.contains("Timing")) {
                continue;
            }
            String LabelFileName1 = LabelFileName.replace(labelFileDirectory + "\\", "");
            String coreFileName = MDSandClusterDirectory + File.separatorChar + "MDSManxcat-" + LabelFileName1;

            if ((new java.io.File(coreFileName)).isFile()) {
                continue;
            }
            int NumberOfLabels = 0;
            String[] CategoryLabel = new String[SALSAUtility.PointCount_Global];
            RefObject<Integer> tempRef_NumberOfLabels = new RefObject<Integer>(NumberOfLabels);
            ReadClusterLabel_File(LabelFileName, CategoryLabel, tempRef_NumberOfLabels);
            NumberOfLabels = tempRef_NumberOfLabels.argValue;

            if (NumberOfLabels != SALSAUtility.PointCount_Global) {
                SALSAUtility.printAndThrowRuntimeException(
                        " Illegal Count " + NumberOfLabels + " in File " + LabelFileName);


            }

            int[] CategoryIndex = new int[SALSAUtility.PointCount_Global];

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++) {
                CategoryIndex[PointIndex] = Integer.parseInt(CategoryLabel[PointIndex]);
            }
            WriteColor_Cluster(coreFileName, InitialMDSString, CategoryIndex, SALSAUtility.PointCount_Global, false);
            SALSAUtility.SALSAPrint(0, "Traditional Directory MDS Info Added to " + coreFileName);
        }
    } // End UpdateManxcatMDS_Option12_TraditionalDirectory

    public static boolean ReadMDSCluster_File(String MDSClusterFileName, String[] InitialString, int[] ColorValue,
                                              RefObject<Integer> NumberofPoints) {
        NumberofPoints.argValue = 0;
        Path MDSClusterFilePath = Paths.get(MDSClusterFileName);
        if (!Files.exists(MDSClusterFilePath)) {
            SALSAUtility.printAndThrowRuntimeException("File " + MDSClusterFileName + " does not exists");
        } else {
            try (BufferedReader reader = Files.newBufferedReader(MDSClusterFilePath)) {
                // Read contents of a file, line by line, into a string
                String inputLineStr;
                Pattern pattern = Pattern.compile("[\t ]");
                while ((inputLineStr = reader.readLine()) != null) {
                    if (inputLineStr.length() < 2) {
                        continue; //replace empty line
                    }
                    inputLineStr = inputLineStr.trim();

                    inputLineStr = inputLineStr.replace("\t\t", "\t");
                    String[] split = pattern.split(inputLineStr);

                    if (split.length != Hotsun.ParameterVectorDimension + 2) {
                        System.out.println(
                                " Bad Line " + split.length + " " + NumberofPoints.argValue + " " + inputLineStr);
                    }
                    inputLineStr = split[0];

                    for (int i = 1; i < (split.length - 1); i++) {
                        inputLineStr += "\t" + split[i];
                    }
                    InitialString[NumberofPoints.argValue] = inputLineStr;

                    ColorValue[NumberofPoints.argValue] = (int) Double.parseDouble(split[split.length - 1]);
                    ++NumberofPoints.argValue;
                }
                reader.close();
            } catch (IOException e) {
                SALSAUtility.printAndThrowRuntimeException(
                        "Failed to read data from " + MDSClusterFileName + " " + e.toString());

            }

            if (NumberofPoints.argValue != SALSAUtility.PointCount_Global) {
                SALSAUtility.printAndThrowRuntimeException(
                        "Incorrect #Points in File " + NumberofPoints.argValue + " Expected " + SALSAUtility
                                .PointCount_Global);


            }
        }

        return true;
    }

    public static boolean ReadClusterLabel_File(String LabelFileName, String[] LabelArray,
                                                RefObject<Integer> NumberofLabels) {
        NumberofLabels.argValue = 0;
        int startnumber = 1;

        Path LabelFileNamePath = Paths.get(LabelFileName);
        if (!Files.exists(LabelFileNamePath)) {
            SALSAUtility.printAndThrowRuntimeException("File " + LabelFileName + " does not exists");
        } else {
            try (BufferedReader reader = Files.newBufferedReader(LabelFileNamePath)) {
                // Read contents of a file, line by line, into a string
                String inputLineStr;
                Pattern pattern = Pattern.compile("[\t ]");
                while ((inputLineStr = reader.readLine()) != null) {
                    if (inputLineStr.length() < 2) {
                        continue; //replace empty line
                    }
                    inputLineStr = inputLineStr.trim();

                    // Parse each record string
                    String[] split = pattern.split(inputLineStr);
                    int newlabelnumber = Integer.parseInt(split[0]);
                    if (NumberofLabels.argValue == 0) {
                        startnumber = newlabelnumber;
                        if ((startnumber < 0) || (startnumber > 1)) {
                            SALSAUtility.printAndThrowRuntimeException("Unexpected Start Number " + startnumber);
                        }
                    }
                    if (NumberofLabels.argValue != (newlabelnumber + startnumber)) {
                        SALSAUtility.printAndThrowRuntimeException(
                                "Unexpected Label Number " + newlabelnumber + " Expected " + NumberofLabels
                                        .argValue + " + " + startnumber
                        );
                    }
                    if (split[1].length() <= 0) {
                        SALSAUtility.printAndThrowRuntimeException(
                                "Zero length label for point " + NumberofLabels.argValue);
                    }
                    LabelArray[NumberofLabels.argValue] = split[1];
                    ++NumberofLabels.argValue;
                }
                reader.close();
            } catch (IOException e) {
                SALSAUtility.printAndThrowRuntimeException(
                        "Failed to read data from " + LabelFileName + " " + e.toString());
            }
        }
        return true;
    }


    // Write label-cluster results into a file
    public static void WriteColor_Cluster(String fname, String[] labels, int[] ColorValues, int dataPoints,
                                          boolean append) {
        if (!Strings.isNullOrEmpty(fname)) {
            try (PrintWriter writer = new PrintWriter(
                    Files.newBufferedWriter(Paths.get(fname), Charset.defaultCharset()), true)) {
                for (int i = 0; i < dataPoints; i++) {
                    String stripped = labels[i].trim();
                    writer.println(String.format(stripped + "\t" + ColorValues[i] + ".0000000000"));
                }
            } catch (IOException e) {
                SALSAUtility.printAndThrowRuntimeException("Failed writing data on " + fname + " " + e);
            }
        }
    }

    // Write Point-cluster results into a file
    public static void WritePointCluster(String fname, int[] ColorValues, int dataPoints) {
        if (!Strings.isNullOrEmpty(fname)) {
            try (PrintWriter writer = new PrintWriter(
                    Files.newBufferedWriter(Paths.get(fname), Charset.defaultCharset()), true)) {
                for (int i = 0; i < dataPoints; i++) {
                    writer.println(String.format(i + 1 + " " + ColorValues[i]));
                }
            } catch (IOException e) {
                SALSAUtility.printAndThrowRuntimeException("Failed writing data on " + fname + " " + e);
            }
        }
    }

    public static boolean ReadString_File(String fname, String[] LabelArray, RefObject<Integer> NumberofLines) {
        NumberofLines.argValue = 0;
        Path filePath = Paths.get(fname);
        if (!Files.exists(filePath)) {
            SALSAUtility.printAndThrowRuntimeException("File " + fname + " does not exists.");
        }

        try (BufferedReader reader = Files.newBufferedReader(filePath, Charset.defaultCharset())) {
            String inputLineStr;
            while ((inputLineStr = reader.readLine()) != null) {
                if (inputLineStr.length() < 2) {
                    continue; //replace empty line
                }

                LabelArray[NumberofLines.argValue] = inputLineStr.trim();
                ++NumberofLines.argValue;
            }
        } catch (IOException e) {
            SALSAUtility.printAndThrowRuntimeException("Failed to read data from " + fname + " " + e.toString());
        }
        return true;
    }
}