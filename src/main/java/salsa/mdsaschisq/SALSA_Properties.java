package salsa.mdsaschisq;

import com.google.common.base.Strings;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.regex.Pattern;

public class SALSA_Properties {
    /**
     * Read data point file
     *
     * @param mdsClusterFileName
     * @param FileType           FileType = 0 simple integers and positions:  PointNumber x y z LABEL or Group Number (LABEL if file name contains label)
     *                           FileType = 1 Just integers -- Point Number and Group Number: PointNumber  Group Number
     *                           FileType = 2 Integer and Label: Point LABEL (LABEL if file name contains label)
     *                           FileType = 3 Pure colon style -- with positions
     *                           FileType = 4 Pure Colon Style -- no positions
     *                           FileType = 5 Pure Colon Style -- Labels
     *                           FileType = 6 Hybrid Style -- with positions
     *                           FileType = 7 Hybrid Style -- no positions
     *                           FileType = 8 Hybrid Style -- Labels
     * @param FileProperties
     * @param DataPoints
     * @param NumberofPoints
     * @return
     */
    public static boolean ReadDataPointFile(String mdsClusterFileName, tangible.RefObject<Integer> FileType,
                                            SALSAFileProperties FileProperties,
                                            tangible.RefObject<SALSADataPointProperties[]> DataPoints,
                                            tangible.RefObject<Integer> NumberofPoints) {

        boolean positionSet = false;
        boolean colonsInput = false;
        boolean nonColonsInput = false;
        boolean labelsFile = false;
        String lowerFileName = mdsClusterFileName.toLowerCase();
        if (lowerFileName.contains("label")) {
            labelsFile = true;
        }

        NumberofPoints.argValue = 0;

        // Check if file exists
        Path mdsClusterFilePath = Paths.get(mdsClusterFileName);
        if (!Files.exists(mdsClusterFilePath)) {
            SALSAUtility.printAndThrowRuntimeException("File " + mdsClusterFileName + " does not exists");
        }

        try (BufferedReader br = Files.newBufferedReader(mdsClusterFilePath, Charset.defaultCharset())) {
            // Read contents of a file
            String inputLineStr;
            int newLabelNumber = -1;
            int localStart = FileProperties.LocalPointStartIndex;
            int originalStart = FileProperties.OriginalPointStartIndex;
            Pattern pattern = Pattern.compile("[\t ]");
            while ((inputLineStr = br.readLine()) != null) {
                if (inputLineStr.length() < 2) {
                    continue; //replace empty line
                }

                inputLineStr = tangible.DotNetToJavaStringHelper.trim(inputLineStr, ' ', '\t');
                try {
                    // Parse each record string
                    inputLineStr = inputLineStr.replace("\t\t", "\t");

                    if (inputLineStr.contains("FileProperties")) { // Not a Data Point
                        String endofInput = inputLineStr.replace("FileProperties", "");
                        ReadFileProperties(FileProperties, endofInput);
                        colonsInput = true;
                        localStart = FileProperties.LocalPointStartIndex;
                        originalStart = FileProperties.OriginalPointStartIndex;
                    } else { // Must be a Data Point
                        ArrayInitializer(DataPoints, SALSAUtility.NumberOriginalPoints,
                                FileProperties.NumberPointsinFile);
                        DataPoints.argValue[NumberofPoints.argValue] = new SALSADataPointProperties();
                        boolean incrementCount = false;
                        int PointDataStarts = inputLineStr.indexOf("PointProperties");
                        if (PointDataStarts >= 0) { // Some Colon Information
                            String EndofInput = inputLineStr.substring(PointDataStarts);
                            EndofInput = EndofInput.replace("PointProperties", "");
                            ReadDataPointProperties(DataPoints.argValue[NumberofPoints.argValue], EndofInput);
                            colonsInput = true;
                            if (DataPoints.argValue[NumberofPoints.argValue].valuesset) {
                                positionSet = true;
                            }
                            incrementCount = true;
                            DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber = NumberofPoints.argValue +
                                    FileProperties.LocalPointStartIndex;
                        } //  End Processing Colon Point Information

                        if (PointDataStarts < 0) { // traditional bare line
                            PointDataStarts = inputLineStr.length();
                        }
                        if (PointDataStarts > 0) { // Process number information
                            String startofString = inputLineStr.substring(0, PointDataStarts);
                            startofString = tangible.DotNetToJavaStringHelper.trim(startofString, ' ', '\t');

                            if (startofString.length() > 0) {
                                // You come here only for traditional bare line of x,y,z coordinates.

                                String[] split = pattern.split(startofString);

                                if ((split.length != 5) && (split.length != 2) && (split.length != 4)) {
                                    SALSAUtility.printAndThrowRuntimeException(
                                            " Bad Line " + split.length + " " + NumberofPoints.argValue + " " +
                                                    inputLineStr
                                    );

                                }
                                newLabelNumber = Integer.parseInt(split[0]);

                                if ((NumberofPoints.argValue + localStart) != newLabelNumber) {
                                    SALSAUtility.printAndThrowRuntimeException(
                                            "Unexpected Label Number " + newLabelNumber + " Expected " +
                                                    NumberofPoints.argValue + " + " + localStart
                                    );

                                }
                                if (DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber >= 0) {
                                    if ((DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber - localStart)
                                            != NumberofPoints.argValue) {
                                        SALSAUtility.printAndThrowRuntimeException(
                                                "Unexpected Local Number " + DataPoints.argValue[NumberofPoints
                                                        .argValue].LocalPointNumber + " Expected " +
                                                        NumberofPoints.argValue + " + " + localStart
                                        );

                                    }
                                }
                                DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber = NumberofPoints.argValue;
                                if (DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber >= 0) {
                                    if ((DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber -
                                            originalStart) < 0) {
                                        SALSAUtility.printAndThrowRuntimeException(
                                                "Unexpected Original Number " + DataPoints
                                                        .argValue[NumberofPoints.argValue].OriginalPointNumber +
                                                        " Local Point " + NumberofPoints.argValue + " Original " +
                                                        "Increment " + originalStart
                                        );

                                    }
                                    DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber -= originalStart;
                                } else {
                                    DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber = newLabelNumber;
                                }

                                if (labelsFile) {
                                    DataPoints.argValue[NumberofPoints.argValue].pointlabel = split[split.length - 1];
                                } else {
                                    DataPoints.argValue[NumberofPoints.argValue].group = Integer
                                            .parseInt(split[split.length - 1]);
                                }

                                if (split.length >= 4) {
                                    DataPoints.argValue[NumberofPoints.argValue].valuesset = true;
                                    DataPoints.argValue[NumberofPoints.argValue].x = Double.parseDouble(split[1]);
                                    DataPoints.argValue[NumberofPoints.argValue].y = Double.parseDouble(split[2]);
                                    positionSet = true;
                                }
                                if (split.length == 5) {
                                    DataPoints.argValue[NumberofPoints.argValue].z = Double.parseDouble(split[3]);
                                }
                                incrementCount = true;
                                nonColonsInput = true;
                            } // End Processing non colon Point information
                        }
                        if (incrementCount) {
                            ++NumberofPoints.argValue;
                        }
                    }
                } catch (RuntimeException e) {
                    SALSAUtility.printAndThrowRuntimeException(
                            "Failed to load data array " + inputLineStr + " " + " " + NumberofPoints.argValue + " " +
                                    newLabelNumber + " " + e
                                    .getMessage()
                    );
                }
            }

            FileType.argValue = 1;
            if (positionSet) {
                FileType.argValue = 0;
            }
            if (labelsFile) {
                FileType.argValue = 2;
            }
            if (colonsInput) {
                if (nonColonsInput) {
                    FileType.argValue += 6;
                } else {
                    FileType.argValue += 3;
                }
            }
            if (FileProperties.NumberOriginalPoints == 0) {
                FileProperties.NumberOriginalPoints = NumberofPoints.argValue;
            }

            if (FileProperties.NumberPointsinFile == 0) {
                FileProperties.NumberPointsinFile = NumberofPoints.argValue;
            }

            if (FileProperties.NumberPointsinFile != NumberofPoints.argValue) {
                SALSAUtility.printAndThrowRuntimeException(
                        "Unexpected Number of Points in File " + NumberofPoints.argValue + " Read but Expected "
                                + FileProperties.NumberPointsinFile
                );

            }
        } catch (IOException e) {
            SALSAUtility.printAndThrowRuntimeException("Failed to read data from " + mdsClusterFileName + " " + e);
        }
        return true;
    }

    public static void ArrayInitializer(tangible.RefObject<SALSADataPointProperties[]> DataArray, int sizemax,
                                        int sizereadin) {
        if (DataArray.argValue != null) {
            if (DataArray.argValue.length < sizereadin) {
                SALSAUtility.printAndThrowRuntimeException(
                        " Data Array too small for file Length " + DataArray.argValue.length + " Needs " +
                                sizereadin
                );

            }
            return;
        }
        int size = sizereadin;
        if (size == 0) {
            size = sizemax;
        }
        DataArray.argValue = new SALSADataPointProperties[size];

    }

    // Write label-cluster results into a file
    public static void WriteDataPointFile(String CoreOutputFileName, boolean write2Das3D, String OutputStyles,
                                          SALSAFileProperties FileProperties, SALSADataPointProperties[] DataPoints,
                                          int NumberofDataPoints) {
        boolean[] DothisOutputStyle = new boolean[5];
        for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++) {
            DothisOutputStyle[StyleIndex] = OutputStyles.contains("all");
        }
        if (OutputStyles.contains("colon")) {
            DothisOutputStyle[0] = true;
        }
        if (OutputStyles.contains("family1")) {
            DothisOutputStyle[1] = true;
        }
        if (OutputStyles.contains("family2")) {
            DothisOutputStyle[2] = true;
        }
        if (OutputStyles.contains("cluster")) {
            DothisOutputStyle[3] = true;
        }
        if (OutputStyles.contains("group")) {
            DothisOutputStyle[4] = true;
        }

        boolean setgroup = false;
        boolean OutputValues = false;

        SALSADataPointProperties FirstOne = DataPoints[0];
        if (FirstOne.family1 < 0) {
            DothisOutputStyle[1] = false;
        }
        if (FirstOne.family2 < 0) {
            DothisOutputStyle[2] = false;
        }
        if (FirstOne.cluster < 0) {
            DothisOutputStyle[3] = false;
        }
        if (FirstOne.group < 0) {
            DothisOutputStyle[4] = false;
        }
        if ((FirstOne.family1 < 0) && (FirstOne.family2 < 0) && (FirstOne.cluster < 0) && (FirstOne.group < 0)) {
            DothisOutputStyle[4] = true;
            setgroup = true;
        }
        if (FirstOne.valuesset) {
            OutputValues = true;
        }
        int LocalVectorDimension = FileProperties.LocalVectorDimension;
        int LocalPointIncrement = FileProperties.LocalPointStartIndex;

        int numberoffiles = 0;
        for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++) {
            if (DothisOutputStyle[StyleIndex]) {
                ++numberoffiles;
            }
        }
        if (numberoffiles == 0) {
            SALSAUtility.printAndThrowRuntimeException("No files output for core name " + CoreOutputFileName);
            return;
        }
        if (numberoffiles > 1) {
            if (OutputStyles.contains("SameFileName")) {
                SALSAUtility.printAndThrowRuntimeException(
                        "Attempt to generate multiple outputs to same file " + CoreOutputFileName);

            }
        }
        for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++) {
            if (!DothisOutputStyle[StyleIndex]) {
                continue;
            }
            String OutputFileName = "";
            if (!OutputStyles.contains("SameFileName")) {
                if (StyleIndex == 0) {
                    OutputFileName = CoreOutputFileName.replace(".txt", "Colon.txt");
                }
                if (StyleIndex == 1) {
                    OutputFileName = CoreOutputFileName.replace(".txt", "Family1.txt");
                }
                if (StyleIndex == 2) {
                    OutputFileName = CoreOutputFileName.replace(".txt", "Family2.txt");
                }
                if (StyleIndex == 3) {
                    OutputFileName = CoreOutputFileName.replace(".txt", "Cluster.txt");
                }
                if (StyleIndex == 4) {
                    OutputFileName = CoreOutputFileName.replace(".txt", "Group.txt");
                }
            } else {
                OutputFileName = CoreOutputFileName;
            }


            if (Strings.isNullOrEmpty(OutputFileName)) continue;

            try (PrintWriter sw = new PrintWriter(Files.newBufferedWriter(Paths.get(OutputFileName),
                    Charset.defaultCharset()))) {

                if (StyleIndex == 0) {
                    WriteFileProperties(FileProperties, sw); // Write Header of a Colon File
                }

                for (int GlobalDataPoint = 0; GlobalDataPoint < NumberofDataPoints; GlobalDataPoint++) {
                    SALSADataPointProperties ThesePointProperties = DataPoints[GlobalDataPoint];
                    if ((LocalVectorDimension == 2) && write2Das3D && OutputValues) {
                        ThesePointProperties.z = 0.0;
                        ThesePointProperties.zerr = 0.0;
                    }
                    if (StyleIndex == 0) {
                        String OutputLine = "";
                        tangible.RefObject<String> tempRef_OutputLine = new tangible.RefObject<String>(OutputLine);
                        AppendDataPointProperties(ThesePointProperties, tempRef_OutputLine);
                        OutputLine = tempRef_OutputLine.argValue;
                        sw.println(OutputLine);
                        continue;
                    }
                    int IntegerIndex = 0;
                    if (StyleIndex == 1) {
                        IntegerIndex = ThesePointProperties.family1;
                    }
                    if (StyleIndex == 2) {
                        IntegerIndex = ThesePointProperties.family2;
                    }
                    if (StyleIndex == 3) {
                        IntegerIndex = ThesePointProperties.cluster;
                    }
                    if (StyleIndex == 4) {
                        IntegerIndex = ThesePointProperties.group;
                        if (setgroup) {
                            IntegerIndex = 1;
                        }
                    }
                    String Coordinates = "";
                    if (OutputValues) {
                        Coordinates = String.format("%.4E", ThesePointProperties.x) + "\t" + String.format("%.4E",
                                ThesePointProperties.y) + "\t";
                        if ((LocalVectorDimension == 3) || write2Das3D) {
                            Coordinates += String.format("%.4E", ThesePointProperties.z) + "\t";
                        }
                    }
                    sw.println((GlobalDataPoint + LocalPointIncrement) + "\t" + Coordinates + IntegerIndex);
                }
                sw.close();
            } catch (IOException e) {
                SALSAUtility.printAndThrowRuntimeException(" Failed to Write Properties File " + CoreOutputFileName);

            }
        }
    }

    public static void ReadFileProperties(SALSAFileProperties FileProps, String inputLine) {
        inputLine = tangible.DotNetToJavaStringHelper.trim(inputLine, ' ', '\t');
        Pattern pattern = Pattern.compile("[:]");
        String[] colonSplits = pattern.split(inputLine);
        if (colonSplits.length < 2) {
            System.out.println("No deliminator in Line " + inputLine);
            return;
        }
        colonSplits[0] = colonSplits[0].trim();
        colonSplits[1] = colonSplits[1].trim();
        if (colonSplits[0].equals("LocalVectorDimension")) {
            FileProps.LocalVectorDimension = Integer.parseInt(colonSplits[1]);
        }
        if (colonSplits[0].equals("ClusterStartIndex")) {
            FileProps.ClusterStartIndex = Integer.parseInt(colonSplits[1]);
        }
        if (colonSplits[0].equals("OriginalPointStartIndex")) {
            FileProps.OriginalPointStartIndex = Integer.parseInt(colonSplits[1]);
        }
        if (colonSplits[0].equals("LocalPointStartIndex")) {
            FileProps.LocalPointStartIndex = Integer.parseInt(colonSplits[1]);
        }
        if (colonSplits[0].equals("FileGenerationType")) {
            FileProps.FileGenerationType = Integer.parseInt(colonSplits[1]);
        }
        if (colonSplits[0].equals("FamilyName1")) {
            FileProps.FamilyName1 = colonSplits[1];
        }
        if (colonSplits[0].equals("FamilyName2")) {
            FileProps.FamilyName2 = colonSplits[1];
        }
        if (colonSplits[0].equals("GroupName")) {
            FileProps.GroupName = colonSplits[1];
        }
        if (colonSplits[0].equals("ClusterName")) {
            FileProps.ClusterName = colonSplits[1];
        }
        if (colonSplits[0].equals("NumberOriginalPoints")) {
            FileProps.NumberOriginalPoints = Integer.parseInt(colonSplits[1]);
        }
        if (colonSplits[0].equals("NumberPointsinFile")) {
            FileProps.NumberPointsinFile = Integer.parseInt(colonSplits[1]);
        }
        if (colonSplits[0].equals("RotationParameters")) {
            FileProps.RotationParameters = colonSplits[1];
        }
        if (colonSplits[0].equals("NumberRotationParameters")) {
            FileProps.NumberRotationParameters = Integer.parseInt(colonSplits[1]);
        }
        if (colonSplits[0].equals("Comment")) {
            if (!FileProps.Comment.equals("")) {
                FileProps.Comment += "\n";
            }
            FileProps.Comment += colonSplits[1];
        }
    }

    public static void WriteFileProperties(SALSAFileProperties FileProps, PrintWriter writer) {
        writer.println(String.format("FileProperties\tLocalVectorDimension:%1$s", FileProps.LocalVectorDimension));
        writer.println(String.format("FileProperties\tClusterStartIndex:%1$s", FileProps.ClusterStartIndex));
        writer.println(
                String.format("FileProperties\tOriginalPointStartIndex:%1$s", FileProps.OriginalPointStartIndex));
        writer.println(String.format("FileProperties\tLocalPointStartIndex:%1$s", FileProps.LocalPointStartIndex));
        writer.println(String.format("FileProperties\tFileGenerationType:%1$s", FileProps.FileGenerationType));
        writer.println(String.format("FileProperties\tFamilyName1:%1$s", FileProps.FamilyName1));
        writer.println(String.format("FileProperties\tFamilyName2:%1$s", FileProps.FamilyName2));
        writer.println(String.format("FileProperties\tGroupName:%1$s", FileProps.GroupName));
        writer.println(String.format("FileProperties\tClusterName:%1$s", FileProps.ClusterName));
        writer.println(String.format("FileProperties\tNumberOriginalPoints:%1$s", FileProps.NumberOriginalPoints));
        writer.println(String.format("FileProperties\tNumberPointsinFile:%1$s", FileProps.NumberPointsinFile));
        writer.println(String.format("FileProperties\tRotationParameters:%1$s", FileProps.RotationParameters));
        writer.println(
                String.format("FileProperties\tNumberRotationParameters:%1$s", FileProps.NumberRotationParameters));

        String FileComments = FileProps.Comment;
        if (!FileComments.equals("")) {
            String[] split = FileComments.split("[\\n]", -1);
            for (String aSplit : split) {
                if (!aSplit.equals("")) {
                    writer.println(String.format("FileProperties\tComment:%1$s", aSplit));
                }
            }
        }
    } // End WriteFileProperties

    public static void ReadDataPointProperties(SALSADataPointProperties DataPointProps, String inputLine) {
        inputLine = tangible.DotNetToJavaStringHelper.trim(inputLine, ' ', '\t');
        Pattern pattern = Pattern.compile("[\t ]");
        Pattern colonPattern = Pattern.compile("[:]");
        String[] splits = colonPattern.split(inputLine);
        DataPointProps.valuesset = false;
        DataPointProps.errorsset = false;

        for (int itemBackwards = 0; itemBackwards < splits.length; itemBackwards++) {
            int item = splits.length - itemBackwards - 1;
            if (splits[item].equals("PointProperties")) {
                continue;
            }
            String[] colonSplits = colonPattern.split(splits[item]);
            if (colonSplits[0].equals("PointProperties")) {
                continue;
            }
            if (colonSplits.length < 2) {
                if (item == 0) {
                    System.out.println("No deliminator in item " + splits[item] + " in Line " + inputLine);
                } else {
                    splits[item - 1] += " " + splits[item];
                    continue;
                }
            }
            colonSplits[0] = colonSplits[0].trim();
            colonSplits[1] = colonSplits[1].trim();
            if (colonSplits[0].equals("x")) {
                DataPointProps.x = Double.parseDouble(colonSplits[1]);
                DataPointProps.valuesset = true;
            }
            if (colonSplits[0].equals("y")) {
                DataPointProps.y = Double.parseDouble(colonSplits[1]);
            }
            if (colonSplits[0].equals("z")) {
                DataPointProps.z = Double.parseDouble(colonSplits[1]);
            }
            if (colonSplits[0].equals("xerr")) {
                DataPointProps.xerr = Double.parseDouble(colonSplits[1]);
                DataPointProps.errorsset = true;
            }
            if (colonSplits[0].equals("source")) {
                DataPointProps.source = colonSplits[1];
            }
            if (colonSplits[0].equals("yerr")) {
                DataPointProps.yerr = Double.parseDouble(colonSplits[1]);
            }
            if (colonSplits[0].equals("zerr")) {
                DataPointProps.zerr = Double.parseDouble(colonSplits[1]);
            }
            if (colonSplits[0].equals("family1")) {
                DataPointProps.family1 = Integer.parseInt(colonSplits[1]);
            }
            if (colonSplits[0].equals("familylabel1")) {
                DataPointProps.familylabel1 = colonSplits[1];
            }
            if (colonSplits[0].equals("family2")) {
                DataPointProps.family2 = Integer.parseInt(colonSplits[1]);
            }
            if (colonSplits[0].equals("familylabel2")) {
                DataPointProps.familylabel2 = colonSplits[1];
            }
            if (colonSplits[0].equals("cluster")) {
                DataPointProps.cluster = Integer.parseInt(colonSplits[1]);
            }
            if (colonSplits[0].equals("clusterlabel")) {
                DataPointProps.clusterlabel = colonSplits[1];
            }
            if (colonSplits[0].equals("group")) {
                DataPointProps.group = Integer.parseInt(colonSplits[1]);
            }
            if (colonSplits[0].equals("grouplabel")) {
                DataPointProps.grouplabel = colonSplits[1];
            }
            if (colonSplits[0].equals("pointlabel")) {
                DataPointProps.pointlabel = colonSplits[1];
            }
            if (colonSplits[0].equals("FixedorVaried")) {
                DataPointProps.FixedorVaried = Integer.parseInt(colonSplits[1]);
            }
            if (colonSplits[0].equals("PointType")) {
                DataPointProps.PointType = Integer.parseInt(colonSplits[1]);
            }
            if (colonSplits[0].equals("LocalPointNumber")) {
                DataPointProps.LocalPointNumber = Integer.parseInt(colonSplits[1]);
            }
            if (colonSplits[0].equals("OriginalPointNumber")) {
                DataPointProps.OriginalPointNumber = Integer.parseInt(colonSplits[1]);
            }
        }
    }

    public static void AppendDataPointProperties(SALSADataPointProperties DataPointProps,
                                                 tangible.RefObject<String> InputLine) { //Calling Routine must supply any needed delimiter before this is appended
        InputLine.argValue += "PointProperties";
        if (DataPointProps.valuesset) { // x y z are set

            InputLine.argValue += "\tvaluesset:true\tsource:" + DataPointProps.source + "\tx:" + String.format("%.4E",
                    DataPointProps.x);
            InputLine.argValue += "\ty:" + String.format("%.4E", DataPointProps.y);
            InputLine.argValue += "\tz:" + String.format("%.4E", DataPointProps.z);
        }
        if (DataPointProps.errorsset) { // x y z errors are set
            InputLine.argValue += "\terrorsset:true\txerr:" + String.format("%.4E", DataPointProps.xerr);
            InputLine.argValue += "\tyerr:" + String.format("%.4E", DataPointProps.yerr);
            InputLine.argValue += "\tzerr:" + String.format("%.4E", DataPointProps.zerr);
        }
        if (DataPointProps.family1 >= 0) {
            InputLine.argValue += "\tfamily1:" + DataPointProps.family1;
            if (DataPointProps.familylabel1.length() > 0) {
                InputLine.argValue += "\tfamilylabel1:" + DataPointProps.familylabel1;
            }
        }
        if (DataPointProps.family2 >= 0) {
            InputLine.argValue += "\tfamily2:" + DataPointProps.family2;
            if (DataPointProps.familylabel2.length() > 0) {
                InputLine.argValue += "\tfamilylabel2:" + DataPointProps.familylabel2;
            }
        }
        if (DataPointProps.cluster >= 0) {
            InputLine.argValue += "\tcluster:" + DataPointProps.cluster;
            if (DataPointProps.clusterlabel.length() > 0) {
                InputLine.argValue += "\tclusterlabel:" + DataPointProps.clusterlabel;
            }
        }
        if (DataPointProps.group >= 0) {
            InputLine.argValue += "\tgroup:" + DataPointProps.group;
            if (DataPointProps.grouplabel.length() > 0) {
                InputLine.argValue += "\tgrouplabel:" + DataPointProps.grouplabel;
            }
        }
        if (DataPointProps.pointlabel.length() > 0) {
            InputLine.argValue += "\tpointlabel:" + DataPointProps.pointlabel;
        }
        InputLine.argValue += "\tFixedorVaried:" + DataPointProps.FixedorVaried;
        InputLine.argValue += "\tPointType:" + DataPointProps.PointType;
        if (DataPointProps.LocalPointNumber >= 0) {
            InputLine.argValue += "\tLocalPointNumber:" + DataPointProps.LocalPointNumber;
        }
        if (DataPointProps.OriginalPointNumber >= 0) {
            InputLine.argValue += "\tOriginalPointNumber:" + DataPointProps.OriginalPointNumber;
        }

    }
}