package salsa.configuration.sections;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Properties;

public class MDSasChisqSection {
    private int threadCount;
    private Properties p;

    public MDSasChisqSection(String configurationFilePath) {
        p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
            BaseResultDirectoryName = getProperty(p,"BaseResultDirectoryName","");
            ControlDirectoryName = getProperty(p,"ControlDirectoryName","");
            ClusterDirectory = getProperty(p,"ClusterDirectory","");
            DistanceMatrixFile = getProperty(p,"DistanceMatrixFile","");
            DataLabelsFileName = getProperty(p,"DataLabelsFileName","");
            ReducedVectorOutputFileName = getProperty(p,"ReducedVectorOutputFileName","");
            ResultDirectoryExtension = getProperty(p,"ResultDirectoryExtension","");
            TimingOutputFileName = getProperty(p,"TimingOutputFileName","");
            SummaryOutputFileName = getProperty(p,"SummaryOutputFileName","");
            InitializationFileName = getProperty(p,"InitializationFileName","");
            WeightingFileName = getProperty(p,"WeightingFileName","");
            Selectedvariedpointfile = getProperty(p,"Selectedvariedpointfile","");
            Selectedfixedpointfile = getProperty(p,"Selectedfixedpointfile","");
            RotationLabelsFileName = getProperty(p,"RotationLabelsFileName","");
            finalRotationFileName = getProperty(p,"finalRotationFileName");
            finalRotationPointCount = Integer.parseInt(getProperty(p,"finalRotationPointCount"));
            CoordinateWriteFrequency = Integer.parseInt(getProperty(p,"CoordinateWriteFrequency","0"));
            IndexFile = getProperty(p,"IndexFile","");
            DataPoints = Integer.parseInt(getProperty(p,"DataPoints","0"));
            CalcFixedCrossFixed = Boolean.parseBoolean(getProperty(p,"CalcFixedCrossFixed","false"));
            StoredDistanceOption = Integer.parseInt(getProperty(p,"StoredDistanceOption","0"));
            DiskDistanceOption = Integer.parseInt(getProperty(p,"DiskDistanceOption","0"));
            ProcessingOption = Integer.parseInt(getProperty(p,"ProcessingOption","0"));
            DistanceProcessingOption = Integer.parseInt(getProperty(p,"DistanceProcessingOption","0"));
            InitializationOption = Integer.parseInt(getProperty(p,"InitializationOption","0"));
            WeightingOption = Integer.parseInt(getProperty(p,"WeightingOption","0"));
            Write2Das3D = Boolean.parseBoolean(getProperty(p,"Write2Das3D","false"));
            LocalVectorDimension = Integer.parseInt(getProperty(p,"LocalVectorDimension","0"));
            Selectedvariedpoints = getProperty(p,"Selectedvariedpoints","");
            VariedPointCriterion = getProperty(p,"VariedPointCriterion","");
            Selectedfixedpoints = getProperty(p,"Selectedfixedpoints","");
            FixedPointCriterion = getProperty(p,"FixedPointCriterion","");
            ConversionOption = getProperty(p,"ConversionOption","");
            ConversionInformation = getProperty(p,"ConversionInformation","");
            RotationOption = Integer.parseInt(getProperty(p,"RotationOption","0"));
            InitializationLoops = Integer.parseInt(getProperty(p,"InitializationLoops","0"));
            Chisqnorm = Integer.parseInt(getProperty(p,"Chisqnorm","0"));
            DistanceFormula = Integer.parseInt(getProperty(p,"DistanceFormula","0"));
            FullSecondDerivativeOption = Integer.parseInt(getProperty(p,"FullSecondDerivativeOption","0"));
            MinimumDistance = Double.parseDouble(getProperty(p,"MinimumDistance","0.0"));
            FunctionErrorCalcMultiplier = Integer.parseInt(getProperty(p,"FunctionErrorCalcMultiplier","0"));
            ChisqPrintConstant = Integer.parseInt(getProperty(p,"ChisqPrintConstant","0"));
            Maxit = Integer.parseInt(getProperty(p,"Maxit","0"));
            Nbadgo = Integer.parseInt(getProperty(p,"Nbadgo","0"));
            ChisqChangePerPoint = Double.parseDouble(getProperty(p,"ChisqChangePerPoint","0.0"));
            FletcherRho = Double.parseDouble(getProperty(p,"FletcherRho","0.0"));
            FletcherSigma = Double.parseDouble(getProperty(p,"FletcherSigma","0.0"));
            Omega = Double.parseDouble(getProperty(p,"Omega","0.0"));
            OmegaOption = Integer.parseInt(getProperty(p,"OmegaOption","0"));
            QHighInitialFactor = Double.parseDouble(getProperty(p,"QHighInitialFactor","0.0"));
            QgoodReductionFactor = Double.parseDouble(getProperty(p,"QgoodReductionFactor","0.0"));
            QLimitscalculationInterval = Integer.parseInt(getProperty(p,"QLimitscalculationInterval","0"));
            Extraprecision = Double.parseDouble(getProperty(p,"Extraprecision","0.0"));
            AddonforQcomputation = Integer.parseInt(getProperty(p,"AddonforQcomputation","0"));
            InitialSteepestDescents = Integer.parseInt(getProperty(p,"InitialSteepestDescents","0"));
            TimeCutmillisec = Integer.parseInt(getProperty(p,"TimeCutmillisec","0"));
            CGResidualLimit = Double.parseDouble(getProperty(p,"CGResidualLimit","0.0"));
            PowerIterationLimit = Integer.parseInt(getProperty(p,"PowerIterationLimit","0"));
            Eigenvaluechange = Double.parseDouble(getProperty(p,"Eigenvaluechange","0.0"));
            Eigenvectorchange = Double.parseDouble(getProperty(p,"Eigenvectorchange","0.0"));
            Derivtest = Boolean.parseBoolean(getProperty(p,"Derivtest","false"));
            RunNumber = Integer.parseInt(getProperty(p,"RunNumber","0"));
            RunSetLabel = getProperty(p,"RunSetLabel","");
            Pattern = getProperty(p,"Pattern","");
            ThreadCount = Integer.parseInt(getProperty(p,"ThreadCount","0"));
            NodeCount = Integer.parseInt(getProperty(p,"NodeCount","0"));
            MPIperNodeCount = Integer.parseInt(getProperty(p,"MPIperNodeCount","0"));
            MPIIOStrategy = Integer.parseInt(getProperty(p,"MPIIOStrategy","0"));
            HistogramBinCount = Integer.parseInt(getProperty(p,"HistogramBinCount","0"));
            Extradata1 = getProperty(p,"Extradata1","");
            Extradata2 = getProperty(p,"Extradata2","");
            Extradata3 = getProperty(p,"Extradata3","");
            Extradata4 = getProperty(p,"Extradata4","");
            ExtraOption1 = Integer.parseInt(getProperty(p,"ExtraOption1","0"));
            DebugPrintOption = Integer.parseInt(getProperty(p,"DebugPrintOption","0"));
            ConsoleDebugOutput = Boolean.parseBoolean(getProperty(p,"ConsoleDebugOutput","false"));
            Comment = getProperty(p,"Comment","");
            UndefinedDistanceValue = Double.parseDouble(getProperty(p,"UndefinedDistanceValue","0.0"));
            DistanceWeightsCuts = getProperty(p,"DistanceWeightsCuts","");
            DistanceCut = Float.parseFloat(getProperty(p,"DistanceCut","0.0f"));
            LinkCut = Integer.parseInt(getProperty(p,"LinkCut","0"));
            TransformMethod = Integer.parseInt(getProperty(p,"TransformMethod","0"));
            TransformParameter = Float.parseFloat(getProperty(p,"TransformParameter","0.0f"));
            ManxcatRunDescription = getProperty(p,"ManxcatRunDescription","");
            ManxcatRunName = getProperty(p,"ManxcatRunName","");
            SelectedClusters = getProperty(p,"SelectedClusters","");
            ClusterFile = getProperty(p,"ClusterFile","");
            Pcutf = Double.parseDouble(getProperty(p,"Pcutf","0.0"));
            Alpha = Double.parseDouble(getProperty(p,"Alpha","0.0"));
            Yres = Integer.parseInt(getProperty(p,"Yres","0"));
            Xres = Integer.parseInt(getProperty(p,"Xres","0"));
            YmaxBound = Double.parseDouble(getProperty(p,"YmaxBound","0.0"));
            XmaxBound = Double.parseDouble(getProperty(p,"XmaxBound","0.0"));
            Normalize = Boolean.parseBoolean(getProperty(p,"Normalize","false"));
            ServerUrlPrefix = getProperty(p,"ServerUrlPrefix","");
            dataTypeSize = Integer.parseInt(getProperty(p,"DataTypeSize", "2"));
            isBigEndian = Boolean.parseBoolean(getProperty(p,"IsBigEndian", "false"));
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    public void saveAs(String file) throws IOException {
        p.setProperty("BaseResultDirectoryName",BaseResultDirectoryName);
        p.setProperty("ControlDirectoryName",ControlDirectoryName);
        p.setProperty("ClusterDirectory",ClusterDirectory);
        p.setProperty("DistanceMatrixFile",DistanceMatrixFile);
        p.setProperty("DataLabelsFileName",DataLabelsFileName);
        p.setProperty("ReducedVectorOutputFileName",ReducedVectorOutputFileName);
        p.setProperty("ResultDirectoryExtension",ResultDirectoryExtension);
        p.setProperty("TimingOutputFileName",TimingOutputFileName);
        p.setProperty("SummaryOutputFileName",SummaryOutputFileName);
        p.setProperty("InitializationFileName",InitializationFileName);
        p.setProperty("WeightingFileName",WeightingFileName);
        p.setProperty("Selectedvariedpointfile",Selectedvariedpointfile);
        p.setProperty("Selectedfixedpointfile",Selectedfixedpointfile);
        p.setProperty("RotationLabelsFileName",RotationLabelsFileName);
        p.setProperty("CoordinateWriteFrequency",String.valueOf(CoordinateWriteFrequency));
        p.setProperty("IndexFile",IndexFile);
        p.setProperty("DataPoints",String.valueOf(DataPoints));
        p.setProperty("CalcFixedCrossFixed",String.valueOf(CalcFixedCrossFixed));
        p.setProperty("StoredDistanceOption",String.valueOf(StoredDistanceOption));
        p.setProperty("DiskDistanceOption",String.valueOf(DiskDistanceOption));
        p.setProperty("ProcessingOption",String.valueOf(ProcessingOption));
        p.setProperty("DistanceProcessingOption",String.valueOf(DistanceProcessingOption));
        p.setProperty("InitializationOption",String.valueOf(InitializationOption));
        p.setProperty("WeightingOption",String.valueOf(WeightingOption));
        p.setProperty("Write2Das3D",String.valueOf(Write2Das3D));
        p.setProperty("LocalVectorDimension",String.valueOf(LocalVectorDimension));
        p.setProperty("Selectedvariedpoints",Selectedvariedpoints);
        p.setProperty("VariedPointCriterion",VariedPointCriterion);
        p.setProperty("Selectedfixedpoints",Selectedfixedpoints);
        p.setProperty("FixedPointCriterion",FixedPointCriterion);
        p.setProperty("ConversionOption",ConversionOption);
        p.setProperty("ConversionInformation",ConversionInformation);
        p.setProperty("RotationOption",String.valueOf(RotationOption));
        p.setProperty("InitializationLoops",String.valueOf(InitializationLoops));
        p.setProperty("Chisqnorm",String.valueOf(Chisqnorm));
        p.setProperty("DistanceFormula",String.valueOf(DistanceFormula));
        p.setProperty("FullSecondDerivativeOption",String.valueOf(FullSecondDerivativeOption));
        p.setProperty("MinimumDistance",String.valueOf(MinimumDistance));
        p.setProperty("FunctionErrorCalcMultiplier",String.valueOf(FunctionErrorCalcMultiplier));
        p.setProperty("ChisqPrintConstant",String.valueOf(ChisqPrintConstant));
        p.setProperty("Maxit",String.valueOf(Maxit));
        p.setProperty("Nbadgo",String.valueOf(Nbadgo));
        p.setProperty("ChisqChangePerPoint",String.valueOf(ChisqChangePerPoint));
        p.setProperty("FletcherRho",String.valueOf(FletcherRho));
        p.setProperty("FletcherSigma",String.valueOf(FletcherSigma));
        p.setProperty("Omega",String.valueOf(Omega));
        p.setProperty("OmegaOption",String.valueOf(OmegaOption));
        p.setProperty("QHighInitialFactor",String.valueOf(QHighInitialFactor));
        p.setProperty("QgoodReductionFactor",String.valueOf(QgoodReductionFactor));
        p.setProperty("QLimitscalculationInterval",String.valueOf(QLimitscalculationInterval));
        p.setProperty("Extraprecision",String.valueOf(Extraprecision));
        p.setProperty("AddonforQcomputation",String.valueOf(AddonforQcomputation));
        p.setProperty("InitialSteepestDescents",String.valueOf(InitialSteepestDescents));
        p.setProperty("TimeCutmillisec",String.valueOf(TimeCutmillisec));
        p.setProperty("CGResidualLimit",String.valueOf(CGResidualLimit));
        p.setProperty("PowerIterationLimit",String.valueOf(PowerIterationLimit));
        p.setProperty("Eigenvaluechange",String.valueOf(Eigenvaluechange));
        p.setProperty("Eigenvectorchange",String.valueOf(Eigenvectorchange));
        p.setProperty("Derivtest",String.valueOf(Derivtest));
        p.setProperty("RunNumber",String.valueOf(RunNumber));
        p.setProperty("RunSetLabel",RunSetLabel);
        p.setProperty("Pattern",Pattern);
        p.setProperty("ThreadCount",String.valueOf(ThreadCount));
        p.setProperty("NodeCount",String.valueOf(NodeCount));
        p.setProperty("MPIperNodeCount",String.valueOf(MPIperNodeCount));
        p.setProperty("MPIIOStrategy",String.valueOf(MPIIOStrategy));
        p.setProperty("HistogramBinCount",String.valueOf(HistogramBinCount));
        p.setProperty("Extradata1",Extradata1);
        p.setProperty("Extradata2",Extradata2);
        p.setProperty("Extradata3",Extradata3);
        p.setProperty("Extradata4",Extradata4);
        p.setProperty("ExtraOption1",String.valueOf(ExtraOption1));
        p.setProperty("DebugPrintOption",String.valueOf(DebugPrintOption));
        p.setProperty("ConsoleDebugOutput",String.valueOf(ConsoleDebugOutput));
        p.setProperty("Comment",Comment);
        p.setProperty("UndefinedDistanceValue",String.valueOf(UndefinedDistanceValue));
        p.setProperty("DistanceWeightsCuts",DistanceWeightsCuts);
        p.setProperty("DistanceCut",String.valueOf(DistanceCut));
        p.setProperty("LinkCut",String.valueOf(LinkCut));
        p.setProperty("TransformMethod",String.valueOf(TransformMethod));
        p.setProperty("TransformParameter",String.valueOf(TransformParameter));
        p.setProperty("ManxcatRunDescription",ManxcatRunDescription);
        p.setProperty("ManxcatRunName",ManxcatRunName);
        p.setProperty("SelectedClusters",SelectedClusters);
        p.setProperty("ClusterFile",ClusterFile);
        p.setProperty("Pcutf",String.valueOf(Pcutf));
        p.setProperty("Alpha",String.valueOf(Alpha));
        p.setProperty("Yres",String.valueOf(Yres));
        p.setProperty("Xres",String.valueOf(Xres));
        p.setProperty("YmaxBound",String.valueOf(YmaxBound));
        p.setProperty("XmaxBound",String.valueOf(XmaxBound));
        p.setProperty("Normalize",String.valueOf(Normalize));
        p.setProperty("ServerUrlPrefix",ServerUrlPrefix);
        p.setProperty("DataTypeSize",String.valueOf(dataTypeSize));
        p.setProperty("IsBigEndian",String.valueOf(isBigEndian));

        p.store(new FileOutputStream(file), "");
    }

    private static String getProperty(Properties p, String name, String def) {
        String val = System.getProperty(name);
        if (val == null) {
            if (def != null) {
                val = p.getProperty(name, def);
            } else {
                val = p.getProperty(name);
            }
        }
        return val;
    }

    private static String getProperty(Properties p, String name) {
        return getProperty(p, name, null);
    }

    public String BaseResultDirectoryName;
    public String ControlDirectoryName;
    public String ClusterDirectory;
    public String DistanceMatrixFile;
    public String DataLabelsFileName;
    public String ReducedVectorOutputFileName;
    public String ResultDirectoryExtension;
    public String TimingOutputFileName;
    public String SummaryOutputFileName;
    public String InitializationFileName;
    public String WeightingFileName;
    public String Selectedvariedpointfile;
    public String Selectedfixedpointfile;
    public String RotationLabelsFileName;
    public String finalRotationFileName;
    public int finalRotationPointCount;
    public int CoordinateWriteFrequency;
    public String IndexFile;
    public int DataPoints;
    public boolean CalcFixedCrossFixed;
    public int StoredDistanceOption;
    public int DiskDistanceOption;
    public int ProcessingOption;
    public int DistanceProcessingOption;
    public int InitializationOption;
    public int WeightingOption;
    public boolean Write2Das3D;
    public int LocalVectorDimension;
    public String Selectedvariedpoints;
    public String VariedPointCriterion;
    public String Selectedfixedpoints;
    public String FixedPointCriterion;
    public String ConversionOption;
    public String ConversionInformation;
    public int RotationOption;
    public int InitializationLoops;
    public int Chisqnorm;
    public int DistanceFormula;
    public int FullSecondDerivativeOption;
    public double MinimumDistance;
    public int FunctionErrorCalcMultiplier;
    public int ChisqPrintConstant;
    public int Maxit;
    public int Nbadgo;
    public double ChisqChangePerPoint;
    public double FletcherRho;
    public double FletcherSigma;
    public double Omega;
    public int OmegaOption;
    public double QHighInitialFactor;
    public double QgoodReductionFactor;
    public int QLimitscalculationInterval;
    public double Extraprecision;
    public int AddonforQcomputation;
    public int InitialSteepestDescents;
    public int TimeCutmillisec;
    public double CGResidualLimit;
    public int PowerIterationLimit;
    public double Eigenvaluechange;
    public double Eigenvectorchange;
    public boolean Derivtest;
    public int RunNumber;
    public String RunSetLabel;
    public String Pattern;
    public int ThreadCount;
    public int NodeCount;
    public int MPIperNodeCount;
    public int MPIIOStrategy;
    public int HistogramBinCount;
    public String Extradata1;
    public String Extradata2;
    public String Extradata3;
    public String Extradata4;
    public int ExtraOption1;
    public int DebugPrintOption;
    public boolean ConsoleDebugOutput;
    public String Comment;
    public double UndefinedDistanceValue;
    public String DistanceWeightsCuts;
    public float DistanceCut;
    public int LinkCut;
    public int TransformMethod;
    public float TransformParameter;
    public String ManxcatRunDescription;
    public String ManxcatRunName;
    public String SelectedClusters;
    public String ClusterFile;
    public double Pcutf;
    public double Alpha;
    public int Yres;
    public int Xres;
    public double YmaxBound;
    public double XmaxBound;
    public boolean Normalize;
    public String ServerUrlPrefix;
    public int dataTypeSize = 2; // 2 for short
    public boolean isBigEndian = false; // true for Java style binary data and false for C# style binary data


    public int getNodeCount() {
        return NodeCount;
    }

    public void setNodeCount(int nodeCount) {
        this.NodeCount = nodeCount;
    }

    public int getThreadCount() {
        return threadCount;
    }

    public void setThreadCount(int threadCount) {
        this.threadCount = threadCount;
    }

    public int getDataTypeSize() {
        return dataTypeSize;
    }

    public boolean isBigEndian() {
        return isBigEndian;
    }
}
