package salsa.configuration.sections;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

public class MDSasChisqSection {
    private int threadCount;

    public MDSasChisqSection(String configurationFilePath) {
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
            BaseResultDirectoryName = p.getProperty("BaseResultDirectoryName","");
            ControlDirectoryName = p.getProperty("ControlDirectoryName","");
            ClusterDirectory = p.getProperty("ClusterDirectory","");
            DistanceMatrixFile = p.getProperty("DistanceMatrixFile","");
            DataLabelsFileName = p.getProperty("DataLabelsFileName","");
            ReducedVectorOutputFileName = p.getProperty("ReducedVectorOutputFileName","");
            ResultDirectoryExtension = p.getProperty("ResultDirectoryExtension","");
            TimingOutputFileName = p.getProperty("TimingOutputFileName","");
            SummaryOutputFileName = p.getProperty("SummaryOutputFileName","");
            InitializationFileName = p.getProperty("InitializationFileName","");
            WeightingFileName = p.getProperty("WeightingFileName","");
            Selectedvariedpointfile = p.getProperty("Selectedvariedpointfile","");
            Selectedfixedpointfile = p.getProperty("Selectedfixedpointfile","");
            RotationLabelsFileName = p.getProperty("RotationLabelsFileName","");
            CoordinateWriteFrequency = Integer.parseInt(p.getProperty("CoordinateWriteFrequency","0"));
            IndexFile = p.getProperty("IndexFile","");
            DataPoints = Integer.parseInt(p.getProperty("DataPoints","0"));
            CalcFixedCrossFixed = Boolean.parseBoolean(p.getProperty("CalcFixedCrossFixed","false"));
            StoredDistanceOption = Integer.parseInt(p.getProperty("StoredDistanceOption","0"));
            DiskDistanceOption = Integer.parseInt(p.getProperty("DiskDistanceOption","0"));
            ProcessingOption = Integer.parseInt(p.getProperty("ProcessingOption","0"));
            DistanceProcessingOption = Integer.parseInt(p.getProperty("DistanceProcessingOption","0"));
            InitializationOption = Integer.parseInt(p.getProperty("InitializationOption","0"));
            WeightingOption = Integer.parseInt(p.getProperty("WeightingOption","0"));
            Write2Das3D = Boolean.parseBoolean(p.getProperty("Write2Das3D","false"));
            LocalVectorDimension = Integer.parseInt(p.getProperty("LocalVectorDimension","0"));
            Selectedvariedpoints = p.getProperty("Selectedvariedpoints","");
            VariedPointCriterion = p.getProperty("VariedPointCriterion","");
            Selectedfixedpoints = p.getProperty("Selectedfixedpoints","");
            FixedPointCriterion = p.getProperty("FixedPointCriterion","");
            ConversionOption = p.getProperty("ConversionOption","");
            ConversionInformation = p.getProperty("ConversionInformation","");
            RotationOption = Integer.parseInt(p.getProperty("RotationOption","0"));
            InitializationLoops = Integer.parseInt(p.getProperty("InitializationLoops","0"));
            Chisqnorm = Integer.parseInt(p.getProperty("Chisqnorm","0"));
            DistanceFormula = Integer.parseInt(p.getProperty("DistanceFormula","0"));
            FullSecondDerivativeOption = Integer.parseInt(p.getProperty("FullSecondDerivativeOption","0"));
            MinimumDistance = Double.parseDouble(p.getProperty("MinimumDistance","0.0"));
            FunctionErrorCalcMultiplier = Integer.parseInt(p.getProperty("FunctionErrorCalcMultiplier","0"));
            ChisqPrintConstant = Integer.parseInt(p.getProperty("ChisqPrintConstant","0"));
            Maxit = Integer.parseInt(p.getProperty("Maxit","0"));
            Nbadgo = Integer.parseInt(p.getProperty("Nbadgo","0"));
            ChisqChangePerPoint = Double.parseDouble(p.getProperty("ChisqChangePerPoint","0.0"));
            FletcherRho = Double.parseDouble(p.getProperty("FletcherRho","0.0"));
            FletcherSigma = Double.parseDouble(p.getProperty("FletcherSigma","0.0"));
            Omega = Double.parseDouble(p.getProperty("Omega","0.0"));
            OmegaOption = Integer.parseInt(p.getProperty("OmegaOption","0"));
            QHighInitialFactor = Double.parseDouble(p.getProperty("QHighInitialFactor","0.0"));
            QgoodReductionFactor = Double.parseDouble(p.getProperty("QgoodReductionFactor","0.0"));
            QLimitscalculationInterval = Integer.parseInt(p.getProperty("QLimitscalculationInterval","0"));
            Extraprecision = Double.parseDouble(p.getProperty("Extraprecision","0.0"));
            AddonforQcomputation = Integer.parseInt(p.getProperty("AddonforQcomputation","0"));
            InitialSteepestDescents = Integer.parseInt(p.getProperty("InitialSteepestDescents","0"));
            TimeCutmillisec = Integer.parseInt(p.getProperty("TimeCutmillisec","0"));
            CGResidualLimit = Double.parseDouble(p.getProperty("CGResidualLimit","0.0"));
            PowerIterationLimit = Integer.parseInt(p.getProperty("PowerIterationLimit","0"));
            Eigenvaluechange = Double.parseDouble(p.getProperty("Eigenvaluechange","0.0"));
            Eigenvectorchange = Double.parseDouble(p.getProperty("Eigenvectorchange","0.0"));
            Derivtest = Boolean.parseBoolean(p.getProperty("Derivtest","false"));
            RunNumber = Integer.parseInt(p.getProperty("RunNumber","0"));
            RunSetLabel = p.getProperty("RunSetLabel","");
            Pattern = p.getProperty("Pattern","");
            ThreadCount = Integer.parseInt(p.getProperty("ThreadCount","0"));
            NodeCount = Integer.parseInt(p.getProperty("NodeCount","0"));
            MPIperNodeCount = Integer.parseInt(p.getProperty("MPIperNodeCount","0"));
            MPIIOStrategy = Integer.parseInt(p.getProperty("MPIIOStrategy","0"));
            HistogramBinCount = Integer.parseInt(p.getProperty("HistogramBinCount","0"));
            Extradata1 = p.getProperty("Extradata1","");
            Extradata2 = p.getProperty("Extradata2","");
            Extradata3 = p.getProperty("Extradata3","");
            Extradata4 = p.getProperty("Extradata4","");
            ExtraOption1 = Integer.parseInt(p.getProperty("ExtraOption1","0"));
            DebugPrintOption = Integer.parseInt(p.getProperty("DebugPrintOption","0"));
            ConsoleDebugOutput = Boolean.parseBoolean(p.getProperty("ConsoleDebugOutput","false"));
            Comment = p.getProperty("Comment","");
            UndefinedDistanceValue = Double.parseDouble(p.getProperty("UndefinedDistanceValue","0.0"));
            DistanceWeightsCuts = p.getProperty("DistanceWeightsCuts","");
            DistanceCut = Float.parseFloat(p.getProperty("DistanceCut","0.0f"));
            LinkCut = Integer.parseInt(p.getProperty("LinkCut","0"));
            TransformMethod = Integer.parseInt(p.getProperty("TransformMethod","0"));
            TransformParameter = Float.parseFloat(p.getProperty("TransformParameter","0.0f"));
            ManxcatRunDescription = p.getProperty("ManxcatRunDescription","");
            ManxcatRunName = p.getProperty("ManxcatRunName","");
            SelectedClusters = p.getProperty("SelectedClusters","");
            ClusterFile = p.getProperty("ClusterFile","");
            Pcutf = Double.parseDouble(p.getProperty("Pcutf","0.0"));
            Alpha = Double.parseDouble(p.getProperty("Alpha","0.0"));
            Yres = Integer.parseInt(p.getProperty("Yres","0"));
            Xres = Integer.parseInt(p.getProperty("Xres","0"));
            YmaxBound = Double.parseDouble(p.getProperty("YmaxBound","0.0"));
            XmaxBound = Double.parseDouble(p.getProperty("XmaxBound","0.0"));
            Normalize = Boolean.parseBoolean(p.getProperty("Normalize","false"));
            ServerUrlPrefix = p.getProperty("ServerUrlPrefix","");
            dataTypeSize = Integer.parseInt(p.getProperty("DataTypeSize", "2"));
            isBigEndian = Boolean.parseBoolean(p.getProperty("IsBigEndian", "fasle"));
        }catch (IOException e){
            e.printStackTrace();
        }
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
}
