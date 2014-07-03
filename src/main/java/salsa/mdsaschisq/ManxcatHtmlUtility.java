package salsa.mdsaschisq;

import com.google.common.base.Strings;
import salsa.configuration.sections.MDSasChisqSection;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;

public class ManxcatHtmlUtility {


    /**
     * Generates a Web page linking to output files
     */
    public static void WriteHTML() {
        MDSasChisqSection config = ManxcatCentral.config;
        String directoryPath = (new java.io.File(config.ReducedVectorOutputFileName)).getParent();
        directoryPath = tangible.DotNetToJavaStringHelper.isNullOrEmpty(directoryPath) ? "" : directoryPath;

        String stylesFile = Paths.get(directoryPath, ManxcatConstants.StyleFileNameTag) + ".css";
        CopyStylesCSS(stylesFile);

        String htmlFile = Paths.get(directoryPath, ManxcatConstants.HtmlFileNameTag) + ".html";
        GenerateHtmlContent(htmlFile, config);
    }

    private static void GenerateHtmlContent(String htmlFile, MDSasChisqSection config) {
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(Paths.get(htmlFile)))) {

            try (InputStream stream = ManxcatHtmlUtility.class.getClassLoader().getResourceAsStream("template")) {
                if (stream != null) {
                    try (BufferedReader reader = new BufferedReader(new InputStreamReader(stream))) {
                        String line;
                        while ((line = reader.readLine()) != null) {
                            if (!Strings.isNullOrEmpty(line)) {
                                writer.println(
                                        !line.contains(ManxcatConstants.Tag) ? line : ProcessTagLine(line, config));
                            }
                        }
                    }
                } else {
                    SALSAUtility.printAndThrowRuntimeException(ManxcatErrorMessages.TemplateHtmlNotAvailable);
                }
            }
        } catch (IOException e) {
            SALSAUtility.printAndThrowRuntimeException(ManxcatErrorMessages.HtmlGenerationFailed);
        }
    }

    private static String ProcessTagLine(String line, MDSasChisqSection config) {
        int tagStartIdx = line.indexOf(ManxcatConstants.Tag);
        int tagEndIdx = line.indexOf(ManxcatConstants.Tag,
                tagStartIdx + ManxcatConstants.Tag.length()) + (ManxcatConstants.Tag.length() - 1);
        String prefix = line.substring(0, tagStartIdx);
        String suffix = line.substring(tagEndIdx + 1);
        String meat = line.substring(tagStartIdx + ManxcatConstants.Tag.length(),
                tagStartIdx + ManxcatConstants.Tag.length() + ((tagEndIdx + 1) - (tagStartIdx + 2 * ManxcatConstants
                        .Tag.length())));

        switch (meat) {
            case ManxcatConstants.TagManxcatName:
                return prefix + config.ManxcatRunName + suffix;
            case ManxcatConstants.TagManxcatDescription:
                return prefix + config.ManxcatRunDescription + suffix;
            case ManxcatConstants.TagManxcatConfig:
                return GenerateConfig(config);
            case ManxcatConstants.TagManxcatLinks:
                return prefix + GenerateLinks(config) + suffix;
            default:
                throw new RuntimeException(ManxcatErrorMessages.UnidentifiedTagInHTMLGeneration);

        }

    }

    private static String GenerateConfig(MDSasChisqSection config) {
        StringBuilder sb = new StringBuilder();
        sb.append("I/O").append("\n");
        sb.append("\tCoordinateWriteFrequency:      ").append(config.CoordinateWriteFrequency).append("\n");
        sb.append("\tDistanceMatrixFile:            ").append(config.DistanceMatrixFile).append("\n");
        sb.append("\n\nManxcatCore").append("\n");
        sb.append("\tAddonforQcomputation:          ").append(config.AddonforQcomputation).append("\n");
        sb.append("\tCalcFixedCrossFixed:           ").append(config.CalcFixedCrossFixed).append("\n");
        sb.append("\tCGResidualLimit:               ").append(config.CGResidualLimit).append("\n");
        sb.append("\tChisqChangePerPoint:           ").append(config.ChisqChangePerPoint).append("\n");
        sb.append("\tChisqnorm:                     ").append(config.Chisqnorm).append("\n");
        sb.append("\tChisqPrintConstant:            ").append(config.ChisqPrintConstant).append("\n");
        sb.append("\tConversionInformation:         ").append(config.ConversionInformation).append("\n");
        sb.append("\tConversionOption:              ").append(config.ConversionOption).append("\n");
        sb.append("\tDataPoints:                    ").append(config.DataPoints).append("\n");
        sb.append("\tDerivtest:                     ").append(config.Derivtest).append("\n");
        sb.append("\tDiskDistanceOption:            ").append(config.DiskDistanceOption).append("\n");
        sb.append("\tDistanceCut:                   ").append(config.DistanceCut).append("\n");
        sb.append("\tDistanceFormula:               ").append(config.DistanceFormula).append("\n");
        sb.append("\tDistanceProcessingOption:      ").append(config.DistanceProcessingOption).append("\n");
        sb.append("\tDistanceWeigthsCuts:           ").append(config.DistanceWeightsCuts).append("\n");
        sb.append("\tEigenvaluechange:              ").append(config.Eigenvaluechange).append("\n");
        sb.append("\tEigenvectorchange:             ").append(config.Eigenvectorchange).append("\n");
        sb.append("\tExtraOption1:                  ").append(config.ExtraOption1).append("\n");
        sb.append("\tExtraprecision:                ").append(config.Extraprecision).append("\n");
        sb.append("\tFixedPointCriterion:           ").append(config.FixedPointCriterion).append("\n");
        sb.append("\tFletcherRho:                   ").append(config.FletcherRho).append("\n");
        sb.append("\tFletcherSigma:                 ").append(config.FletcherSigma).append("\n");
        sb.append("\tFullSecondDerivativeOption:    ").append(config.FullSecondDerivativeOption).append("\n");
        sb.append("\tFunctionErrorCalcMultiplier:   ").append(config.FunctionErrorCalcMultiplier).append("\n");
        sb.append("\tHistogramBinCount              ").append(config.HistogramBinCount).append("\n");
        sb.append("\tInitializationLoops:           ").append(config.InitializationLoops).append("\n");
        sb.append("\tInitializationOption:          ").append(config.InitializationOption).append("\n");
        sb.append("\tInitialSteepestDescents:       ").append(config.InitialSteepestDescents).append("\n");
        sb.append("\tLinkCut:                       ").append(config.LinkCut).append("\n");
        sb.append("\tLocalVectorDimension:          ").append(config.LocalVectorDimension).append("\n");
        sb.append("\tMaxit:                         ").append(config.Maxit).append("\n");
        sb.append("\tMinimumDistance:               ").append(config.MinimumDistance).append("\n");
        sb.append("\tMPIIOStrategy:                 ").append(config.MPIIOStrategy).append("\n");
        sb.append("\tNbadgo:                        ").append(config.Nbadgo).append("\n");
        sb.append("\tOmega:                         ").append(config.Omega).append("\n");
        sb.append("\tOmegaOption:                   ").append(config.OmegaOption).append("\n");
        sb.append("\tPowerIterationLimit:           ").append(config.PowerIterationLimit).append("\n");
        sb.append("\tProcessingOption:              ").append(config.ProcessingOption).append("\n");
        sb.append("\tQgoodReductionFactor:          ").append(config.QgoodReductionFactor).append("\n");
        sb.append("\tQHighInitialFactor:            ").append(config.QHighInitialFactor).append("\n");
        sb.append("\tQLimiscalecalculationInterval: ").append(config.QLimitscalculationInterval).append("\n");
        sb.append("\tRotationOption:                ").append(config.RotationOption).append("\n");
        sb.append("\tSelectedfixedpoints:           ").append(config.Selectedfixedpoints).append("\n");
        sb.append("\tSelectedvariedpoints:          ").append(config.Selectedvariedpoints).append("\n");
        sb.append("\tStoredDistanceOption:          ").append(config.StoredDistanceOption).append("\n");
        sb.append("\tTimeCutmillisec:               ").append(config.TimeCutmillisec).append("\n");
        sb.append("\tTransformMethod:               ").append(config.TransformMethod).append("\n");
        sb.append("\tTransformParameter:            ").append(config.TransformParameter).append("\n");
        sb.append("\tUndefindDistanceValue:         ").append(config.UndefinedDistanceValue).append("\n");
        sb.append("\tVariedPointCriterion:          ").append(config.VariedPointCriterion).append("\n");
        sb.append("\tWeightingOption:               ").append(config.WeightingOption).append("\n");
        sb.append("\tWrite2Das3D:                   ").append(config.Write2Das3D).append("\n");
        sb.append("\n\nDensity").append("\n");
        sb.append("\tAlpha:                         ").append(config.Alpha).append("\n");
        sb.append("\tPcutf:                         ").append(config.Pcutf).append("\n");
        sb.append("\tSelectedClusters:              ").append(config.SelectedClusters).append("\n");
        sb.append("\tXmaxBound:                     ").append(config.XmaxBound).append("\n");
        sb.append("\tXres:                          ").append(config.Xres).append("\n");
        sb.append("\tYmaxBound:                     ").append(config.YmaxBound).append("\n");
        sb.append("\tYres:                          ").append(config.Yres).append("\n");
        return sb.toString();
    }

    private static String GenerateLinks(MDSasChisqSection config) {
        final String liTemplate = "<li><a href=\"%1$s\">%2$s</a></li>";
        final String txtExt = ".txt";
        String dotslash = tangible.DotNetToJavaStringHelper.isNullOrEmpty(
                config.ServerUrlPrefix) ? "./" : (config.ServerUrlPrefix + "/" + config.ManxcatRunName + "/");
        StringBuilder sb = new StringBuilder("<ul>");

        String reducedVectorFileNameTag = com.google.common.io.Files.getNameWithoutExtension(
                config.ReducedVectorOutputFileName);

        String link = dotslash + ManxcatConstants.SimplePointsPrefix + reducedVectorFileNameTag + txtExt;
        String name = "Simple Points";
        sb.append(String.format(liTemplate, link, name));

        link = dotslash + reducedVectorFileNameTag + ManxcatConstants.ColonPointsSuffix + txtExt;
        name = "Colon Points";
        sb.append(String.format(liTemplate, link, name));

        link = dotslash + reducedVectorFileNameTag + ManxcatConstants.GroupPointsSuffix + txtExt;
        name = "Group Points";
        sb.append(String.format(liTemplate, link, name));

        link = dotslash + (new java.io.File(config.SummaryOutputFileName)).getName();
        name = "Summary File";
        sb.append(String.format(liTemplate, link, name));

        link = dotslash + (new java.io.File(config.TimingOutputFileName)).getName();
        name = "Timing File";
        sb.append(String.format(liTemplate, link, name));

        link = dotslash + ManxcatConstants.OutFileNameTag + txtExt;
        name = "Manxcat Ouput";
        sb.append(String.format(liTemplate, link, name));

        link = dotslash + "whole-plot.png";
        name = "Density Graph";
        sb.append(String.format(liTemplate, link, name));

        if (SALSAUtility.IsClustersSelected) {
            link = dotslash + "selected-plot.png";
            name = "Density Graph for Selected Clusters (Intra Cluster Distances)";
            sb.append(String.format(liTemplate, link, name));

            link = dotslash + "selected-inter-plot.png";
            name = "Density Graph for Selected Clusters (Inter Cluster Distances)";
            sb.append(String.format(liTemplate, link, name));
        }

        return sb.toString();
    }

    private static void CopyStylesCSS(String toPath) {
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(Paths.get(toPath)))) {
            try (InputStream stream = ManxcatHtmlUtility.class.getClassLoader().getResourceAsStream("style")) {
                if (stream != null) {
                    try (BufferedReader reader = new BufferedReader(new InputStreamReader(stream))) {
                        String line;
                        while ((line = reader.readLine()) != null) {
                            writer.println(reader.readLine());
                        }
                    }
                } else {
                    SALSAUtility.printAndThrowRuntimeException(ManxcatErrorMessages.TemplateStyleSheetNotAvailable);
                }
            }
        } catch (IOException e) {
            SALSAUtility.printAndThrowRuntimeException(ManxcatErrorMessages.CopyingCSSFailed);
        }
    }
}