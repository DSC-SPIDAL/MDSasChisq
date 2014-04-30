package salsa.mdsaschisq;

import MDS.*;
import Salsa.Core.Configuration.Sections.*;
import SALSALibrary.*;

public class ManxcatHtmlUtility
{


	/** 
	 Generates a Web page linking to output files
	*/
	public static void WriteHTML()
	{
		ManxcatSection config = ManxcatCentral.config;
		String directoryPath = (new java.io.File(config.ReducedVectorOutputFileName)).getParent();
		directoryPath = tangible.DotNetToJavaStringHelper.isNullOrEmpty(directoryPath) ? "" : directoryPath;

		String stylesFile = Path.Combine(directoryPath, ManxcatConstants.StyleFileNameTag) + ".css";
		CopyStylesCSS(stylesFile);

		String htmlFile = Path.Combine(directoryPath, ManxcatConstants.HtmlFileNameTag) + ".html";
		GenerateHtmlContent(htmlFile, config);
	}

	private static void GenerateHtmlContent(String htmlFile, ManxcatSection config)
	{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//		using (StreamWriter writer = new StreamWriter(htmlFile))
		StreamWriter writer = new StreamWriter(htmlFile);
		try
		{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//			using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Manxcat.Web.template"))
			Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Manxcat.Web.template");
			try
			{
				if (stream != null)
				{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//					using (StreamReader reader = new StreamReader(stream))
					StreamReader reader = new StreamReader(stream);
					try
					{
						while (!reader.EndOfStream)
						{
							String line = reader.ReadLine();
							if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(line))
							{
								writer.WriteLine(!line.contains(ManxcatConstants.Tag) ? line : ProcessTagLine(line, config));
							}
						}
					}
					finally
					{
						reader.dispose();
					}
				}
				else
				{
					throw new RuntimeException(ManxcatErrorMessages.TemplateHtmlNotAvailable);
				}
			}
			finally
			{
				stream.dispose();
			}
		}
		finally
		{
			writer.dispose();
		}
	}

	private static String ProcessTagLine(String line, ManxcatSection config)
	{
		int tagStartIdx = line.indexOf(ManxcatConstants.Tag);
		int tagEndIdx = line.indexOf(ManxcatConstants.Tag, tagStartIdx + ManxcatConstants.Tag.length()) + (ManxcatConstants.Tag.length() - 1);
		String prefix = line.substring(0, tagStartIdx);
		String suffix = line.substring(tagEndIdx + 1);
		String meat = line.substring(tagStartIdx + ManxcatConstants.Tag.length(), tagStartIdx + ManxcatConstants.Tag.length() + ((tagEndIdx + 1) - (tagStartIdx + 2 * ManxcatConstants.Tag.length())));
//C# TO JAVA CONVERTER NOTE: The following 'switch' operated on a string member and was converted to Java 'if-else' logic:
//		switch (meat)
//ORIGINAL LINE: case ManxcatConstants.TagManxcatName:
		if (ManxcatConstants.TagManxcatName.equals(meat))
		{
		return prefix + config.ManxcatRunName + suffix;
		}
//ORIGINAL LINE: case ManxcatConstants.TagManxcatDescription:
		else if (ManxcatConstants.TagManxcatDescription.equals(meat))
		{
		return prefix + config.ManxcatRunDescription + suffix;
		}
//ORIGINAL LINE: case ManxcatConstants.TagManxcatConfig:
		else if (ManxcatConstants.TagManxcatConfig.equals(meat))
		{
		return GenerateConfig(config);
		}
//ORIGINAL LINE: case ManxcatConstants.TagManxcatLinks:
		else if (ManxcatConstants.TagManxcatLinks.equals(meat))
		{
		return prefix + GenerateLinks(config) + suffix;
		}
		else
		{
		throw new RuntimeException(ManxcatErrorMessages.UnidentifiedTagInHTMLGeneration);

		}

	}

	private static String GenerateConfig(ManxcatSection config)
	{
		StringBuilder sb = new StringBuilder();
		sb.AppendLine("I/O");
		sb.AppendLine("\tCoordinateWriteFrequency:      " + config.CoordinateWriteFrequency);
		sb.AppendLine("\tDistanceMatrixFile:            " + config.DistanceMatrixFile);
		sb.AppendLine("\n\nManxcatCore");
		sb.AppendLine("\tAddonforQcomputation:          " + config.AddonforQcomputation);
		sb.AppendLine("\tCalcFixedCrossFixed:           " + config.CalcFixedCrossFixed);
		sb.AppendLine("\tCGResidualLimit:               " + config.CGResidualLimit);
		sb.AppendLine("\tChisqChangePerPoint:           " + config.ChisqChangePerPoint);
		sb.AppendLine("\tChisqnorm:                     " + config.Chisqnorm);
		sb.AppendLine("\tChisqPrintConstant:            " + config.ChisqPrintConstant);
		sb.AppendLine("\tConversionInformation:         " + config.ConversionInformation);
		sb.AppendLine("\tConversionOption:              " + config.ConversionOption);
		sb.AppendLine("\tDataPoints:                    " + config.DataPoints);
		sb.AppendLine("\tDerivtest:                     " + config.Derivtest);
		sb.AppendLine("\tDiskDistanceOption:            " + config.DiskDistanceOption);
		sb.AppendLine("\tDistanceCut:                   " + config.DistanceCut);
		sb.AppendLine("\tDistanceFormula:               " + config.DistanceFormula);
		sb.AppendLine("\tDistanceProcessingOption:      " + config.DistanceProcessingOption);
		sb.AppendLine("\tDistanceWeigthsCuts:           " + config.DistanceWeightsCuts);
		sb.AppendLine("\tEigenvaluechange:              " + config.Eigenvaluechange);
		sb.AppendLine("\tEigenvectorchange:             " + config.Eigenvectorchange);
		sb.AppendLine("\tExtraOption1:                  " + config.ExtraOption1);
		sb.AppendLine("\tExtraprecision:                " + config.Extraprecision);
		sb.AppendLine("\tFixedPointCriterion:           " + config.FixedPointCriterion);
		sb.AppendLine("\tFletcherRho:                   " + config.FletcherRho);
		sb.AppendLine("\tFletcherSigma:                 " + config.FletcherSigma);
		sb.AppendLine("\tFullSecondDerivativeOption:    " + config.FullSecondDerivativeOption);
		sb.AppendLine("\tFunctionErrorCalcMultiplier:   " + config.FunctionErrorCalcMultiplier);
		sb.AppendLine("\tHistogramBinCount              " + config.HistogramBinCount);
		sb.AppendLine("\tInitializationLoops:           " + config.InitializationLoops);
		sb.AppendLine("\tInitializationOption:          " + config.InitializationOption);
		sb.AppendLine("\tInitialSteepestDescents:       " + config.InitialSteepestDescents);
		sb.AppendLine("\tLinkCut:                       " + config.LinkCut);
		sb.AppendLine("\tLocalVectorDimension:          " + config.LocalVectorDimension);
		sb.AppendLine("\tMaxit:                         " + config.Maxit);
		sb.AppendLine("\tMinimumDistance:               " + config.MinimumDistance);
		sb.AppendLine("\tMPIIOStrategy:                 " + config.MPIIOStrategy);
		sb.AppendLine("\tNbadgo:                        " + config.Nbadgo);
		sb.AppendLine("\tOmega:                         " + config.Omega);
		sb.AppendLine("\tOmegaOption:                   " + config.OmegaOption);
		sb.AppendLine("\tPowerIterationLimit:           " + config.PowerIterationLimit);
		sb.AppendLine("\tProcessingOption:              " + config.ProcessingOption);
		sb.AppendLine("\tQgoodReductionFactor:          " + config.QgoodReductionFactor);
		sb.AppendLine("\tQHighInitialFactor:            " + config.QHighInitialFactor);
		sb.AppendLine("\tQLimiscalecalculationInterval: " + config.QLimitscalculationInterval);
		sb.AppendLine("\tRotationOption:                " + config.RotationOption);
		sb.AppendLine("\tSelectedfixedpoints:           " + config.Selectedfixedpoints);
		sb.AppendLine("\tSelectedvariedpoints:          " + config.Selectedvariedpoints);
		sb.AppendLine("\tStoredDistanceOption:          " + config.StoredDistanceOption);
		sb.AppendLine("\tTimeCutmillisec:               " + config.TimeCutmillisec);
		sb.AppendLine("\tTransformMethod:               " + config.TransformMethod);
		sb.AppendLine("\tTransformParameter:            " + config.TransformParameter);
		sb.AppendLine("\tUndefindDistanceValue:         " + config.UndefinedDistanceValue);
		sb.AppendLine("\tVariedPointCriterion:          " + config.VariedPointCriterion);
		sb.AppendLine("\tWeightingOption:               " + config.WeightingOption);
		sb.AppendLine("\tWrite2Das3D:                   " + config.Write2Das3D);
		sb.AppendLine("\n\nDensity");
		sb.AppendLine("\tAlpha:                         " + config.Alpha);
		sb.AppendLine("\tPcutf:                         " + config.Pcutf);
		sb.AppendLine("\tSelectedClusters:              " + config.SelectedClusters);
		sb.AppendLine("\tXmaxBound:                     " + config.XmaxBound);
		sb.AppendLine("\tXres:                          " + config.Xres);
		sb.AppendLine("\tYmaxBound:                     " + config.YmaxBound);
		sb.AppendLine("\tYres:                          " + config.Yres);
		return sb.toString();
	}

	private static String GenerateLinks(ManxcatSection config)
	{
		final String liTemplate = "<li><a href=\"{0}\">{1}</a></li>";
		final String txtExt = ".txt";
		String dotslash = tangible.DotNetToJavaStringHelper.isNullOrEmpty(config.ServerUrlPrefix) ? "./" : (config.ServerUrlPrefix + "/" + config.ManxcatRunName + "/");
		StringBuilder sb = new StringBuilder("<ul>");

		String reducedVectorFileNameTag = Path.GetFileNameWithoutExtension(config.ReducedVectorOutputFileName);

		String link = dotslash + ManxcatConstants.SimplePointsPrefix + reducedVectorFileNameTag + txtExt;
		String name = "Simple Points";
		sb.AppendLine(String.format(liTemplate, link, name));

		link = dotslash + reducedVectorFileNameTag + ManxcatConstants.ColonPointsSuffix + txtExt;
		name = "Colon Points";
		sb.AppendLine(String.format(liTemplate, link, name));

		link = dotslash + reducedVectorFileNameTag + ManxcatConstants.GroupPointsSuffix + txtExt;
		name = "Group Points";
		sb.AppendLine(String.format(liTemplate, link, name));

		link = dotslash + (new java.io.File(config.SummaryOutputFileName)).getName();
		name = "Summary File";
		sb.AppendLine(String.format(liTemplate, link, name));

		link = dotslash + (new java.io.File(config.TimingOutputFileName)).getName();
		name = "Timing File";
		sb.AppendLine(String.format(liTemplate, link, name));

		link = dotslash + ManxcatConstants.OutFileNameTag + txtExt;
		name = "Manxcat Ouput";
		sb.AppendLine(String.format(liTemplate, link, name));

		link = dotslash + "whole-plot.png";
		name = "Density Graph";
		sb.AppendLine(String.format(liTemplate, link, name));

		if (SALSAUtility.IsClustersSelected)
		{
			link = dotslash + "selected-plot.png";
			name = "Density Graph for Selected Clusters (Intra Cluster Distances)";
			sb.AppendLine(String.format(liTemplate, link, name));

			link = dotslash + "selected-inter-plot.png";
			name = "Density Graph for Selected Clusters (Inter Cluster Distances)";
			sb.AppendLine(String.format(liTemplate, link, name));
		}

		return sb.toString();
	}

	private static void CopyStylesCSS(String toPath)
	{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//		using (StreamWriter writer = new StreamWriter(toPath))
		StreamWriter writer = new StreamWriter(toPath);
		try
		{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//			using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Manxcat.Web.style"))
			Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Manxcat.Web.style");
			try
			{
				if (stream != null)
				{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//					using (StreamReader reader = new StreamReader(stream))
					StreamReader reader = new StreamReader(stream);
					try
					{
						while (!reader.EndOfStream)
						{
							writer.WriteLine(reader.ReadLine());
						}
					}
					finally
					{
						reader.dispose();
					}
				}
				else
				{
					throw new RuntimeException(ManxcatErrorMessages.TemplateStyleSheetNotAvailable);
				}
			}
			finally
			{
				stream.dispose();
			}
		}
		finally
		{
			writer.dispose();
		}
	}

	public static void TestMethod()
	{
		String directoryName = (new java.io.File(ManxcatCentral.config.ReducedVectorOutputFileName)).getParent();
		String testfile = !tangible.DotNetToJavaStringHelper.isNullOrEmpty(directoryName) ? Path.Combine(directoryName, "TestMethod_rank_" + SALSAUtility.MPI_Rank + ".txt") : "TestMethod_rank_" + SALSAUtility.MPI_Rank + ".txt";
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//		using (StreamWriter writer = new StreamWriter(testfile))
		StreamWriter writer = new StreamWriter(testfile);
		try
		{
			for (int globalPointIndex = 0; globalPointIndex < SALSAUtility.PointCount_Global; globalPointIndex++)
			{
				if (ManxcatMDS.PointStatus[globalPointIndex] == -1)
				{
					continue;
				}
				int originalPointIndex = SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex];

				String coordinates = "";
				final int singleCluster = 1;
				int usedPointIndex = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex];

			   /* for (int globalPointIndex2 = 0; globalPointIndex2 < SALSAUtility.PointCount_Global; globalPointIndex2++)
			    {
			        if (ManxcatMDS.PointStatus[globalPointIndex2] == -1)
			            continue;
			        int originalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex2];
			        int usedPointIndex2 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex2];
			    }*/

				for (int localVectorIndex = 0; localVectorIndex < 3; localVectorIndex++)
				{
					coordinates += String.format("%.4E", Hotsun.GlobalParameter[usedPointIndex][localVectorIndex]) + "\t";
				}

				writer.WriteLine(String.format(globalPointIndex + "\t" + coordinates + singleCluster));

			}
		}
		finally
		{
			writer.dispose();
		}
	}
}