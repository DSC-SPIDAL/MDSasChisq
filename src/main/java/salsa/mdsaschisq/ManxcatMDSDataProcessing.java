package salsa.mdsaschisq;

import Manxcat.*;
import SALSALibrary.*;

//C# TO JAVA CONVERTER NOTE: There is no Java equivalent to C# namespace aliases:
//using LMatrix = MathNet.Numerics.LinearAlgebra.Matrix;
//C# TO JAVA CONVERTER NOTE: There is no Java equivalent to C# namespace aliases:
//using LU = MathNet.Numerics.LinearAlgebra.LUDecomposition;

public class ManxcatMDSDataProcessing
{
	public static void UpdateManxcatMDS_Option12()
	{
		if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("traditionaldirectory"))
		{
			UpdateManxcatMDS_Option12_TraditionalDirectory();
			return;
		}

		if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("cluster"))
		{
			UpdateManxcatMDS_Option12_Cluster();
			return;
		}

		if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("cutfile"))
		{
			UpdateManxcatMDS_Option12_FamilybyCuts();
			return;
		}

		if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("function"))
		{
			if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("statistics"))
			{
				UpdateManxcatMDS_Option12_Functions(true);
			}
			else
			{
				UpdateManxcatMDS_Option12_Functions(false);
			}
			return;
		}
	} // End  UpdateManxcatMDS_Option12()

	public static void UpdateManxcatMDS_Option12_Functions(boolean statistics)
	{
		//  Setup weight: Only used in statistics
		double[] weights = new double[SALSAUtility.PointCount_Global];
		ManxcatMDS.WeightingOption = ManxcatCentral.config.WeightingOption;
		ManxcatMDS.SetupWeightings(weights);
		String TypicalMDSFileName = "";

		//  Setup MDS: Only used in statistics
		double[][] MDSPoints = new double[SALSAUtility.PointCount_Global][Hotsun.ParameterVectorDimension];

		if (statistics)
		{ // Set up statistics by reading MDS file
			TypicalMDSFileName = ManxcatCentral.ResultDirectoryName + "\\SIMPLE" + ManxcatCentral.config.ReducedVectorOutputFileName;

			if (!(new java.io.File(TypicalMDSFileName)).isFile())
			{
				RuntimeException e = SALSAUtility.SALSAError(" File " + TypicalMDSFileName + " Does Not Exist");

				throw (e);
			}

			int NumberMDSPoints = 0;
			String[] InitialMDSString = new String[SALSAUtility.PointCount_Global];
			int[] ColorValue = new int[SALSAUtility.PointCount_Global];

			tangible.RefObject<Integer> tempRef_NumberMDSPoints = new tangible.RefObject<Integer>(NumberMDSPoints);
			boolean tempVar = ReadMDSCluster_File(TypicalMDSFileName, InitialMDSString, ColorValue, tempRef_NumberMDSPoints);
				NumberMDSPoints = tempRef_NumberMDSPoints.argValue;
			if (tempVar)
			{
				SALSAUtility.SALSAPrint(0, "MDS File " + TypicalMDSFileName + " Read Successfully");
			}

			for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
			{ // Extract MDS POints
				String inputLineStr = InitialMDSString[PointIndex].trim();
				String[] split = inputLineStr.split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);

				if (split.length < Hotsun.ParameterVectorDimension)
				{
					RuntimeException e = SALSAUtility.SALSAError(" Point " + (new Integer(PointIndex)).toString() + " has wrong number of Entries " + inputLineStr + " Entries" + split.length);

					throw (e);
				}

				for (int VectorIndex = 0; VectorIndex < Hotsun.ParameterVectorDimension; VectorIndex++)
				{
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

		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			if (PointIndex >= Limit)
			{
				++DivisionCount;
				Limit += BasicSize;

				if (LeftOver >= DivisionCount)
				{
					++Limit;
				}
			}
			ClusterLabel[PointIndex] = DivisionCount;
		}

		String FunctionFileName = ManxcatCentral.config.ConversionInformation;

		if (!FunctionFileName.contains(":"))
		{
			FunctionFileName = ManxcatCentral.config.ControlDirectoryName + "\\" + FunctionFileName;
		}

		double[] functionkeys = new double[SALSAUtility.PointCount_Global];
		int[] functionlabels = new int[SALSAUtility.PointCount_Global];

		int pickout = -1;

		while (true)
		{
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

			if (statistics)
			{
				for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
				{
					meanxyz[VectorIndex1] = 0.0;

					for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
					{
						varxyz[VectorIndex1][VectorIndex2] = 0.0;
					}
					correlxyzgamma[VectorIndex1] = 0.0;
				}
			}

			try
			{
				// Check if file exists
				if (!(new java.io.File(FunctionFileName)).isFile())
				{
					RuntimeException e = SALSAUtility.SALSAError("File " + FunctionFileName + " does not exists.");

					throw (e);
				}

				// Create a new stream to read from a file
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//				using (StreamReader sr = File.OpenText(FunctionFileName))
				StreamReader sr = File.OpenText(FunctionFileName);
				try
				{
					// Read contents of a file, line by line, into a string
					String inputLineStr;

					while ((inputLineStr = sr.ReadLine()) != null)
					{
						inputLineStr = inputLineStr.trim();
						String[] split = inputLineStr.split(new char[] {' '}, StringSplitOptions.RemoveEmptyEntries);

						if (split.length <= pickout)
						{
							break;
						}
						endofdata = false;
						String usethis = split[pickout].trim();

						if (usethis.length() < 1)
						{
							continue; //replace empty line
						}
						double gamma = Double.parseDouble(usethis);
						functionkeys[NumberofLines] = gamma;
						functionlabels[NumberofLines] = NumberofLines;
						sumabs += Math.abs(gamma);

						if (statistics)
						{
							double wgt = weights[NumberofLines];
							meangamma += gamma * wgt;
							vargamma += wgt * gamma * gamma;
							totalweight += wgt;

							for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
							{
								meanxyz[VectorIndex1] += MDSPoints[NumberofLines][VectorIndex1] * wgt;
								correlxyzgamma[VectorIndex1] += MDSPoints[NumberofLines][VectorIndex1] * gamma * wgt;

								for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
								{
									varxyz[VectorIndex1][VectorIndex2] += MDSPoints[NumberofLines][VectorIndex1] * MDSPoints[NumberofLines][VectorIndex2] * wgt;
								}
							}
						}
						++NumberofLines;
					}
					sr.Close();
				}
				finally
				{
					sr.dispose();
				}
			}
			catch (RuntimeException e)
			{
				SALSAUtility.SALSAError("Failed to read data from " + FunctionFileName + " " + e.toString());

				throw (e);
			}

			if (endofdata)
			{
				break;
			}

			if (NumberofLines != SALSAUtility.PointCount_Global)
			{
				SALSAUtility.SALSAError("Incorrect Function count read " + NumberofLines + " Expected " + SALSAUtility.PointCount_Global);
			}

			// Set Statistics
			if (statistics)
			{
				double[] alpha = new double[Hotsun.ParameterVectorDimension];
				double alphanorm = 0.0;
				meangamma = meangamma / totalweight;
				vargamma = (vargamma / totalweight) - meangamma * meangamma;

				for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
				{
					meanxyz[VectorIndex1] = meanxyz[VectorIndex1] / totalweight;
				}

				for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
				{
					for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
					{
						varxyz[VectorIndex1][VectorIndex2] = (varxyz[VectorIndex1][VectorIndex2] / totalweight) - meanxyz[VectorIndex1] * meanxyz[VectorIndex2];
					}
					correlxyzgamma[VectorIndex1] = correlxyzgamma[VectorIndex1] / totalweight - meangamma * meanxyz[VectorIndex1];
					alpha[VectorIndex1] = correlxyzgamma[VectorIndex1] / varxyz[VectorIndex1][VectorIndex1];
					alphanorm += alpha[VectorIndex1] * alpha[VectorIndex1];
				}

				// invert Matrix to find best alpha
				double[][] ConventionalFirst = new double[Hotsun.ParameterVectorDimension][1];
				double[][] ConventionalAnswer = new double[Hotsun.ParameterVectorDimension][];
				double[][] ConventionalMatrix = new double[Hotsun.ParameterVectorDimension][Hotsun.ParameterVectorDimension];

				for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
				{
					ConventionalAnswer[VectorIndex1] = new double[1];
					ConventionalFirst[VectorIndex1][0] = correlxyzgamma[VectorIndex1] / Math.sqrt(varxyz[VectorIndex1][VectorIndex1]);

					for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
					{
						ConventionalMatrix[VectorIndex1][VectorIndex2] = varxyz[VectorIndex1][VectorIndex2] / Math.sqrt(varxyz[VectorIndex1][VectorIndex1] * varxyz[VectorIndex2][VectorIndex2]);
					}
				}
				MathNet.Numerics.LinearAlgebra.Matrix cMatrix = MathNet.Numerics.LinearAlgebra.Matrix.Create(ConventionalMatrix);
				MathNet.Numerics.LinearAlgebra.Matrix RightHandSide = MathNet.Numerics.LinearAlgebra.Matrix.Create(ConventionalFirst);
				MathNet.Numerics.LinearAlgebra.Matrix MatrixAnswer = cMatrix.Solve(RightHandSide);
				ConventionalAnswer = (double[][])MatrixAnswer;

				alphanorm = 0.0;

				for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
				{
					alpha[VectorIndex1] = ConventionalAnswer[VectorIndex1][0] / Math.sqrt(varxyz[VectorIndex1][VectorIndex1]);
					alphanorm += alpha[VectorIndex1] * alpha[VectorIndex1];
				}

				double Fullcorrelation = 0.0;
				double varalphaxyz = 0.0;

				for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
				{
					alpha[VectorIndex1] = alpha[VectorIndex1] / Math.sqrt(alphanorm);
					Fullcorrelation += alpha[VectorIndex1] * correlxyzgamma[VectorIndex1];
				}

				for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
				{
					for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
					{
						varalphaxyz += alpha[VectorIndex1] * alpha[VectorIndex2] * varxyz[VectorIndex1][VectorIndex2];
					}
				}
				System.out.println(varalphaxyz + " " + Fullcorrelation + " " + vargamma);
				Fullcorrelation = Fullcorrelation / (Math.sqrt(vargamma * varalphaxyz));
				SALSAUtility.SALSAPrint(0, "Column " + (new Integer(pickout)).toString() + " File " + FunctionFileName + " MDS File " + TypicalMDSFileName + " Total Weight " + (new Double(totalweight)).toString() + " Correlation " + String.format("%0.4f", Fullcorrelation) + " Direction " + String.format("%0.4f", alpha[0]) + " " + String.format("%0.4f", alpha[1]) + " " + String.format("%0.4f", alpha[2]));
			}

			// Set Divisions
			Arrays.sort(functionkeys, functionlabels);
			int[] PointColors = new int[SALSAUtility.PointCount_Global];

			for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
			{
				int UnsortedPosition = functionlabels[PointIndex];
				PointColors[UnsortedPosition] = ClusterLabel[PointIndex];
			}
			String OutputFileExtension = ManxcatCentral.config.ConversionInformation;
			OutputFileExtension = OutputFileExtension.replace(".dat", "");
			OutputFileExtension = OutputFileExtension.replace(".txt", "");
			OutputFileExtension = OutputFileExtension + "-" + (new Integer(pickout)).toString() + ".txt";
			String labelFileDirectory = ManxcatCentral.config.ClusterDirectory; // Directory for Colors attached to Points
			String OutputFileName = labelFileDirectory + "\\" + OutputFileExtension;
			WritePointCluster(OutputFileName, PointColors, SALSAUtility.PointCount_Global);

			double test = 0.1 * sumabs / NumberofLines;

			for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
			{
				int UnsortedPosition = functionlabels[PointIndex];

				if (Math.abs(functionkeys[UnsortedPosition]) < test)
				{
					PointColors[UnsortedPosition] = 0;
				}
			}

			OutputFileName = OutputFileName.replace(".txt", "NZ.txt");
			WritePointCluster(OutputFileName, PointColors, SALSAUtility.PointCount_Global);
		} // End while over pickout (columns in data)
		return;
	} // End UpdateManxcatMDS_Option12_Functions

	public static void UpdateManxcatMDS_Option12_FamilybyCuts()
	{
		String CutFileName = ManxcatCentral.config.ConversionInformation;

		if (!CutFileName.contains(":"))
		{
			CutFileName = ManxcatCentral.config.ControlDirectoryName + "\\" + CutFileName;
		}

		String[] inputdata = new String[SALSAUtility.PointCount_Global];
		int NumberofCuts = 0;

		String FamilyLabel = "Family";
		int FamilyNumber = 1;

		if (ManxcatCentral.config.ConversionOption.toLowerCase().contains("family2"))
		{
			FamilyNumber = 2;
		}

		tangible.RefObject<Integer> tempRef_NumberofCuts = new tangible.RefObject<Integer>(NumberofCuts);
		ReadString_File(CutFileName, inputdata, tempRef_NumberofCuts);
		NumberofCuts = tempRef_NumberofCuts.argValue;

		if (NumberofCuts <= 0)
		{
			RuntimeException e = SALSAUtility.SALSAError("Too few points " + (new Integer(NumberofCuts)).toString() + " in Cut File " + CutFileName);

			throw (e);
		}

		int[] CutStart = new int[NumberofCuts];
		String[] CutLabel = new String[NumberofCuts];
		int CountFamilies = 0;
		int itmpold = -2;

		for (int CutPos = 0; CutPos < NumberofCuts; CutPos++)
		{
			String[] split = inputdata[CutPos].split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
			int itmp = Integer.parseInt(split[0]);

			if (itmp < 0)
			{
				FamilyLabel = split[1];
				continue;
			}

			if (itmp <= itmpold)
			{
				continue;
			}
			CutStart[CountFamilies] = itmp;
			CutLabel[CountFamilies] = FamilyLabel + "-" + split[1];
			++CountFamilies;
			itmpold = itmp;
		}

		if (FamilyNumber == 1)
		{
			SALSAUtility.GlobalFileProperties.FamilyName1 = FamilyLabel + " from " + CutFileName;
		}

		if (FamilyNumber == 2)
		{
			SALSAUtility.GlobalFileProperties.FamilyName2 = FamilyLabel + " from " + CutFileName;
		}

		int minlabel = 1;
		int CurrentFamily = 0;
		int[] CategoryIndex = new int[SALSAUtility.PointCount_Global];

		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			if (CurrentFamily < CountFamilies - 1)
			{
				if ((PointIndex + 1) >= CutStart[CurrentFamily + 1])
				{
					++CurrentFamily;
				}
			}
			CategoryIndex[PointIndex] = minlabel + CurrentFamily;
			continue;
		}

		int[] FamilyCounts = new int[CountFamilies];

		for (int Icount = 0; Icount < CountFamilies; Icount++)
		{
			FamilyCounts[Icount] = 0;
		}

		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			++FamilyCounts[CategoryIndex[PointIndex] - minlabel];
		}

		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			int Icount = CategoryIndex[PointIndex];

			if (FamilyNumber == 1)
			{
				SALSAUtility.GlobalPointProperties[PointIndex].family1 = Icount;
				SALSAUtility.GlobalPointProperties[PointIndex].familylabel1 = CutLabel[Icount - minlabel] + " " + (new Integer(Icount)).toString() + " out of " + (new Integer(CountFamilies)).toString() + " Count " + (new Integer(FamilyCounts[Icount - minlabel])).toString();
			}

			if (FamilyNumber == 2)
			{
				SALSAUtility.GlobalPointProperties[PointIndex].family2 = Icount;
				SALSAUtility.GlobalPointProperties[PointIndex].familylabel2 = CutLabel[Icount - minlabel] + " " + (new Integer(Icount)).toString() + " out of " + (new Integer(CountFamilies)).toString() + " Count " + (new Integer(FamilyCounts[Icount - minlabel])).toString();
			}
		}

		String cleandate = SALSAUtility.startTime.ToLocalTime().toString();
		cleandate = cleandate.replace(":", ".");
		cleandate = cleandate.replace(" ", "-");
		String OldComment = SALSAUtility.GlobalFileProperties.Comment;

		if (!OldComment.equals(""))
		{
			OldComment += "\n";
		}
		OldComment += "Family" + (new Integer(FamilyNumber)).toString() + " " + FamilyLabel + " Information added from " + CutFileName + " Time " + cleandate;
		SALSAUtility.GlobalFileProperties.Comment = OldComment;

		//  Write new label file
		String labelfilename = ManxcatCentral.config.DataLabelsFileName;

		if (!labelfilename.contains(":"))
		{
			labelfilename = ManxcatCentral.config.ControlDirectoryName + "\\" + labelfilename;
		}
		SALSALibrary.SALSA_Properties.WriteDataPointFile(labelfilename, ManxcatCentral.config.Write2Das3D, "colon,SameFileName ", SALSAUtility.GlobalFileProperties, SALSAUtility.GlobalPointProperties, SALSAUtility.NumberOriginalPoints);
		SALSAUtility.SALSAPrint(0, "Family" + (new Integer(FamilyNumber)).toString() + " " + FamilyLabel + " Info Added to " + labelfilename + " from " + CutFileName);
		ManxcatCentral.config.Comment += "\nFamily" + (new Integer(FamilyNumber)).toString() + " " + FamilyLabel + " Info Added to " + labelfilename + " from " + CutFileName;
	} // End UpdateManxcatMDS_Option12_FamilybyCuts

	public static void UpdateManxcatMDS_Option12_Cluster()
	{
		//  Read Cluster Information with File specified in ConversionInformation
		String ClusterFileName = ManxcatCentral.config.ConversionInformation;
		int NumberOfClusterLines = 0;
		String[] CategoryLabel = new String[SALSAUtility.PointCount_Global];
		tangible.RefObject<Integer> tempRef_NumberOfClusterLines = new tangible.RefObject<Integer>(NumberOfClusterLines);
		ReadClusterLabel_File(ClusterFileName, CategoryLabel, tempRef_NumberOfClusterLines);
		NumberOfClusterLines = tempRef_NumberOfClusterLines.argValue;

		if (NumberOfClusterLines != SALSAUtility.PointCount_Global)
		{
			RuntimeException e = SALSAUtility.SALSAError(" Illegal Count " + (new Integer(NumberOfClusterLines)).toString() + " in File " + ClusterFileName);

			throw (e);
		}

		int[] CategoryIndex = new int[SALSAUtility.PointCount_Global];
		int minlabel = SALSAUtility.PointCount_Global;
		int maxlabel = 0;

		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			CategoryIndex[PointIndex] = Integer.parseInt(CategoryLabel[PointIndex]);

			if (minlabel > CategoryIndex[PointIndex])
			{
				minlabel = CategoryIndex[PointIndex];
			}

			if (maxlabel < CategoryIndex[PointIndex])
			{
				maxlabel = CategoryIndex[PointIndex];
			}
		}

		int TotalNumberClusters = maxlabel - minlabel + 1;
		int[] ClusterCounts = new int[TotalNumberClusters];

		for (int Icount = 0; Icount < TotalNumberClusters; Icount++)
		{
			ClusterCounts[Icount] = 0;
		}

		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			++ClusterCounts[CategoryIndex[PointIndex] - minlabel];
		}

		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			int Icount = CategoryIndex[PointIndex];
			SALSAUtility.GlobalPointProperties[PointIndex].cluster = Icount;
			SALSAUtility.GlobalPointProperties[PointIndex].clusterlabel = "Cluster " + (new Integer(Icount)).toString() + " out of " + (new Integer(TotalNumberClusters)).toString() + " Count " + (new Integer(ClusterCounts[Icount - minlabel])).toString();
		}
		SALSAUtility.GlobalFileProperties.ClusterName = "Clusters from " + ClusterFileName;
		SALSAUtility.GlobalFileProperties.ClusterStartIndex = minlabel;

		String cleandate = SALSAUtility.startTime.ToLocalTime().toString();
		cleandate = cleandate.replace(":", ".");
		cleandate = cleandate.replace(" ", "-");
		String OldComment = SALSAUtility.GlobalFileProperties.Comment;

		if (!OldComment.equals(""))
		{
			OldComment += "\n";
		}
		OldComment += "Cluster Information added from " + ClusterFileName + " Time " + cleandate;
		SALSAUtility.GlobalFileProperties.Comment = OldComment;

		//  Write new label file
		String labelfilename = ManxcatCentral.config.DataLabelsFileName;

		if (!labelfilename.contains(":") && !labelfilename.contains("$"))
		{
			labelfilename = ManxcatCentral.config.ControlDirectoryName + "\\" + labelfilename;
		}
		SALSALibrary.SALSA_Properties.WriteDataPointFile(labelfilename, ManxcatCentral.config.Write2Das3D, "colon,SameFileName ", SALSAUtility.GlobalFileProperties, SALSAUtility.GlobalPointProperties, SALSAUtility.NumberOriginalPoints);
		SALSAUtility.SALSAPrint(0, "Cluster Info Added to " + labelfilename + " from " + ClusterFileName);
		ManxcatCentral.config.Comment += "\nCluster Info Added to " + labelfilename + " from " + ClusterFileName;
	} // End UpdateManxcatMDS_Option12_Cluster()

	public static void UpdateManxcatMDS_Option12_TraditionalDirectory()
	{
		String labelFileDirectory = ManxcatCentral.config.ClusterDirectory; // Input Directory
		String MDSandClusterDirectory = labelFileDirectory + "\\" + ManxcatCentral.config.RunSetLabel + "-R" + ManxcatCentral.config.RunNumber.toString() + "-ManxcatMDS"; // Output Directory

		if ((new java.io.File(MDSandClusterDirectory)).isDirectory())
		{
			SALSAUtility.SALSAPrint(0, "The directory " + MDSandClusterDirectory + " exists");
		}
		else
		{
			java.io.File di = (new java.io.File(MDSandClusterDirectory)).mkdir();
			SALSAUtility.SALSAPrint(0, "The directory " + MDSandClusterDirectory + " was created successfully at " + Directory.GetCreationTime(MDSandClusterDirectory));
		}

		String[] SMACOFfileEntries = Directory.GetFiles(MDSandClusterDirectory);
		String TypicalMDSFileName = ManxcatCentral.config.ReducedVectorOutputFileName;

		if (!TypicalMDSFileName.contains(":"))
		{
			TypicalMDSFileName = ManxcatCentral.ResultDirectoryName + "\\SIMPLE" + TypicalMDSFileName;
		}

		if (!(new java.io.File(TypicalMDSFileName)).isFile())
		{
			RuntimeException e = SALSAUtility.SALSAError(" File " + TypicalMDSFileName + " Does Not Exist");

			throw (e);
		}

		int NumberMDSPoints = 0;
		String[] InitialMDSString = new String[SALSAUtility.PointCount_Global];
		int[] ColorValue = new int[SALSAUtility.PointCount_Global];

		tangible.RefObject<Integer> tempRef_NumberMDSPoints = new tangible.RefObject<Integer>(NumberMDSPoints);
		boolean tempVar = ReadMDSCluster_File(TypicalMDSFileName, InitialMDSString, ColorValue, tempRef_NumberMDSPoints);
			NumberMDSPoints = tempRef_NumberMDSPoints.argValue;
		if (tempVar)
		{
			System.out.println(TypicalMDSFileName + " Read Successfully");
		}

		String[] fileEntries = Directory.GetFiles(labelFileDirectory);

		for (String LabelFileName : fileEntries)
		{
			if (LabelFileName.contains(".dot"))
			{
				continue;
			}

			if (LabelFileName.contains("Summary"))
			{
				continue;
			}

			if (LabelFileName.contains("Timing"))
			{
				continue;
			}
			String LabelFileName1 = LabelFileName.replace(labelFileDirectory + "\\", "");
			String coreFileName = MDSandClusterDirectory + "\\" + "MDSManxcat-" + LabelFileName1;

			if ((new java.io.File(coreFileName)).isFile())
			{
				continue;
			}
			int NumberOfLabels = 0;
			String[] CategoryLabel = new String[SALSAUtility.PointCount_Global];
			tangible.RefObject<Integer> tempRef_NumberOfLabels = new tangible.RefObject<Integer>(NumberOfLabels);
			ReadClusterLabel_File(LabelFileName, CategoryLabel, tempRef_NumberOfLabels);
			NumberOfLabels = tempRef_NumberOfLabels.argValue;

			if (NumberOfLabels != SALSAUtility.PointCount_Global)
			{
				RuntimeException e = SALSAUtility.SALSAError(" Illegal Count " + (new Integer(NumberOfLabels)).toString() + " in File " + LabelFileName);

				throw (e);
			}

			int[] CategoryIndex = new int[SALSAUtility.PointCount_Global];

			for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
			{
				CategoryIndex[PointIndex] = Integer.parseInt(CategoryLabel[PointIndex]);
			}
			WriteColor_Cluster(coreFileName, InitialMDSString, CategoryIndex, SALSAUtility.PointCount_Global, false);
			SALSAUtility.SALSAPrint(0, "Traditional Directory MDS Info Added to " + coreFileName);
		}
	} // End UpdateManxcatMDS_Option12_TraditionalDirectory

	public static boolean ReadMDSCluster_File(String MDSClusterFileName, String[] InitialString, int[] ColorValue, tangible.RefObject<Integer> NumberofPoints)
	{
		NumberofPoints.argValue = 0;

		try
		{
			// Check if file exists
			if (!(new java.io.File(MDSClusterFileName)).isFile())
			{
				RuntimeException e = SALSAUtility.SALSAError("File " + MDSClusterFileName + " does not exists");

				throw (e);
			}

			// Create a new stream to read from a file
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//			using (StreamReader sr = File.OpenText(MDSClusterFileName))
			StreamReader sr = File.OpenText(MDSClusterFileName);
			try
			{
				// Read contents of a file, line by line, into a string
				String inputLineStr;
				int newlabelnumber = -1;

				while ((inputLineStr = sr.ReadLine()) != null)
				{
					if (inputLineStr.length() < 2)
					{
						continue; //replace empty line
					}
					inputLineStr = tangible.DotNetToJavaStringHelper.trim(inputLineStr, new char[] {' ', '\t'});

					try
					{
						// Parse each record string
						inputLineStr = inputLineStr.replace("\t\t", "\t");
						String[] split = inputLineStr.split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
						newlabelnumber = Integer.parseInt(split[0]);

						if (split.length != Hotsun.ParameterVectorDimension + 2)
						{
							System.out.println(" Bad Line " + split.length + " " + NumberofPoints.argValue + " " + inputLineStr);
						}

						// Begin Changes saliya 3/22/11
						// Note. This seems unnecessary. Let's agree on zero based indices
//                            if ((NumberofPoints + 1) != newlabelnumber)
//                            {
//                                Exception e = SALSAUtility.SALSAError("Unexpected Label Number "
//                                    + newlabelnumber + " Expected " + NumberofPoints + "+1");
//
//                                throw (e);
//                            }
						// End Changes saliya 3/22/11

						inputLineStr = split[0];

						for (int i = 1; i < (split.length - 1); i++)
						{
							inputLineStr += "\t" + split[i];
						}
						InitialString[NumberofPoints.argValue] = inputLineStr;

						ColorValue[NumberofPoints.argValue] = (int)Double.parseDouble(split[split.length - 1]);
						++NumberofPoints.argValue;
					}
					catch (RuntimeException e)
					{
						SALSAUtility.SALSAError("Failed to load data array " + inputLineStr + " " + " " + (new Integer(NumberofPoints.argValue)).toString() + " " + (new Integer(newlabelnumber)).toString() + " " + e.toString());

						throw (e);
					}
				}

				sr.Close();
			}
			finally
			{
				sr.dispose();
			}
		}
		catch (RuntimeException e)
		{
			SALSAUtility.SALSAError("Failed to read data from " + MDSClusterFileName + " " + e.toString());

			throw (e);
		}

		if (NumberofPoints.argValue != SALSAUtility.PointCount_Global)
		{
			RuntimeException e = SALSAUtility.SALSAError("Incorrect #Points in File " + NumberofPoints.argValue + " Expected " + SALSAUtility.PointCount_Global);

			throw (e);
		}

		return true;
	} // End ReadMDSClusterFile

	public static boolean ReadClusterLabel_File(String LabelFileName, String[] LabelArray, tangible.RefObject<Integer> NumberofLabels)
	{
		NumberofLabels.argValue = 0;
		int startnumber = 1;

		try
		{
			// Check if file exists
			if (!(new java.io.File(LabelFileName)).isFile())
			{
				RuntimeException e = SALSAUtility.SALSAError("File " + LabelFileName + " does not exists.");

				throw (e);
			}

			// Create a new stream to read from a file
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//			using (StreamReader sr = File.OpenText(LabelFileName))
			StreamReader sr = File.OpenText(LabelFileName);
			try
			{
				// Read contents of a file, line by line, into a string
				String inputLineStr;

				while ((inputLineStr = sr.ReadLine()) != null)
				{
					if (inputLineStr.length() < 2)
					{
						continue; //replace empty line
					}

					inputLineStr = inputLineStr.trim();

					try
					{
						// Parse each record string
						String[] split = inputLineStr.split(new char[] {' ', '\t'}, 2, StringSplitOptions.RemoveEmptyEntries);
						int newlabelnumber = Integer.parseInt(split[0]);

						if (NumberofLabels.argValue == 0)
						{
							startnumber = newlabelnumber;

							if ((startnumber < 0) || (startnumber > 1))
							{
								RuntimeException e = SALSAUtility.SALSAError("Unexpected Start Number " + (new Integer(startnumber)).toString());

								throw (e);
							}
						}

						if (NumberofLabels.argValue != (newlabelnumber + startnumber))
						{
							RuntimeException e = SALSAUtility.SALSAError("Unexpected Label Number " + (new Integer(newlabelnumber)).toString() + " Expected " + (new Integer(NumberofLabels.argValue)).toString() + " + " + (new Integer(startnumber)).toString());

							throw (e);
						}

						if (split[1].length() <= 0)
						{
							RuntimeException e = SALSAUtility.SALSAError("Zero length label for point " + (new Integer(NumberofLabels.argValue)).toString());

							throw (e);
						}

						LabelArray[NumberofLabels.argValue] = split[1];
						++NumberofLabels.argValue;
					}
					catch (RuntimeException e)
					{
						SALSAUtility.SALSAError("Failed to load data from " + LabelFileName + " " + e.toString());

						throw (e);
					}
				}

				sr.Close();
			}
			finally
			{
				sr.dispose();
			}
		}
		catch (RuntimeException e)
		{
			SALSAUtility.SALSAError("Failed to read data from " + LabelFileName + " " + e.toString());

			throw (e);
		}

		return true;
	} // End ReadClusterLabel_File


	// Write label-cluster results into a file
	public static void WriteColor_Cluster(String fname, String[] labels, int[] ColorValues, int dataPoints, boolean append)
	{
		try
		{
			StreamWriter sw = null;

			if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(fname))
			{
				sw = new StreamWriter(fname, append, Encoding.UTF8);
			}

			if (sw != null)
			{
				for (int i = 0; i < dataPoints; i++)
				{
					String stripped = tangible.DotNetToJavaStringHelper.trim(labels[i], new char[] {' ', '\t'});
					sw.WriteLine(String.format(stripped + "\t" + (new Integer(ColorValues[i])).toString() + ".0000000000"));
				}
			}

			sw.Flush();
			sw.Close();
		}
		catch (RuntimeException e)
		{
			SALSAUtility.SALSAError("Failed writing data on " + fname + " " + e.toString());

			throw (e);
		}
	} // End WriteColor_Cluster

	// Write Point-cluster results into a file
	public static void WritePointCluster(String fname, int[] ColorValues, int dataPoints)
	{
		try
		{
			StreamWriter sw = null;

			if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(fname))
			{
				sw = new StreamWriter(fname, false, Encoding.UTF8);
			}

			if (sw != null)
			{
				for (int i = 0; i < dataPoints; i++)
				{
					sw.WriteLine(String.format((new Integer(i + 1)).toString() + " " + (new Integer(ColorValues[i])).toString()));
				}
			}
			sw.Flush();
			sw.Close();
		}
		catch (RuntimeException e)
		{
			SALSAUtility.SALSAError("Failed writing data on " + fname + " " + e.toString());

			throw (e);
		}
	} // End WritePointCluster

	public static boolean ReadString_File(String StringFileName, String[] LabelArray, tangible.RefObject<Integer> NumberofLines)
	{
		NumberofLines.argValue = 0;

		try
		{
			// Check if file exists
			if (!(new java.io.File(StringFileName)).isFile())
			{
				RuntimeException e = SALSAUtility.SALSAError("File " + StringFileName + " does not exists.");

				throw (e);
			}

			// Create a new stream to read from a file
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//			using (StreamReader sr = File.OpenText(StringFileName))
			StreamReader sr = File.OpenText(StringFileName);
			try
			{
				// Read contents of a file, line by line, into a string
				String inputLineStr;

				while ((inputLineStr = sr.ReadLine()) != null)
				{
					if (inputLineStr.length() < 2)
					{
						continue; //replace empty line
					}

					LabelArray[NumberofLines.argValue] = inputLineStr.trim();
					++NumberofLines.argValue;
				}
				sr.Close();
			}
			finally
			{
				sr.dispose();
			}
		}
		catch (RuntimeException e)
		{
			SALSAUtility.SALSAError("Failed to read data from " + StringFileName + " " + e.toString());

			throw (e);
		}

		return true;
	} // End ReadString_File
}