package salsa.mdsaschisq;

public class SALSA_Properties
{


	public static boolean ReadDataPointFile(String MDSClusterFileName, tangible.RefObject<Integer> FileType, SALSAFileProperties FileProperties, tangible.RefObject<SALSADataPointProperties[]> DataPoints, tangible.RefObject<Integer> NumberofPoints)
	{
		// FileType = 0 simple integers and positions:  PointNumber x y z LABEL or Group Number (LABEL if file name contains label)
		// FileType = 1 Just integers -- Point Number and Group Number: PointNumber  Group Number 
		// File Type = 2 Integer and Label: Point LABEL (LABEL if file name contains label)
		// File Type = 3 Pure colon style -- with  positions
		// File Type = 4 Pure Colon Style -- no    positions
		// File Type = 5 Pure Colon Style --    Labels
		// File Type = 6 Hybrid Style --     with positions
		// File Type = 7 Hybrid Style --     no positions
		// File Type = 8 Hybrid Style --     Labels


		boolean positionsset = false;
		boolean colonsinput = false;
		boolean NonColonsInput = false;
		boolean labelfile = false;
		String LowerFileName = MDSClusterFileName.toLowerCase();
		if (LowerFileName.contains("label"))
		{
			labelfile = true;
		}

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
				// Read contents of a file
				String inputLineStr;
				int newlabelnumber = -1;
				int LocalStart = FileProperties.LocalPointStartIndex;
				int OriginalStart = FileProperties.OriginalPointStartIndex;
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

						if (inputLineStr.contains("FileProperties"))
						{ // Not a Data Point
							String EndofInput = inputLineStr.replace("FileProperties", "");
							ReadFileProperties(FileProperties, EndofInput);
							colonsinput = true;
							LocalStart = FileProperties.LocalPointStartIndex;
							OriginalStart = FileProperties.OriginalPointStartIndex;
						}
						else
						{ // Must be a Data Point
							ArrayInitializer(DataPoints, SALSAUtility.NumberOriginalPoints, FileProperties.NumberPointsinFile);
							DataPoints.argValue[NumberofPoints.argValue] = new SALSADataPointProperties();
							boolean incrementcount = false;
							int PointDataStarts = inputLineStr.indexOf("PointProperties");
							if (PointDataStarts >= 0)
							{ // Some Colon Information
								String EndofInput = inputLineStr.substring(PointDataStarts);
								EndofInput = EndofInput.replace("PointProperties", "");
								ReadDataPointProperties(DataPoints.argValue[NumberofPoints.argValue], EndofInput);
								colonsinput = true;
								if (DataPoints.argValue[NumberofPoints.argValue].valuesset)
								{
									positionsset = true;
								}
								incrementcount = true;
								DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber = NumberofPoints.argValue + FileProperties.LocalPointStartIndex;
							} //  End Processing Colon Point Information

							if (PointDataStarts < 0)
							{ // traditional bare line
								PointDataStarts = inputLineStr.length();
							}
							if (PointDataStarts > 0)
							{ // Process number information
								String StartofString = inputLineStr.substring(0, PointDataStarts);
								StartofString = tangible.DotNetToJavaStringHelper.trim(StartofString, new char[] {' ', '\t'});

								if (StartofString.length() > 0)
								{
									// You come here only for traditional bare line of x,y,z coordinates.

									String[] split = StartofString.split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);

									if ((split.length != 5) && (split.length != 2) && (split.length != 4))
									{
										RuntimeException e = SALSAUtility.SALSAError(" Bad Line " + split.length.toString() + " " + (new Integer(NumberofPoints.argValue)).toString() + " " + inputLineStr);
										throw (e);
									}
									newlabelnumber = Integer.parseInt(split[0]);

									if ((NumberofPoints.argValue + LocalStart) != newlabelnumber)
									{
										RuntimeException e = SALSAUtility.SALSAError("Unexpected Label Number " + newlabelnumber + " Expected " + (new Integer(NumberofPoints.argValue)).toString() + " + " + (new Integer(LocalStart)).toString());
										throw (e);
									}
									if (DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber >= 0)
									{
										if ((DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber - LocalStart) != NumberofPoints.argValue)
										{
											RuntimeException e = SALSAUtility.SALSAError("Unexpected Local Number " + DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber + " Expected " + (new Integer(NumberofPoints.argValue)).toString() + " + " + (new Integer(LocalStart)).toString());
											throw (e);
										}
									}
									DataPoints.argValue[NumberofPoints.argValue].LocalPointNumber = NumberofPoints.argValue;
									if (DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber >= 0)
									{
										if ((DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber - OriginalStart) < 0)
										{
											RuntimeException e = SALSAUtility.SALSAError("Unexpected Original Number " + DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber + " Local Point " + (new Integer(NumberofPoints.argValue)).toString() + " Original Increment " + (new Integer(OriginalStart)).toString());
											throw (e);
										}
										DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber -= OriginalStart;
									}
									else
									{
										DataPoints.argValue[NumberofPoints.argValue].OriginalPointNumber = newlabelnumber;
									}

									if (labelfile)
									{
										DataPoints.argValue[NumberofPoints.argValue].pointlabel = split[split.length - 1];
									}
									else
									{
										DataPoints.argValue[NumberofPoints.argValue].group = Integer.parseInt(split[split.length - 1]);
									}

									if (split.length >= 4)
									{
										DataPoints.argValue[NumberofPoints.argValue].valuesset = true;
										DataPoints.argValue[NumberofPoints.argValue].x = Double.parseDouble(split[1]);
										DataPoints.argValue[NumberofPoints.argValue].y = Double.parseDouble(split[2]);
										positionsset = true;
									}
									if (split.length == 5)
									{
										DataPoints.argValue[NumberofPoints.argValue].z = Double.parseDouble(split[3]);
									}
									incrementcount = true;
									NonColonsInput = true;
								} // End Processing non colon Point information
							}
							if (incrementcount)
							{
								++NumberofPoints.argValue;
							}
						}
					}
					catch (RuntimeException e)
					{
						SALSAUtility.SALSAError("Failed to load data array " + inputLineStr + " " + " " + (new Integer(NumberofPoints.argValue)).toString() + " " + (new Integer(newlabelnumber)).toString() + " " + e.toString());
						throw (e);
					}
				}

				FileType.argValue = 1;
				if (positionsset)
				{
					FileType.argValue = 0;
				}
				if (labelfile)
				{
					FileType.argValue = 2;
				}
				if (colonsinput)
				{
					if (NonColonsInput)
					{
						FileType.argValue += 6;
					}
					else
					{
						FileType.argValue += 3;
					}
				}
				if (FileProperties.NumberOriginalPoints == 0)
				{
					FileProperties.NumberOriginalPoints = NumberofPoints.argValue;
				}

				if (FileProperties.NumberPointsinFile == 0)
				{
					FileProperties.NumberPointsinFile = NumberofPoints.argValue;
				}

				if (FileProperties.NumberPointsinFile != NumberofPoints.argValue)
				{
					RuntimeException e = SALSAUtility.SALSAError("Unexpected Number of Points in File " + (new Integer(NumberofPoints.argValue)).toString() + " Read but Expected " + (new Integer(FileProperties.NumberPointsinFile)).toString());
					throw (e);
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
		return true;

	} // End ReadDataPointFile

	public static void ArrayInitializer(tangible.RefObject<SALSADataPointProperties[]> DataArray, int sizemax, int sizereadin)
	{
		if (DataArray.argValue != null)
		{
			if (DataArray.argValue.length < sizereadin)
			{
				RuntimeException e = SALSAUtility.SALSAError(" Data Array too small for file Length " + DataArray.argValue.length.toString() + " Needs " + (new Integer(sizereadin)).toString());
				throw (e);
			}
			return;
		}
		int size = sizereadin;
		if (size == 0)
		{
			size = sizemax;
		}
		DataArray.argValue = new SALSADataPointProperties[size];
		return;

	} // End ArrayInitializer

	// Write label-cluster results into a file
	public static void WriteDataPointFile(String CoreOutputFileName, boolean write2Das3D, String OutputStyles, SALSAFileProperties FileProperties, SALSADataPointProperties[] DataPoints, int NumberofDataPoints)
	{
		boolean[] DothisOutputStyle = new boolean[5];
		for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++)
		{
			DothisOutputStyle[StyleIndex] = false;
			if (OutputStyles.contains("all"))
			{
				DothisOutputStyle[StyleIndex] = true;
			}
		}
		if (OutputStyles.contains("colon"))
		{
			DothisOutputStyle[0] = true;
		}
		if (OutputStyles.contains("family1"))
		{
			DothisOutputStyle[1] = true;
		}
		if (OutputStyles.contains("family2"))
		{
			DothisOutputStyle[2] = true;
		}
		if (OutputStyles.contains("cluster"))
		{
			DothisOutputStyle[3] = true;
		}
		if (OutputStyles.contains("group"))
		{
			DothisOutputStyle[4] = true;
		}

		boolean setgroup = false;
		boolean OutputValues = false;

		SALSADataPointProperties FirstOne = DataPoints[0];
		if (FirstOne.family1 < 0)
		{
			DothisOutputStyle[1] = false;
		}
		if (FirstOne.family2 < 0)
		{
			DothisOutputStyle[2] = false;
		}
		if (FirstOne.cluster < 0)
		{
			DothisOutputStyle[3] = false;
		}
		if (FirstOne.group < 0)
		{
			DothisOutputStyle[4] = false;
		}
		if ((FirstOne.family1 < 0) && (FirstOne.family2 < 0) && (FirstOne.cluster < 0) && (FirstOne.group < 0))
		{
			DothisOutputStyle[4] = true;
			setgroup = true;
		}
		if (FirstOne.valuesset)
		{
			OutputValues = true;
		}
		int LocalVectorDimension = FileProperties.LocalVectorDimension;
		int LocalPointIncrement = FileProperties.LocalPointStartIndex;

		int numberoffiles = 0;
		for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++)
		{
			if (DothisOutputStyle[StyleIndex])
			{
				++numberoffiles;
			}
		}
		if (numberoffiles == 0)
		{
			SALSAUtility.SALSAError("No files output for core name " + CoreOutputFileName);
			return;
		}
		if (numberoffiles > 1)
		{
			if (OutputStyles.contains("SameFileName"))
			{
				RuntimeException e = SALSAUtility.SALSAError("Attempt to generate multiple outputs to same file " + CoreOutputFileName);
				throw (e);
			}
		}
		for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++)
		{
			if (!DothisOutputStyle[StyleIndex])
			{
				continue;
			}
			String OutputFileName = "";
			if (!OutputStyles.contains("SameFileName"))
			{
				if (StyleIndex == 0)
				{
					OutputFileName = CoreOutputFileName.replace(".txt", "Colon.txt");
				}
				if (StyleIndex == 1)
				{
					OutputFileName = CoreOutputFileName.replace(".txt", "Family1.txt");
				}
				if (StyleIndex == 2)
				{
					OutputFileName = CoreOutputFileName.replace(".txt", "Family2.txt");
				}
				if (StyleIndex == 3)
				{
					OutputFileName = CoreOutputFileName.replace(".txt", "Cluster.txt");
				}
				if (StyleIndex == 4)
				{
					OutputFileName = CoreOutputFileName.replace(".txt", "Group.txt");
				}
			}
			else
			{
				OutputFileName = CoreOutputFileName;
			}

			try
			{

				StreamWriter sw = null;
				if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(OutputFileName))
				{
					sw = new StreamWriter(OutputFileName, false, Encoding.UTF8);
				}
				if (sw != null)
				{
					if (StyleIndex == 0)
					{
						WriteFileProperties(FileProperties, sw); // Write Header of a Colon File
					}

					for (int GlobalDataPoint = 0; GlobalDataPoint < NumberofDataPoints; GlobalDataPoint++)
					{
						SALSADataPointProperties ThesePointProperties = DataPoints[GlobalDataPoint];
						if ((LocalVectorDimension == 2) && write2Das3D && OutputValues)
						{
							ThesePointProperties.z = 0.0;
							ThesePointProperties.zerr = 0.0;
						}
						if (StyleIndex == 0)
						{
							String OutputLine = "";
							tangible.RefObject<String> tempRef_OutputLine = new tangible.RefObject<String>(OutputLine);
							AppendDataPointProperties(ThesePointProperties, tempRef_OutputLine);
							OutputLine = tempRef_OutputLine.argValue;
							sw.WriteLine(OutputLine);
							continue;
						}
						int IntegerIndex = 0;
						if (StyleIndex == 1)
						{
							IntegerIndex = ThesePointProperties.family1;
						}
						if (StyleIndex == 2)
						{
							IntegerIndex = ThesePointProperties.family2;
						}
						if (StyleIndex == 3)
						{
							IntegerIndex = ThesePointProperties.cluster;
						}
						if (StyleIndex == 4)
						{
							IntegerIndex = ThesePointProperties.group;
							if (setgroup)
							{
								IntegerIndex = 1;
							}
						}
						String Coordinates = "";
						if (OutputValues)
						{
							Coordinates = String.format("%.4E", ThesePointProperties.x) + "\t" + String.format("%.4E", ThesePointProperties.y) + "\t";
							if ((LocalVectorDimension == 3) || write2Das3D)
							{
								Coordinates += String.format("%.4E", ThesePointProperties.z) + "\t";
							}
						}
						sw.WriteLine(String.format((GlobalDataPoint + LocalPointIncrement).toString() + "\t" + Coordinates + (new Integer(IntegerIndex)).toString()));
					}

					sw.Flush();
					sw.Close();
				}
			}
			catch (java.lang.Exception e)
			{
				RuntimeException e = SALSAUtility.SALSAError(" Failed to Write Properties File " + CoreOutputFileName);
				throw (e);
			}

		} // End Loop over File Types
		return;

	} // End WriteDataPointFile


	public static void ReadFileProperties(SALSAFileProperties FileProps, String InputLine)
	{
		InputLine = tangible.DotNetToJavaStringHelper.trim(InputLine, new char[] {' ', '\t'});
		String[] colonsplit = InputLine.split(new char[] {':'}, 2);
		if (colonsplit.length < 2)
		{
			System.out.println("No deliminator in Line " + InputLine);
			return;
		}
		colonsplit[0] = colonsplit[0].trim();
		colonsplit[1] = colonsplit[1].trim();
		if (colonsplit[0].equals("LocalVectorDimension"))
		{
			FileProps.LocalVectorDimension = Integer.parseInt(colonsplit[1]);
		}
		if (colonsplit[0].equals("ClusterStartIndex"))
		{
			FileProps.ClusterStartIndex = Integer.parseInt(colonsplit[1]);
		}
		if (colonsplit[0].equals("OriginalPointStartIndex"))
		{
			FileProps.OriginalPointStartIndex = Integer.parseInt(colonsplit[1]);
		}
		if (colonsplit[0].equals("LocalPointStartIndex"))
		{
			FileProps.LocalPointStartIndex = Integer.parseInt(colonsplit[1]);
		}
		if (colonsplit[0].equals("FileGenerationType"))
		{
			FileProps.FileGenerationType = Integer.parseInt(colonsplit[1]);
		}
		if (colonsplit[0].equals("FamilyName1"))
		{
			FileProps.FamilyName1 = colonsplit[1];
		}
		if (colonsplit[0].equals("FamilyName2"))
		{
			FileProps.FamilyName2 = colonsplit[1];
		}
		if (colonsplit[0].equals("GroupName"))
		{
			FileProps.GroupName = colonsplit[1];
		}
		if (colonsplit[0].equals("ClusterName"))
		{
			FileProps.ClusterName = colonsplit[1];
		}
		if (colonsplit[0].equals("NumberOriginalPoints"))
		{
			FileProps.NumberOriginalPoints = Integer.parseInt(colonsplit[1]);
		}
		if (colonsplit[0].equals("NumberPointsinFile"))
		{
			FileProps.NumberPointsinFile = Integer.parseInt(colonsplit[1]);
		}
		if (colonsplit[0].equals("RotationParameters"))
		{
			FileProps.RotationParameters = colonsplit[1];
		}
		if (colonsplit[0].equals("NumberRotationParameters"))
		{
			FileProps.NumberRotationParameters = Integer.parseInt(colonsplit[1]);
		}
		if (colonsplit[0].equals("Comment"))
		{
			if (!FileProps.Comment.equals(""))
			{
				FileProps.Comment += "\n";
			}
			FileProps.Comment += colonsplit[1];
		}
		return;

	} // End ReadFileProperties

	public static void WriteFileProperties(SALSAFileProperties FileProps, StreamWriter sw)
	{
		sw.WriteLine(String.format("FileProperties\tLocalVectorDimension:%1$s", FileProps.LocalVectorDimension));
		sw.WriteLine(String.format("FileProperties\tClusterStartIndex:%1$s", FileProps.ClusterStartIndex));
		sw.WriteLine(String.format("FileProperties\tOriginalPointStartIndex:%1$s", FileProps.OriginalPointStartIndex));
		sw.WriteLine(String.format("FileProperties\tLocalPointStartIndex:%1$s", FileProps.LocalPointStartIndex));
		sw.WriteLine(String.format("FileProperties\tFileGenerationType:%1$s", FileProps.FileGenerationType));
		sw.WriteLine(String.format("FileProperties\tFamilyName1:%1$s", FileProps.FamilyName1));
		sw.WriteLine(String.format("FileProperties\tFamilyName2:%1$s", FileProps.FamilyName2));
		sw.WriteLine(String.format("FileProperties\tGroupName:%1$s", FileProps.GroupName));
		sw.WriteLine(String.format("FileProperties\tClusterName:%1$s", FileProps.ClusterName));
		sw.WriteLine(String.format("FileProperties\tNumberOriginalPoints:%1$s", FileProps.NumberOriginalPoints));
		sw.WriteLine(String.format("FileProperties\tNumberPointsinFile:%1$s", FileProps.NumberPointsinFile));
		sw.WriteLine(String.format("FileProperties\tRotationParameters:%1$s", FileProps.RotationParameters));
		sw.WriteLine(String.format("FileProperties\tNumberRotationParameters:%1$s", FileProps.NumberRotationParameters));

		String FileComments = FileProps.Comment;
		if (!FileComments.equals(""))
		{
			String[] split = FileComments.split("[\\n]", -1);
			for (int runthroughcomments = 0; runthroughcomments < split.length; runthroughcomments++)
			{
				if (!split[runthroughcomments].equals(""))
				{
					sw.WriteLine(String.format("FileProperties\tComment:%1$s", split[runthroughcomments]));
				}
			}
		}
		return;

	} // End WriteFileProperties

	public static void ReadDataPointProperties(SALSADataPointProperties DataPointProps, String InputLine)
	{
		InputLine = tangible.DotNetToJavaStringHelper.trim(InputLine, new char[] {' ', '\t'});
		String[] split = InputLine.split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
		DataPointProps.valuesset = false;
		DataPointProps.errorsset = false;

		for (int itembackwards = 0; itembackwards < split.length; itembackwards++)
		{
			int item = split.length - itembackwards - 1;
			if (split[item].equals("PointProperties"))
			{
				continue;
			}
			String[] colonsplit = split[item].split(new char[] {':'}, 2);
			if (colonsplit[0].equals("PointProperties"))
			{
				continue;
			}
			if (colonsplit.length < 2)
			{
				if (item == 0)
				{
					System.out.println("No deliminator in item " + split[item] + " in Line " + InputLine);
				}
				else
				{
					split[item - 1] += " " + split[item];
					continue;
				}
			}
			colonsplit[0] = colonsplit[0].trim();
			colonsplit[1] = colonsplit[1].trim();
			if (colonsplit[0].equals("x"))
			{
				DataPointProps.x = Double.parseDouble(colonsplit[1]);
				DataPointProps.valuesset = true;
			}
			if (colonsplit[0].equals("y"))
			{
				DataPointProps.y = Double.parseDouble(colonsplit[1]);
			}
			if (colonsplit[0].equals("z"))
			{
				DataPointProps.z = Double.parseDouble(colonsplit[1]);
			}
			if (colonsplit[0].equals("xerr"))
			{
				DataPointProps.xerr = Double.parseDouble(colonsplit[1]);
				DataPointProps.errorsset = true;
			}
			if (colonsplit[0].equals("source"))
			{
				DataPointProps.source = colonsplit[1];
			}
			if (colonsplit[0].equals("yerr"))
			{
				DataPointProps.yerr = Double.parseDouble(colonsplit[1]);
			}
			if (colonsplit[0].equals("zerr"))
			{
				DataPointProps.zerr = Double.parseDouble(colonsplit[1]);
			}
			if (colonsplit[0].equals("family1"))
			{
				DataPointProps.family1 = Integer.parseInt(colonsplit[1]);
			}
			if (colonsplit[0].equals("familylabel1"))
			{
				DataPointProps.familylabel1 = colonsplit[1];
			}
			if (colonsplit[0].equals("family2"))
			{
				DataPointProps.family2 = Integer.parseInt(colonsplit[1]);
			}
			if (colonsplit[0].equals("familylabel2"))
			{
				DataPointProps.familylabel2 = colonsplit[1];
			}
			if (colonsplit[0].equals("cluster"))
			{
				DataPointProps.cluster = Integer.parseInt(colonsplit[1]);
			}
			if (colonsplit[0].equals("clusterlabel"))
			{
				DataPointProps.clusterlabel = colonsplit[1];
			}
			if (colonsplit[0].equals("group"))
			{
				DataPointProps.group = Integer.parseInt(colonsplit[1]);
			}
			if (colonsplit[0].equals("grouplabel"))
			{
				DataPointProps.grouplabel = colonsplit[1];
			}
			if (colonsplit[0].equals("pointlabel"))
			{
				DataPointProps.pointlabel = colonsplit[1];
			}
			if (colonsplit[0].equals("FixedorVaried"))
			{
				DataPointProps.FixedorVaried = Integer.parseInt(colonsplit[1]);
			}
			if (colonsplit[0].equals("PointType"))
			{
				DataPointProps.PointType = Integer.parseInt(colonsplit[1]);
			}
			if (colonsplit[0].equals("LocalPointNumber"))
			{
				DataPointProps.LocalPointNumber = Integer.parseInt(colonsplit[1]);
			}
			if (colonsplit[0].equals("OriginalPointNumber"))
			{
				DataPointProps.OriginalPointNumber = Integer.parseInt(colonsplit[1]);
			}
		}
		return;

	} // End ReadDataPointProperties

	public static void AppendDataPointProperties(SALSADataPointProperties DataPointProps, tangible.RefObject<String> InputLine)
	{ //Calling Routine must supply any needed delimiter before this is appended
		InputLine.argValue += "PointProperties";
		if (DataPointProps.valuesset == true)
		{ // x y z are set

			InputLine.argValue += "\tvaluesset:true\tsource:" + DataPointProps.source + "\tx:" + String.format("%.4E", DataPointProps.x);
			InputLine.argValue += "\ty:" + String.format("%.4E", DataPointProps.y);
			InputLine.argValue += "\tz:" + String.format("%.4E", DataPointProps.z);
		}
		if (DataPointProps.errorsset == true)
		{ // x y z errors are set
			InputLine.argValue += "\terrorsset:true\txerr:" + String.format("%.4E", DataPointProps.xerr);
			InputLine.argValue += "\tyerr:" + String.format("%.4E", DataPointProps.yerr);
			InputLine.argValue += "\tzerr:" + String.format("%.4E", DataPointProps.zerr);
		}
		if (DataPointProps.family1 >= 0)
		{
			InputLine.argValue += "\tfamily1:" + (new Integer(DataPointProps.family1)).toString();
			if (DataPointProps.familylabel1.length() > 0)
			{
				InputLine.argValue += "\tfamilylabel1:" + DataPointProps.familylabel1;
			}
		}
		if (DataPointProps.family2 >= 0)
		{
			InputLine.argValue += "\tfamily2:" + (new Integer(DataPointProps.family2)).toString();
			if (DataPointProps.familylabel2.length() > 0)
			{
				InputLine.argValue += "\tfamilylabel2:" + DataPointProps.familylabel2;
			}
		}
		if (DataPointProps.cluster >= 0)
		{
			InputLine.argValue += "\tcluster:" + (new Integer(DataPointProps.cluster)).toString();
			if (DataPointProps.clusterlabel.length() > 0)
			{
				InputLine.argValue += "\tclusterlabel:" + DataPointProps.clusterlabel;
			}
		}
		if (DataPointProps.group >= 0)
		{
			InputLine.argValue += "\tgroup:" + (new Integer(DataPointProps.group)).toString();
			if (DataPointProps.grouplabel.length() > 0)
			{
				InputLine.argValue += "\tgrouplabel:" + DataPointProps.grouplabel;
			}
		}
		if (DataPointProps.pointlabel.length() > 0)
		{
			InputLine.argValue += "\tpointlabel:" + DataPointProps.pointlabel;
		}
		InputLine.argValue += "\tFixedorVaried:" + (new Integer(DataPointProps.FixedorVaried)).toString();
		InputLine.argValue += "\tPointType:" + (new Integer(DataPointProps.PointType)).toString();
		if (DataPointProps.LocalPointNumber >= 0)
		{
			InputLine.argValue += "\tLocalPointNumber:" + (new Integer(DataPointProps.LocalPointNumber)).toString();
		}
		if (DataPointProps.OriginalPointNumber >= 0)
		{
			InputLine.argValue += "\tOriginalPointNumber:" + (new Integer(DataPointProps.OriginalPointNumber)).toString();
		}

	} // End AppendDataPointProperties

} // End SALSA_Properties