package salsa.mdsaschisq;

import java.io.File;
import java.util.regex.Pattern;

public class SALSA_ProcessVariedandFixed
{
	public static SALSADataPointProperties[] FixedPointProperties;
	public static SALSADataPointProperties[] VariedPointProperties;

	public static void setupFixedandVaried()
	{ // Set up Load Balanced Fixed and Varied Points
		// Assumes NumberOriginalPoints is original number and some ignored, some fixed and some varied
		// Set SALSAUtility.PointCount_Global to be number of fixed PLUS number of varied points which is number of Used points
		// This can be less than  SALSAUtility.NumberOriginalPoints as some points can be not used
		//
		// Manxcat already supports fixed parameters (for 6 fixed points to remove ambiguities in 3D rotation)
		// Use this but fix points specified by user instead of default fixing
		// Before fit starts, these points are rearranged so each thread has equal numbers (up to 1) of varied and fixed points
		//  On input distance matrix is rearranged so parameters are in order of points in thread (varied, fixed), threads in process, process
		// Distances involving not used points are dropped on input

		// Set Initial point and file Properties
		SALSAUtility.OriginalPointProperties = new SALSADataPointProperties[SALSAUtility.NumberOriginalPoints];
		SALSAUtility.OriginalFileProperties = new SALSAFileProperties();

		for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
		{
			SALSAUtility.OriginalPointProperties[OriginalPointIndex] = new SALSADataPointProperties();
		}

		// Set Properties of Points -- if properties file exists it MUST correspond to NumberOriginalPoints
		// If no file specified, default labels are generated
		//  Also initialize File and Point properties 
		String OriginalPropertiesfilename = ManxcatCentral.config.DataLabelsFileName;
		if (!OriginalPropertiesfilename.contains(":"))
		{
			OriginalPropertiesfilename = ManxcatCentral.config.ControlDirectoryName + File.separatorChar + OriginalPropertiesfilename;
		}
		SALSAUtility.OriginalDataLabels = new String[SALSAUtility.NumberOriginalPoints];
		int OriginalPropertiesFileType = -1;
		int InputNumberofPoints = 0;
		if ((new java.io.File(OriginalPropertiesfilename)).isFile())
		{
			tangible.RefObject<Integer> tempRef_OriginalPropertiesFileType = new tangible.RefObject<Integer>(OriginalPropertiesFileType);
			tangible.RefObject<SALSADataPointProperties[]> tempRef_OriginalPointProperties = new tangible.RefObject<SALSADataPointProperties[]>(SALSAUtility.OriginalPointProperties);
			tangible.RefObject<Integer> tempRef_InputNumberofPoints = new tangible.RefObject<Integer>(InputNumberofPoints);
			SALSA_Properties.ReadDataPointFile(OriginalPropertiesfilename, tempRef_OriginalPropertiesFileType, SALSAUtility.OriginalFileProperties, tempRef_OriginalPointProperties, tempRef_InputNumberofPoints);
			OriginalPropertiesFileType = tempRef_OriginalPropertiesFileType.argValue;
			SALSAUtility.OriginalPointProperties = tempRef_OriginalPointProperties.argValue;
			InputNumberofPoints = tempRef_InputNumberofPoints.argValue;
			if (SALSAUtility.NumberOriginalPoints != InputNumberofPoints)
			{
				SALSAUtility.printAndThrowRuntimeException(
                        " Inconsistent Point Counts " + InputNumberofPoints + " Expected is " + SALSAUtility
                                .NumberOriginalPoints);

			}
			for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
			{
				SALSAUtility.OriginalDataLabels[OriginalPointIndex] = SALSAUtility.OriginalPointProperties
                        [OriginalPointIndex].grouplabel + ":" + SALSAUtility
                        .OriginalPointProperties[OriginalPointIndex].group;
			}
		}
		else
		{
			SetDataLabels(SALSAUtility.OriginalDataLabels);
		}

		// Set File properties of Used data array identical to Original data array
		SALSAUtility.GlobalFileProperties = new SALSAFileProperties();
		SALSAUtility.GlobalFileProperties = SALSAUtility.OriginalFileProperties.ShallowCopy();
		// End setting of labels and properties

		// Arrays Mapping "Original" "Used" "Fixed" and "Varied" points
		// Used points are sum of fixed and varied and their total stored in SALSAUtility.PointCount_Global
		SALSAUtility.NumberFixedPoints = 0; // Number of Fixed Points
		SALSAUtility.NumberVariedPoints = 0; // Number of Varied Points
		SALSAUtility.OriginalPointDisposition = new int[SALSAUtility.NumberOriginalPoints]; // = 0 Original Point Not Used; = SALSASHIFT + n This is Varied Point n; = -SALSASHIFT -m This is Fixed Point m
		SALSAUtility.OriginalPointtoUsedPointMap = new int[SALSAUtility.NumberOriginalPoints]; // Value of Used Point Index Corresponding to this Original point, if -1 Not a used point

		// Initialize Original Point Values and their mapping to used
		for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
		{
			int testvalue = SALSAUtility.OriginalPointProperties[OriginalPointIndex].OriginalPointNumber;
			if ((testvalue != -1) && (testvalue != OriginalPointIndex))
			{
				SALSAUtility.printAndThrowRuntimeException(
                        " Position " + OriginalPointIndex + " Inconsistent Point Counts " + testvalue + " Expected is" +
                                " " + OriginalPointIndex);

			}
			SALSAUtility.OriginalPointProperties[OriginalPointIndex].OriginalPointNumber = OriginalPointIndex;
			SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex] = -1;
		}

		String VariedCollection = ManxcatCentral.config.VariedPointCriterion;
		String FixedCollection = ManxcatCentral.config.FixedPointCriterion;
		if (VariedCollection.length() == 0 || FixedCollection.length() == 0)
		{
			VariedCollection = "all";
			FixedCollection = "none";
		}

//  Process special case where all varied -- this is default scenario
		if ((ManxcatCentral.config.FixedPointCriterion.equals("none")) && (ManxcatCentral.config.VariedPointCriterion.equals("all")))
		{ // Nothing Special -- MDS on all points
			// SALSAUtility.FixedPointOriginal does not need to be defined as zero fixed points
			SALSAUtility.VariedPointOriginal = new int[SALSAUtility.NumberOriginalPoints]; // Value of Original Point Index Corresponding to this varied point
			SALSAUtility.UsedPointtoOriginalPointMap = new int[SALSAUtility.NumberOriginalPoints]; // Value of Original Point Index Corresponding to this Used point

			SALSAUtility.GlobalPointProperties = new SALSADataPointProperties[SALSAUtility.NumberOriginalPoints];
			SALSAUtility.NumberVariedPoints = SALSAUtility.NumberOriginalPoints;
			SALSAUtility.ActualtoNaiveUsedOrder = new int[SALSAUtility.NumberOriginalPoints];
			SALSAUtility.NaivetoActualUsedOrder = new int[SALSAUtility.NumberOriginalPoints];
			for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.NumberOriginalPoints; GlobalPointIndex++)
			{
				SALSAUtility.OriginalPointDisposition[GlobalPointIndex] = SALSAUtility.SALSASHIFT + GlobalPointIndex;
				SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex] = GlobalPointIndex;
				SALSAUtility.OriginalPointtoUsedPointMap[GlobalPointIndex] = GlobalPointIndex;
				SALSAUtility.VariedPointOriginal[GlobalPointIndex] = GlobalPointIndex;
				SALSAUtility.ActualtoNaiveUsedOrder[GlobalPointIndex] = GlobalPointIndex;
				SALSAUtility.NaivetoActualUsedOrder[GlobalPointIndex] = GlobalPointIndex;
				SALSAUtility.GlobalPointProperties[GlobalPointIndex] = SALSAUtility.OriginalPointProperties[GlobalPointIndex].ShallowCopy();
			}

			SALSAUtility.PointCount_Global = SALSAUtility.NumberOriginalPoints;
			SALSAUtility.NumberVariedPoints = SALSAUtility.PointCount_Global;
			if (SALSAUtility.MPI_Rank == 0)
			{
				SALSAUtility.SALSAPrint(0, " ALL points varied and none fixed");
			}
			return;
		}

		// Process Varied and Fixed Points in case where point selection non trivial
		for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
		{
			SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = 0;
			SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex] = -1;
		}

		//  Process Varied Points in case that varied selection specific -- cases of VariedPointCriterion = all or rest done later
		//  We need to find if points specified in file or in input string and if specification is a point number or a collection (cluster, family ..) number
		// VariedPointCriterion can be
		// all -- All points are varied (all treated like rest if FixedPointCriterion is not none)
		// rest -- All non fixed points are varied
		// allinfile -- vary all points specified in Selectedvariedpointfile
		// originalpoint -- vary all points in Selectedvariedpoints
		// family1, family2, group, cluster -- vary all points in collection defined in Selectedvariedpoints
		if (!VariedCollection.equals("all") && !VariedCollection.equals("rest"))
		{
			// The case when varied collection is not all or rest
			// Interpret the list of selected points
			String ToDecode = ManxcatCentral.config.Selectedvariedpoints.trim();
			if (SALSAUtility.MPI_Rank == 0)
			{
				SALSAUtility.SALSAPrint(0, " Varied Point Selection Categories: " + ToDecode + " Method: " + VariedCollection + " File " + ManxcatCentral.config.Selectedvariedpointfile);
			}

			if ((ToDecode.length() == 0) && (!VariedCollection.equals("allinfile")))
			{
				SALSAUtility.printAndThrowRuntimeException(" No Varied Points Selected");

			}

			int CollectionMax = 0;

			// Read Label File Containing extra information -- only property referenced by VariedPointCriterion is used
			// Error if no label file unless VariedCollection is originalpoint
			String VariedFilename = ManxcatCentral.config.Selectedvariedpointfile;
			if (VariedFilename.length() <= 0)
			{
				if (!VariedCollection.equals("originalpoint"))
				{
					SALSAUtility.printAndThrowRuntimeException(" No Varied File Specified");

				}
			}
			else if (!VariedCollection.equals("originalpoint"))
			{ // There is a VariedFilename and we need it - cases family1 family2 group cluster allinfile
				SALSAFileProperties VariedFileProperties = new SALSAFileProperties();
				if (!VariedFilename.contains(":"))
				{
					VariedFilename = ManxcatCentral.config.ControlDirectoryName + File.separatorChar + VariedFilename;
				}

				int VariedFileType = -1;
				int InputVariedNumberofPoints = 0;
				tangible.RefObject<Integer> tempRef_VariedFileType = new tangible.RefObject<Integer>(VariedFileType);
				tangible.RefObject<SALSADataPointProperties[]> tempRef_VariedPointProperties = new tangible.RefObject<SALSADataPointProperties[]>(VariedPointProperties);
				tangible.RefObject<Integer> tempRef_InputVariedNumberofPoints = new tangible.RefObject<Integer>(InputVariedNumberofPoints);
				SALSA_Properties.ReadDataPointFile(VariedFilename, tempRef_VariedFileType, VariedFileProperties, tempRef_VariedPointProperties, tempRef_InputVariedNumberofPoints);
				VariedFileType = tempRef_VariedFileType.argValue;
				VariedPointProperties = tempRef_VariedPointProperties.argValue;
				InputVariedNumberofPoints = tempRef_InputVariedNumberofPoints.argValue;
				if ((SALSAUtility.NumberOriginalPoints < InputVariedNumberofPoints) || (VariedFileProperties.NumberOriginalPoints != SALSAUtility.NumberOriginalPoints))
				{
					SALSAUtility.printAndThrowRuntimeException(
                            " Inconsistent Varied File Point Counts " + InputVariedNumberofPoints + " or " +
                                    VariedFileProperties.NumberOriginalPoints + " Expected is " + SALSAUtility
                                    .NumberOriginalPoints);

				}
				if (VariedCollection.equals("family1") || VariedCollection.equals("family2") || VariedCollection.equals("group") || VariedCollection.equals("cluster"))
				{
					if (SALSAUtility.NumberOriginalPoints != InputVariedNumberofPoints)
					{
						SALSAUtility.printAndThrowRuntimeException(
                                " Inconsistent Varied File Point Counts " + InputVariedNumberofPoints + " Expected is" +
                                        " " + SALSAUtility.NumberOriginalPoints);

					}
				}

				SetLabelsfromFile(VariedCollection, VariedFileProperties, SALSAUtility.GlobalFileProperties);
				tangible.RefObject<Integer> tempRef_CollectionMax = new tangible.RefObject<Integer>(CollectionMax);
				SetPointsfromFile(tempRef_CollectionMax, 1, VariedCollection, InputVariedNumberofPoints, VariedPointProperties);
				CollectionMax = tempRef_CollectionMax.argValue;

			} // End reading information from Varied Selection File

			// Above call to SetPointsfromFile finishes Varied processing for allinfile case. 
			//  For originalpoint family1 family2 cluster group still need to process
			if (!VariedCollection.equals("allinfile"))
			{
				int NumberinVariedPointList = 0;
				if (VariedCollection.equals("originalpoint"))
				{
					CollectionMax = SALSAUtility.NumberOriginalPoints - 1;
				}
				int[] VariedPointList = new int[1 + CollectionMax];

				// This sets SALSAUtility.OriginalPointDisposition if originalpoint
				tangible.RefObject<Integer> tempRef_NumberinVariedPointList = new tangible.RefObject<Integer>(NumberinVariedPointList);
				InterpretIntegerSelection(ToDecode, VariedCollection, 1, tempRef_NumberinVariedPointList, VariedPointList, CollectionMax);
				NumberinVariedPointList = tempRef_NumberinVariedPointList.argValue;

				if (VariedCollection.equals("originalpoint"))
				{
					for (int SelectIndex = 0; SelectIndex < NumberinVariedPointList; SelectIndex++)
					{
						int OriginalPointIndex = VariedPointList[SelectIndex];
						SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = SALSAUtility.SALSASHIFT + SelectIndex;
					}
				}
				else
				{
					int LocalVariedCount = 0;
					for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
					{
						int TestIndex = -1;
						if (VariedCollection.equals("family1"))
						{
							TestIndex = VariedPointProperties[OriginalPointIndex].family1;
						}
						else if (VariedCollection.equals("family2"))
						{
							TestIndex = VariedPointProperties[OriginalPointIndex].family2;
						}
						else if (VariedCollection.equals("group"))
						{
							TestIndex = VariedPointProperties[OriginalPointIndex].group;
						}
						else if (VariedCollection.equals("cluster"))
						{
							TestIndex = VariedPointProperties[OriginalPointIndex].cluster;
						}
						for (int SelectIndex = 0; SelectIndex < NumberinVariedPointList; SelectIndex++)
						{
							if (TestIndex != VariedPointList[SelectIndex])
							{
								continue;
							}
							SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = SALSAUtility.SALSASHIFT + LocalVariedCount;
							++LocalVariedCount;
							break;
						}
					} // End Loop over GlobalPoinTindex Setting Varied Character

				} // End cases where some sort of collection specifies Varied Character
			} // End case when VariedCollection is not all in File

		} // End Process Varied Points in cases other that "all" and "rest"
		else
		{
			if (SALSAUtility.MPI_Rank == 0)
			{
				SALSAUtility.SALSAPrint(0, " Varied Point Selection Method: " + VariedCollection);
			}
		} // If VariedCollection all or rest just print out as it will be processed after fixed points

		//  Process Fixed Points
		if (FixedCollection.equals("all") || FixedCollection.equals("rest"))
		{ // Set all the non varied points to be fixed
			int LocalFixedCount = 0;
			for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
			{
				if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex] >= SALSAUtility.SALSASHIFT)
				{
					continue;
				}
				SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = - SALSAUtility.SALSASHIFT - LocalFixedCount;
				++LocalFixedCount;
			}
		}
		else if (!FixedCollection.equals("none"))
		{ // Cases originalpoint family1 family2 cluster group allinfile
			String ToDecode = ManxcatCentral.config.Selectedfixedpoints.trim();
			if (SALSAUtility.MPI_Rank == 0)
			{
				SALSAUtility.SALSAPrint(0, " Fixed Point Selection Categories: " + ToDecode + " Method: " + FixedCollection + " File " + ManxcatCentral.config.Selectedfixedpointfile);
			}
			if ((ToDecode.length() == 0) && (!FixedCollection.equals("allinfile")))
			{
				SALSAUtility.printAndThrowRuntimeException(" No Fixed Points Selected");

			}

			// Read Label or Properties File Containing extra information -- only property referenced by FixedPointCriterion is used
			int CollectionMax = 0;
			String FixedFilename = ManxcatCentral.config.Selectedfixedpointfile;

			if (FixedFilename.length() <= 0)
			{ // Error if no file unless original point option
				if (!FixedCollection.equals("originalPoint"))
				{
					SALSAUtility.printAndThrowRuntimeException(" No Fixed File Specified");

				}
			}
			else
			{ // Process case that file exists. For case of original point file, file read for fixed positions; for other cases read for collection numbers and fixed positions
				if (!FixedFilename.contains(":") && !FixedFilename.contains("$"))
				{
					FixedFilename = ManxcatCentral.config.ControlDirectoryName + File.separatorChar + FixedFilename;
				}

				// There is a FixedFilename 
				int FixedFileType = -1;
				int InputFixedNumberofPoints = 0;
				SALSAFileProperties FixedFileProperties = new SALSAFileProperties();
				tangible.RefObject<Integer> tempRef_FixedFileType = new tangible.RefObject<Integer>(FixedFileType);
				tangible.RefObject<SALSADataPointProperties[]> tempRef_FixedPointProperties = new tangible.RefObject<SALSADataPointProperties[]>(FixedPointProperties);
				tangible.RefObject<Integer> tempRef_InputFixedNumberofPoints = new tangible.RefObject<Integer>(InputFixedNumberofPoints);
				SALSA_Properties.ReadDataPointFile(FixedFilename, tempRef_FixedFileType, FixedFileProperties, tempRef_FixedPointProperties, tempRef_InputFixedNumberofPoints);
				FixedFileType = tempRef_FixedFileType.argValue;
				FixedPointProperties = tempRef_FixedPointProperties.argValue;
				InputFixedNumberofPoints = tempRef_InputFixedNumberofPoints.argValue;

				if ((SALSAUtility.NumberOriginalPoints < InputFixedNumberofPoints) || (FixedFileProperties.NumberOriginalPoints != SALSAUtility.NumberOriginalPoints))
				{
					SALSAUtility.printAndThrowRuntimeException(
                            " Inconsistent Fixed File Point Counts " + InputFixedNumberofPoints + " or " +
                                    FixedFileProperties.NumberOriginalPoints + " Expected is "
                                    + SALSAUtility.NumberOriginalPoints);

				}
				if (FixedCollection.equals("family1") || FixedCollection.equals("family2") || FixedCollection.equals("group") || FixedCollection.equals("cluster"))
				{
					if (SALSAUtility.NumberOriginalPoints != InputFixedNumberofPoints)
					{
						SALSAUtility.printAndThrowRuntimeException(
                                " Inconsistent Fixed File Point Counts " + InputFixedNumberofPoints + " Expected is "
                                        + SALSAUtility.NumberOriginalPoints);

					}
				}

				SetLabelsfromFile(FixedCollection, FixedFileProperties, SALSAUtility.GlobalFileProperties);
				tangible.RefObject<Integer> tempRef_CollectionMax2 = new tangible.RefObject<Integer>(CollectionMax);
				SetPointsfromFile(tempRef_CollectionMax2, -1, FixedCollection, InputFixedNumberofPoints, FixedPointProperties);
				CollectionMax = tempRef_CollectionMax2.argValue;
			} // End processing Selectedfixedpointfile file

			// Above call to SetPointsfromFile finishes Fixed processing for allinfile case. 
			//  For originalpoint family1 family2 cluster group still need to process
			if (!FixedCollection.equals("allinfile"))
			{
				int NumberinFixedPointList = 0;
				if (FixedCollection.equals("originalpoint"))
				{
					CollectionMax = SALSAUtility.NumberOriginalPoints - 1;
				}
				int[] FixedPointList = new int[1 + CollectionMax];

				tangible.RefObject<Integer> tempRef_NumberinFixedPointList = new tangible.RefObject<Integer>(NumberinFixedPointList);
				InterpretIntegerSelection(ToDecode, FixedCollection, -1, tempRef_NumberinFixedPointList, FixedPointList, CollectionMax);
				NumberinFixedPointList = tempRef_NumberinFixedPointList.argValue;

				if (FixedCollection.equals("originalpoint"))
				{
					for (int SelectIndex = 0; SelectIndex < NumberinFixedPointList; SelectIndex++)
					{
						int OriginalPointIndex = FixedPointList[SelectIndex];
						SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = -SALSAUtility.SALSASHIFT - SelectIndex;
					}
				}
				else
				{ // Cases family1 family2 group cluster
					int LocalFixedCount = 0;
					for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
					{
						int TestIndex = -1;
						if (FixedCollection.equals("family1"))
						{
							TestIndex = FixedPointProperties[OriginalPointIndex].family1;
						}
						else if (FixedCollection.equals("family2"))
						{
							TestIndex = FixedPointProperties[OriginalPointIndex].family2;
						}
						else if (FixedCollection.equals("group"))
						{
							TestIndex = FixedPointProperties[OriginalPointIndex].group;
						}
						else if (FixedCollection.equals("cluster"))
						{
							TestIndex = FixedPointProperties[OriginalPointIndex].cluster;
						}
						for (int SelectIndex = 0; SelectIndex < NumberinFixedPointList; SelectIndex++)
						{
							if (TestIndex != FixedPointList[SelectIndex])
							{
								continue;
							}
							SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = -SALSAUtility.SALSASHIFT - LocalFixedCount;
							++LocalFixedCount;
							break;
						}
					} // End Loop over GlobalPoinTindex Setting Fixed Character

				} // End cases where some sort of collection specifies Fixed Character
			}

		} // End Process Fixed Points if they exist
		else
		{ // Case of No Fixed Points
			if (SALSAUtility.MPI_Rank == 0)
			{
				SALSAUtility.SALSAPrint(0, "Fixed Point Selection Method: " + VariedCollection);
			}
		}

		//  Do Special Cases of Varied
		if (VariedCollection.equals("all") || VariedCollection.equals("rest"))
		{
			int LocalVariedCount = 0;
			for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
			{
				if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex] > - SALSAUtility.SALSASHIFT)
				{
					SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = SALSAUtility.SALSASHIFT + LocalVariedCount;
					++LocalVariedCount;
				}
			}
			SALSAUtility.NumberVariedPoints = LocalVariedCount;
			SALSAUtility.NumberFixedPoints = SALSAUtility.NumberOriginalPoints - SALSAUtility.NumberVariedPoints;
		}

		//  Ensure Fixed and Varied points assigned in order
		//  Set SALSAUtility.OriginalPointDisposition
		int VariedCount = 0;
		int FixedCount = 0;
		for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
		{
			if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex] <= -SALSAUtility.SALSASHIFT)
			{
				SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = -SALSAUtility.SALSASHIFT - FixedCount;
				++FixedCount;
			}

			else if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex] >= SALSAUtility.SALSASHIFT)
			{
				SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = SALSAUtility.SALSASHIFT + VariedCount;
				++VariedCount;
			}
			else
			{
				SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = 0;
			}
		}
		SALSAUtility.NumberVariedPoints = VariedCount;
		SALSAUtility.NumberFixedPoints = FixedCount;

		// Set arrays defining mapping between Original Fixed and Varied Points
		if (SALSAUtility.NumberFixedPoints > 0)
		{
			SALSAUtility.FixedPointOriginal = new int[SALSAUtility.NumberFixedPoints]; // Value of Original Point Index Corresponding to this fixed point
		}
		SALSAUtility.VariedPointOriginal = new int[SALSAUtility.NumberVariedPoints]; // Value of Original Point Index Corresponding to this varied point
		SALSAUtility.PointCount_Global = SALSAUtility.NumberFixedPoints + SALSAUtility.NumberVariedPoints;
		SALSAUtility.ActualtoNaiveUsedOrder = new int[SALSAUtility.PointCount_Global];
		SALSAUtility.NaivetoActualUsedOrder = new int[SALSAUtility.PointCount_Global];
		SALSAUtility.GlobalPointProperties = new SALSADataPointProperties[SALSAUtility.PointCount_Global];
		SALSAUtility.UsedPointtoOriginalPointMap = new int[SALSAUtility.PointCount_Global]; // Value of Original Point Index Corresponding to this Used point

		int CountUsed = 0;
		int CountFixed = 0;
		int CountVaried = 0;
		for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
		{
			if ((SALSAUtility.OriginalPointDisposition[OriginalPointIndex] > -SALSAUtility.SALSASHIFT) && (SALSAUtility.OriginalPointDisposition[OriginalPointIndex] < SALSAUtility.SALSASHIFT))
			{
				continue;
			}
			if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex] <= -SALSAUtility.SALSASHIFT)
			{
				SALSAUtility.FixedPointOriginal[CountFixed] = OriginalPointIndex;
				++CountFixed;
			}
			if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex] >= SALSAUtility.SALSASHIFT)
			{
				SALSAUtility.VariedPointOriginal[CountVaried] = OriginalPointIndex;
				++CountVaried;
			}
			SALSAUtility.UsedPointtoOriginalPointMap[CountUsed] = OriginalPointIndex;
			SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex] = CountUsed;
			SALSAUtility.GlobalPointProperties[CountUsed] = SALSAUtility.OriginalPointProperties[OriginalPointIndex].ShallowCopy();
			SALSAUtility.ActualtoNaiveUsedOrder[CountUsed] = CountUsed;
			SALSAUtility.NaivetoActualUsedOrder[CountUsed] = CountUsed;
			++CountUsed;
		}
		// Check Fixed Points do have positions set with values stored in GlobalPointProperties which would have been copied from OriginalPointProperties
		//  There are three possible input files read here: 
		//  DataLabelsFileName Selectedfixedpointfile Selectedvariedpointfile and last two overwrite first file and need not be full set
		int missingdata = 0;
		for (int FixedPointIndex = 0; FixedPointIndex < SALSAUtility.NumberFixedPoints; FixedPointIndex++)
		{
			int OriginalPointIndex = SALSAUtility.FixedPointOriginal[FixedPointIndex];
			int usedpoint = SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex];
			if (!SALSAUtility.GlobalPointProperties[usedpoint].valuesset)
			{
				++missingdata;
			}
		}

		if (SALSAUtility.MPI_Rank == 0)
		{
			if (missingdata > 0)
			{
				SALSAUtility.SALSAPrint(0,
                                        " At this stage missing initial vales for " + missingdata + " Points -- must be in initialization file");
			}

			int NotUsed = SALSAUtility.NumberOriginalPoints - SALSAUtility.NumberVariedPoints - SALSAUtility.NumberFixedPoints;
			SALSAUtility.SALSAPrint(0,
                                    " Number of Varied Points " + SALSAUtility.NumberVariedPoints + " Fixed Points "
                                            + SALSAUtility.NumberFixedPoints + " Not Used " + NotUsed);
		}

    } // End setupFixedandVaried()

	public static void RedistributePoints()
	{ //  Redistribute Points if needed for Load Balancing

		if (SALSAUtility.NumberVariedPoints == SALSAUtility.NumberOriginalPoints)
		{
			return;
		}

		int VariedperProcess = SALSAUtility.NumberVariedPoints / SALSAUtility.MPI_Size;
		int VariedperProcessleft = SALSAUtility.NumberVariedPoints - VariedperProcess * SALSAUtility.MPI_Size;
		int FixedperProcess = SALSAUtility.NumberFixedPoints / SALSAUtility.MPI_Size;

		int VariedPositioninGlobal = 0;
		int FixedPositioninGlobal = 0;
		int CountUsed = 0;
		int CountVaried = 0;
		int CountFixed = 0;

		for (int ProcessIndex = 0; ProcessIndex < SALSAUtility.MPI_Size; ProcessIndex++)
		{
			int LocalProcessCount = SALSAUtility.PointsperProcess[ProcessIndex];

			// Divide Process Points into three Categories
			int VariedthisProcess = VariedperProcess;
			int Processleft = LocalProcessCount - VariedthisProcess - FixedperProcess;

			if ((Processleft > 2) || (Processleft < 0))
			{
				SALSAUtility.printAndThrowRuntimeException(
                        "Internal Error in LoadBalancing I -- Complain to Programmer MPI " + SALSAUtility.MPI_Rank +
                                " Process " + LocalProcessCount);

			}
			if (Processleft >= 1)
			{
				if (VariedperProcessleft > 0)
				{
					++VariedthisProcess;
					--VariedperProcessleft;
				}
			}
			int FixedthisProcess = LocalProcessCount - VariedthisProcess;
			if (ProcessIndex == SALSAUtility.MPI_Rank)
			{
				SALSAUtility.VariedPointCount_Process = VariedthisProcess;
				SALSAUtility.VariedPointStart_Process = CountVaried;
			}

			int VariedperThread = VariedthisProcess / SALSAUtility.ThreadCount;
			int VariedperThreadleft = VariedthisProcess - VariedperThread * SALSAUtility.ThreadCount;
			int FixedperThread = FixedthisProcess / SALSAUtility.ThreadCount;

			//	Now divide points among threads
			for (int ThreadNo = 0; ThreadNo < SALSAUtility.ThreadCount; ThreadNo++)
			{
				int LocalNumberInThread = SALSAUtility.PointsperThreadperProcess[ProcessIndex][ThreadNo];

				// Divide Thread Points into three Categories
				int VariedthisThread = VariedperThread;
				int Threadleft = LocalNumberInThread - VariedthisThread - FixedperThread;

				if ((Threadleft > 2) || (Threadleft < 0))
				{
					SALSAUtility.printAndThrowRuntimeException(
                            "Internal Error in LoadBalancing II -- Complain to Programmer MPI " + SALSAUtility
                                    .MPI_Rank + " Thread " + ThreadNo + " Process " + ProcessIndex + " # in Thread "
                                    + LocalNumberInThread
                                    + " # in Process " + LocalProcessCount);

				}
				if (Threadleft >= 1)
				{
					if (VariedperThreadleft > 0)
					{
						++VariedthisThread;
						--VariedperThreadleft;
					}
				}
				int FixedthisThread = LocalNumberInThread - VariedthisThread;
				for (int LocalVaried = 0; LocalVaried < VariedthisThread; LocalVaried++)
				{
					while (true)
					{
						if (SALSAUtility.OriginalPointDisposition[VariedPositioninGlobal] >= SALSAUtility.SALSASHIFT)
						{
							SALSAUtility.UsedPointtoOriginalPointMap[CountUsed] = VariedPositioninGlobal;
							if (CountVaried != (SALSAUtility.OriginalPointDisposition[VariedPositioninGlobal] - SALSAUtility.SALSASHIFT))
							{
								SALSAUtility.printAndThrowRuntimeException(
                                        "Error in Varied Point Assignment II -- Complain to Programmer MPI " +
                                                SALSAUtility.MPI_Rank + " Thread " + ThreadNo + " Varied " +
                                                VariedPositioninGlobal);

							}
							SALSAUtility.VariedPointOriginal[CountVaried] = VariedPositioninGlobal;
							++VariedPositioninGlobal;
							++CountUsed;
							++CountVaried;
							break;
						}
						++VariedPositioninGlobal;
						if (VariedPositioninGlobal >= SALSAUtility.NumberOriginalPoints)
						{
							SALSAUtility.printAndThrowRuntimeException(
                                    "Error in Varied Point Assignment I -- Complain to Programmer MPI " +
                                            SALSAUtility.MPI_Rank + " Thread " + ThreadNo + " Varied " +
                                            VariedPositioninGlobal);

						}
					}

				}
				for (int LocalFixed = 0; LocalFixed < FixedthisThread; LocalFixed++)
				{
					while (true)
					{
						if (SALSAUtility.OriginalPointDisposition[FixedPositioninGlobal] <= -SALSAUtility.SALSASHIFT)
						{
							SALSAUtility.UsedPointtoOriginalPointMap[CountUsed] = FixedPositioninGlobal;
							if (CountFixed != (-SALSAUtility.OriginalPointDisposition[FixedPositioninGlobal] - SALSAUtility.SALSASHIFT))
							{
								SALSAUtility.printAndThrowRuntimeException(
                                        "Error in Fixed Point Assignment II -- Complain to Programmer MPI " +
                                                SALSAUtility.MPI_Rank + " Thread " + ThreadNo + " Fixed " +
                                                FixedPositioninGlobal);

							}
							SALSAUtility.FixedPointOriginal[CountFixed] = FixedPositioninGlobal;
							++CountFixed;
							++FixedPositioninGlobal;
							++CountUsed;
							break;
						}
						++FixedPositioninGlobal;
						if (FixedPositioninGlobal >= SALSAUtility.NumberOriginalPoints)
						{
							SALSAUtility.printAndThrowRuntimeException(
                                    "Error in Fixed Point Assignment I -- Complain to Programmer MPI " + SALSAUtility
                                            .MPI_Rank + " Thread " + ThreadNo + " Fixed " + FixedPositioninGlobal);

						}
					}

				}
			} // End loop assigning points to threads

		} // End loop over Assigning points to processes
		if (CountUsed != SALSAUtility.PointCount_Global)
		{
			SALSAUtility.printAndThrowRuntimeException(
                    "Error in Used Point Assignment -- Complain to Programmer " + SALSAUtility.MPI_Rank + " Count " +
                            "Used " + CountUsed);

		}

		for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
		{
			int OriginalPointIndex = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
			int oldused = SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex];
			SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex] = GlobalPointIndex;
			SALSAUtility.ActualtoNaiveUsedOrder[GlobalPointIndex] = oldused;
			if (GlobalPointIndex != oldused)
			{
				SALSAUtility.Usedreordered = true;
			}
			SALSAUtility.NaivetoActualUsedOrder[oldused] = GlobalPointIndex;
			SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex] = GlobalPointIndex;
			SALSAUtility.GlobalPointProperties[GlobalPointIndex] = SALSAUtility.OriginalPointProperties[OriginalPointIndex].ShallowCopy();
		}
		for (int OriginalPointIndex = 0; OriginalPointIndex < SALSAUtility.NumberOriginalPoints; OriginalPointIndex++)
		{
			if (SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex] >= 0)
			{
				continue;
			}
			if ((SALSAUtility.OriginalPointDisposition[OriginalPointIndex] <= -SALSAUtility.SALSASHIFT) || (SALSAUtility.OriginalPointDisposition[OriginalPointIndex] >= SALSAUtility.SALSASHIFT))
			{
				SALSAUtility.printAndThrowRuntimeException(
                        "Error in Redistribution -- Complain to Programmer " + SALSAUtility.MPI_Rank + " Original " +
                                "Point " + OriginalPointIndex);

			}
		}
	} // End ReDistributePoints

	//  VariedorFixedFileProperties is properties of local varied/fixed file
	//  UsedPointFileProperties is properties of Used Point File initialized to Original Point File Properties above
	//  UsedPointFileProperties is reset to reflect information from VariedorFixedFileProperties according to value of VariedorFixedCollectionName
	//  VariedorFixedCollectionName is family1, family2, group,  OR cluster
	public static void SetLabelsfromFile(String VariedorFixedCollectionName, SALSAFileProperties VariedorFixedFileProperties, SALSAFileProperties UsedPointFileProperties)
	{
		String newlabel = "";
		if (VariedorFixedCollectionName.equals("group"))
		{
			if ((VariedorFixedFileProperties.GroupName.length() > 0) && !(VariedorFixedFileProperties.GroupName.equals("Undefined Group")))
			{
				newlabel = VariedorFixedFileProperties.GroupName;
			}
			else
			{
				if ((UsedPointFileProperties.GroupName.length() > 0) && !(UsedPointFileProperties.GroupName.equals("Undefined Group")))
				{
					newlabel = UsedPointFileProperties.GroupName;
				}
			}
			if (newlabel.length() == 0)
			{
				newlabel = "Undefined Group";
			}
			UsedPointFileProperties.GroupName = newlabel;
			return;
		}
		if (VariedorFixedCollectionName.equals("cluster"))
		{
			if ((VariedorFixedFileProperties.ClusterName.length() > 0) && !(VariedorFixedFileProperties.ClusterName.equals("Undefined Cluster")))
			{
				newlabel = VariedorFixedFileProperties.ClusterName;
			}
			else
			{
				if ((UsedPointFileProperties.ClusterName.length() > 0) && !(UsedPointFileProperties.ClusterName.equals("Undefined Cluster")))
				{
					newlabel = UsedPointFileProperties.ClusterName;
				}
			}
			if (newlabel.length() == 0)
			{
				newlabel = "Undefined Cluster";
			}
			UsedPointFileProperties.ClusterName = newlabel;
			return;
		}
		if (VariedorFixedCollectionName.equals("family1"))
		{
			if ((VariedorFixedFileProperties.FamilyName1.length() > 0) && !(VariedorFixedFileProperties.FamilyName1.equals("Undefined Family1")))
			{
				newlabel = VariedorFixedFileProperties.FamilyName1;
			}
			else
			{
				if ((UsedPointFileProperties.FamilyName1.length() > 0) && !(UsedPointFileProperties.FamilyName1.equals("Undefined Family1")))
				{
					newlabel = UsedPointFileProperties.FamilyName1;
				}
			}
			if (newlabel.length() == 0)
			{
				newlabel = "Undefined Family1";
			}
			UsedPointFileProperties.FamilyName1 = newlabel;
			return;
		}
		if (VariedorFixedCollectionName.equals("family2"))
		{
			if ((VariedorFixedFileProperties.FamilyName2.length() > 0) && !(VariedorFixedFileProperties.FamilyName2.equals("Undefined Family2")))
			{
				newlabel = VariedorFixedFileProperties.FamilyName2;
			}
			else
			{
				if ((UsedPointFileProperties.FamilyName2.length() > 0) && !(UsedPointFileProperties.FamilyName2.equals("Undefined Family2")))
				{
					newlabel = UsedPointFileProperties.FamilyName2;
				}
			}
			if (newlabel.length() == 0)
			{
				newlabel = "Undefined Family2";
			}
			UsedPointFileProperties.FamilyName2 = newlabel;
        }

	} // End SetLabelsfromFile

	//  VariedFixedFudge = 1 if varied and -1 if Fixed
	//  VariedorFixedCollectionName is specification of Collection allinfile family1 family2 group cluster
	//  InputNumberofPoints is number of points in varied or fixed file
	//  VariedorFixedPointProperties is array of properties of local read in file
	//
	//  Return if allinfile, OriginalPointDisposition set appropriately
	//  If family1 family2 group cluster
	//  OriginalPointProperties set properly for these labels and CollectionMAX holding largest value
	//  Copy x y z xerr yerr zerr to OriginalPointProperties
	public static void SetPointsfromFile(tangible.RefObject<Integer> CollectionMax, int VariedFixedFudge, String VariedorFixedCollectionName, int InputNumberofPoints, SALSADataPointProperties[] VariedorFixedPointProperties)
	{
		CollectionMax.argValue = 0;
		for (int VariedorFixedPointIndex = 0; VariedorFixedPointIndex < InputNumberofPoints; VariedorFixedPointIndex++)
		{
			int OriginalPointIndex = VariedorFixedPointProperties[VariedorFixedPointIndex].OriginalPointNumber;
			if (OriginalPointIndex < 0)
			{
				if (InputNumberofPoints != SALSAUtility.NumberOriginalPoints)
				{
					String calltype = "Varied";
					if (VariedFixedFudge < 0)
					{
						calltype = "Fixed";
					}
					SALSAUtility.printAndThrowRuntimeException(
                            " Unspecified Original Point in " + calltype + " File " + VariedorFixedPointIndex);

				}
				OriginalPointIndex = VariedorFixedPointIndex;
				VariedorFixedPointProperties[VariedorFixedPointIndex].OriginalPointNumber = OriginalPointIndex;
			}

			// Copy Coordinates if unset in basic file
			if ((!SALSAUtility.OriginalPointProperties[OriginalPointIndex].valuesset) && VariedorFixedPointProperties[VariedorFixedPointIndex].valuesset)
			{
				SALSAUtility.OriginalPointProperties[OriginalPointIndex].x = VariedorFixedPointProperties[VariedorFixedPointIndex].x;
				SALSAUtility.OriginalPointProperties[OriginalPointIndex].y = VariedorFixedPointProperties[VariedorFixedPointIndex].y;
				SALSAUtility.OriginalPointProperties[OriginalPointIndex].z = VariedorFixedPointProperties[VariedorFixedPointIndex].z;
				SALSAUtility.OriginalPointProperties[OriginalPointIndex].valuesset = VariedorFixedPointProperties[VariedorFixedPointIndex].valuesset;

				SALSAUtility.OriginalPointProperties[OriginalPointIndex].xerr = VariedorFixedPointProperties[VariedorFixedPointIndex].xerr;
				SALSAUtility.OriginalPointProperties[OriginalPointIndex].yerr = VariedorFixedPointProperties[VariedorFixedPointIndex].yerr;
				SALSAUtility.OriginalPointProperties[OriginalPointIndex].zerr = VariedorFixedPointProperties[VariedorFixedPointIndex].zerr;
				SALSAUtility.OriginalPointProperties[OriginalPointIndex].errorsset = VariedorFixedPointProperties[VariedorFixedPointIndex].errorsset;
			}

			// Set status for allinfile case
			if (VariedorFixedCollectionName.equals("allinfile"))
			{
				SALSAUtility.OriginalPointDisposition[OriginalPointIndex] = VariedFixedFudge * (SALSAUtility.SALSASHIFT + VariedorFixedPointIndex);
			}

			// Set requisite Collection Type
			if (VariedorFixedCollectionName.equals("family1"))
			{
				int basic = SALSAUtility.OriginalPointProperties[OriginalPointIndex].family1;
				int replacement = VariedorFixedPointProperties[VariedorFixedPointIndex].family1;
				if (replacement >= 0)
				{
					SALSAUtility.OriginalPointProperties[OriginalPointIndex].family1 = replacement;
					SALSAUtility.OriginalPointProperties[OriginalPointIndex].familylabel1 = VariedorFixedPointProperties[VariedorFixedPointIndex].familylabel1;
				}
				else
				{
					if (basic < 0)
					{
						SALSAUtility.OriginalPointProperties[OriginalPointIndex].family1 = 0;
					}
				}
				CollectionMax.argValue = Math.max(CollectionMax.argValue, SALSAUtility.OriginalPointProperties[OriginalPointIndex].family1);
				return;
			}
			if (VariedorFixedCollectionName.equals("family2"))
			{
				int basic = SALSAUtility.OriginalPointProperties[OriginalPointIndex].family2;
				int replacement = VariedorFixedPointProperties[VariedorFixedPointIndex].family2;
				if (replacement >= 0)
				{
					SALSAUtility.OriginalPointProperties[OriginalPointIndex].family2 = replacement;
					SALSAUtility.OriginalPointProperties[OriginalPointIndex].familylabel2 = VariedorFixedPointProperties[VariedorFixedPointIndex].familylabel2;
				}
				else
				{
					if (basic < 0)
					{
						SALSAUtility.OriginalPointProperties[OriginalPointIndex].family2 = 0;
					}
				}
				CollectionMax.argValue = Math.max(CollectionMax.argValue, SALSAUtility.OriginalPointProperties[OriginalPointIndex].family2);
				return;
			}
			if (VariedorFixedCollectionName.equals("group"))
			{
				int basic = SALSAUtility.OriginalPointProperties[OriginalPointIndex].group;
				int replacement = VariedorFixedPointProperties[VariedorFixedPointIndex].group;
				if (replacement >= 0)
				{
					SALSAUtility.OriginalPointProperties[OriginalPointIndex].group = replacement;
					SALSAUtility.OriginalPointProperties[OriginalPointIndex].grouplabel = VariedorFixedPointProperties[VariedorFixedPointIndex].grouplabel;
				}
				else
				{
					if (basic < 0)
					{
						SALSAUtility.OriginalPointProperties[OriginalPointIndex].group = 0;
					}
				}
				CollectionMax.argValue = Math.max(CollectionMax.argValue, SALSAUtility.OriginalPointProperties[OriginalPointIndex].group);
				return;
			}
			if (VariedorFixedCollectionName.equals("cluster"))
			{
				int basic = SALSAUtility.OriginalPointProperties[OriginalPointIndex].cluster;
				int replacement = VariedorFixedPointProperties[VariedorFixedPointIndex].cluster;
				if (replacement >= 0)
				{
					SALSAUtility.OriginalPointProperties[OriginalPointIndex].cluster = replacement;
					SALSAUtility.OriginalPointProperties[OriginalPointIndex].clusterlabel = VariedorFixedPointProperties[VariedorFixedPointIndex].clusterlabel;
				}
				else
				{
					if (basic < 0)
					{
						SALSAUtility.OriginalPointProperties[OriginalPointIndex].cluster = 0;
					}
				}
				CollectionMax.argValue = Math.max(CollectionMax.argValue, SALSAUtility.OriginalPointProperties[OriginalPointIndex].cluster);
				return;
			}
		}

	} // end SetPointsfromFile

	// Cases family1 family2 cluster group originalpoint with maximum number of selected points CollectionMax+1
	//  VariedorFixedPointList holds selected points with NumberinVariedorFixedPointList members
	public static void InterpretIntegerSelection(String ToDecode, String CollectionString, int VariedFixedFudge, tangible.RefObject<Integer> NumberinVariedorFixedPointList, int[] VariedorFixedPointList, int CollectionMax)
	{
		String calltype = "Varied";
		if (VariedFixedFudge < 0)
		{
			calltype = "Fixed";
		}
		NumberinVariedorFixedPointList.argValue = 0;
        Pattern pattern = Pattern.compile("[,]");
		String[] ItemsinToDecode = pattern.split(ToDecode);
		if (ItemsinToDecode.length <= 0)
		{
			SALSAUtility.printAndThrowRuntimeException("No " + calltype + " Points Selected for Extraction with string " + ToDecode);

		}

		for (String t : ItemsinToDecode)
		{
			String cleanup = t.trim();
			if (cleanup.length() <= 0)
			{
				continue;
			}
			int Startlist, NuminList;
			if (cleanup.contains("-"))
			{
				String[] RangeofItems = cleanup.split("[-]", -1);
				if (RangeofItems.length != 2)
				{
					SALSAUtility.printAndThrowRuntimeException("Illegal " + calltype + " Point Range " + t);

				}
				Startlist = Integer.parseInt(RangeofItems[0]);
				NuminList = Integer.parseInt(RangeofItems[1]) - Integer.parseInt(RangeofItems[0]) + 1;
			}
			else
			{
				Startlist = Integer.parseInt(t);
				NuminList = 1;
			}

			for (int Runoverlist = 0; Runoverlist < NuminList; Runoverlist++)
			{
				VariedorFixedPointList[NumberinVariedorFixedPointList.argValue] = Startlist;
				if ((Startlist < 0) || (Startlist > CollectionMax))
				{
					if (SALSAUtility.MPI_Rank == 0)
					{
						SALSAUtility.SALSAPrint(0, "Illegal " + calltype + " Point Number " + Startlist);
					}
					continue;
				}
				if (CollectionString.equals("originalpoint"))
				{
					SALSAUtility.OriginalPointDisposition[Startlist] = VariedFixedFudge * (SALSAUtility.SALSASHIFT + NumberinVariedorFixedPointList.argValue);
				}
				++NumberinVariedorFixedPointList.argValue;
				++Startlist;
			}
		}
		if (NumberinVariedorFixedPointList.argValue <= 0)
		{
			SALSAUtility.printAndThrowRuntimeException("No " + calltype + " Points Selected for Extraction with string " + ToDecode);

		}
	} // End InterpretIntegerSelection

	//  Set Labels of Original Points
	public static void SetDataLabels(String[] GlobalDataLabels)
	{
		for (int GlobalIndex = 0; GlobalIndex < SALSAUtility.NumberOriginalPoints; GlobalIndex++)
		{
            GlobalDataLabels[GlobalIndex] = String.valueOf(GlobalIndex);
		}

	} // End SetDataLabels

	public static boolean AreValuesSet(SALSADataPointProperties[] DataPropertiesArray)
	{
		boolean valuesareset = true;
		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			int OriginalIndex = SALSAUtility.UsedPointtoOriginalPointMap[PointIndex];
			if (!DataPropertiesArray[OriginalIndex].valuesset)
			{
				valuesareset = false;
				break;
			}
		}
		return valuesareset;
	} // End bool AreValuesSet

	public static boolean AreErrorsSet(SALSADataPointProperties[] DataPropertiesArray)
	{
		boolean errorsareset = true;
		for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
		{
			int OriginalIndex = SALSAUtility.UsedPointtoOriginalPointMap[PointIndex];
			if (!DataPropertiesArray[OriginalIndex].errorsset)
			{
				errorsareset = false;
				break;
			}
		}
		return errorsareset;
	} // End bool AreValuesSet


} // End SALSA_ProcessVariedandFixed
 // End Namespace SALSALibrary
