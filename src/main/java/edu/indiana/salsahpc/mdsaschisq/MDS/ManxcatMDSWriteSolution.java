﻿package MDS;

import Manxcat.*;
import Salsa.Core.Configuration.Sections.*;
import SALSALibrary.*;

public class ManxcatMDSWriteSolution
{

	// Write label-cluster results into a file
	public static void WriteMDS(String fname, int OutputOption, double[][] param, double[][] perr)
	{
		int LocalVectorDimension = param[0].GetLength(0);

		if (LocalVectorDimension != perr[0].GetLength(0))
		{
			SALSAUtility.SALSAError("Inconsistent Dimensions Labels " + (new Integer(LocalVectorDimension)).toString() + " Perr " + perr[0].GetLength(0).toString());
		}

		ManxcatSection Configuration = ManxcatCentral.Configuration;

		// Write SALSA Properties File Header
		if (OutputOption == 0)
		{
			int filetype = 2;

			if ((SALSAUtility.NumberFixedPoints > 0) || (SALSAUtility.NumberVariedPoints < SALSAUtility.NumberOriginalPoints))
			{
				filetype = 4;
			}

			String cleandate = SALSAUtility.startTime.ToLocalTime().toString();
			cleandate = cleandate.replace(":", ".");
			cleandate = cleandate.replace(" ", "-");
			SALSAUtility.GlobalFileProperties.GroupName = "MDSasChisq-" + Configuration.RunSetLabel + "-Run" + Configuration.RunNumber.toString() + "-Date-" + cleandate;
			SALSAUtility.GlobalFileProperties.FileGenerationType = filetype;
			SALSAUtility.GlobalFileProperties.OriginalPointStartIndex = 0;
			SALSAUtility.GlobalFileProperties.LocalPointStartIndex = 0;
			SALSAUtility.GlobalFileProperties.NumberOriginalPoints = SALSAUtility.NumberOriginalPoints;
			SALSAUtility.GlobalFileProperties.NumberPointsinFile = SALSAUtility.PointCount_Global;

			// Comment should have key features of Run
			String OldComment = SALSAUtility.GlobalFileProperties.Comment;

			if (!OldComment.equals(""))
			{
				OldComment += "\n";
			}
			OldComment += "MDSasChisq " + Configuration.RunNumber.toString() + " StartTime " + cleandate + " ****";
			OldComment += "\n Distance Input Option " + Configuration.DistanceProcessingOption.toString() + " Distance Formula " + Configuration.DistanceFormula.toString() + " Weighting Option " + Configuration.Chisqnorm.toString() + " Minimum Distance Value " + Configuration.MinimumDistance.toString();
			OldComment += "\n" + Hotsun.HotsunComment;
			SALSAUtility.GlobalFileProperties.Comment = OldComment;

			for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
			{
				int OriginalPointIndex = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
				int groupnumber = -1;
				String grouplabel = "Ignored";
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].FixedorVaried = 0;
				if (ManxcatMDS.PointStatus[GlobalPointIndex] == -1)
				{
					groupnumber = 3;
					grouplabel = "Deleted";
				}
				else
				{
					if (SALSAUtility.OriginalPointDisposition[GlobalPointIndex] <= -SALSAUtility.SALSASHIFT)
					{
						groupnumber = 2;
						grouplabel = "Fixed";
						SALSAUtility.GlobalPointProperties[GlobalPointIndex].FixedorVaried = 2;
					}

					if (SALSAUtility.OriginalPointDisposition[GlobalPointIndex] >= SALSAUtility.SALSASHIFT)
					{
						groupnumber = 1;
						grouplabel = "Varied";
						SALSAUtility.GlobalPointProperties[GlobalPointIndex].source = "MDSasChisq-" + Configuration.RunSetLabel + "-Run" + Configuration.RunNumber.toString() + "-Date-" + cleandate;
						SALSAUtility.GlobalPointProperties[GlobalPointIndex].FixedorVaried = 1;
					}
				}
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].group = groupnumber;
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].grouplabel = grouplabel;
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].OriginalPointNumber = OriginalPointIndex;
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].LocalPointNumber = GlobalPointIndex;
			}

			for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
			{
				if (ManxcatMDS.PointStatus[GlobalPointIndex] == -1)
				{
					SALSAUtility.GlobalPointProperties[GlobalPointIndex].valuesset = false;
					SALSAUtility.GlobalPointProperties[GlobalPointIndex].errorsset = false;
					continue;
				}
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].x = param[GlobalPointIndex][0];
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].y = param[GlobalPointIndex][1];

				if (LocalVectorDimension > 2)
				{
					SALSAUtility.GlobalPointProperties[GlobalPointIndex].z = param[GlobalPointIndex][2];
				}
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].xerr = perr[GlobalPointIndex][0];
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].yerr = perr[GlobalPointIndex][1];

				if (LocalVectorDimension > 2)
				{
					SALSAUtility.GlobalPointProperties[GlobalPointIndex].zerr = perr[GlobalPointIndex][2];
				}
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].valuesset = true;
				SALSAUtility.GlobalPointProperties[GlobalPointIndex].errorsset = true;
			}
			if (SALSAUtility.Usedreordered)
			{
				SALSADataPointProperties[] NaivePointProperties = new SALSADataPointProperties[SALSAUtility.PointCount_Global];
				for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
				{
					int NaivePointIndex = SALSAUtility.ActualtoNaiveUsedOrder[GlobalPointIndex];
					NaivePointProperties[NaivePointIndex] = SALSAUtility.GlobalPointProperties[GlobalPointIndex].ShallowCopy();
				}
				SALSALibrary.SALSA_Properties.WriteDataPointFile(fname, ManxcatCentral.Configuration.Write2Das3D, "all ", SALSAUtility.GlobalFileProperties, NaivePointProperties, SALSAUtility.PointCount_Global);
				return;
			}
			SALSALibrary.SALSA_Properties.WriteDataPointFile(fname, ManxcatCentral.Configuration.Write2Das3D, "all ", SALSAUtility.GlobalFileProperties, SALSAUtility.GlobalPointProperties, SALSAUtility.PointCount_Global);
			return;
		}

		//  Option OutputOption = 1
		// Simple output of Used in ORIGINAL not Loadbalanced Order
		try
		{
			StreamWriter sw = null;

			if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(fname))
			{
				sw = new StreamWriter(fname, false, Encoding.UTF8);
			}

			if (sw != null)
			{
				for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
				{
					if (ManxcatMDS.PointStatus[GlobalPointIndex] == -1)
					{
						continue;
					}
					String Coordinates = "";
					int SingleCluster = 1;
					int UsedPointIndex = SALSAUtility.NaivetoActualUsedOrder[GlobalPointIndex];
					for (int LocalVectorIndex = 0; LocalVectorIndex < LocalVectorDimension; LocalVectorIndex++)
					{
						Coordinates += String.format("%0.4E", param[UsedPointIndex][LocalVectorIndex]) + "\t";
					}

					sw.WriteLine(String.format((GlobalPointIndex).toString() + "\t" + Coordinates + (new Integer(SingleCluster)).toString()));
				}
			}

			sw.Flush();
			sw.Close();
		}
		catch (RuntimeException e)
		{
			System.out.println("Failed writing data" + e);
			throw (e);
		}
	} // End WriteMDS

} // End class ManxcatMDSWriteSolution // End Namespace MDS