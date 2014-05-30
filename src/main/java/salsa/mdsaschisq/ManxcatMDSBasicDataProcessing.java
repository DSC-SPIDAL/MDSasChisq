package salsa.mdsaschisq;

import Manxcat.*;
import Salsa.Core.*;
import MPI.*;
import SALSALibrary.*;

public class ManxcatMDSBasicDataProcessing
{
	// Initial distance analysis looped over till no more deletions
	public static void IntialDistanceAnalysis(tangible.RefObject<Double> AverageDistance, tangible.RefObject<Double> MaxDistance, tangible.RefObject<Double> SigmaDistance, tangible.RefObject<Integer> TotalDisconnectedPoints, tangible.RefObject<Double> TotalMissingDistances)
	{

		GlobalReductions.FindDoubleSum FindMissingDistances = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);
		GlobalReductions.FindIntSum FindDisconnectedPoints = new GlobalReductions.FindIntSum(SALSAUtility.ThreadCount);
		GlobalReductions.FindDoubleMax FindSystemMax = new GlobalReductions.FindDoubleMax(SALSAUtility.ThreadCount);
		GlobalReductions.FindMeanSigma FindSystemMeanSigma = new GlobalReductions.FindMeanSigma(SALSAUtility.ThreadCount);
		GlobalReductions.FindDoubleSum FindLinkswithoutCut = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount); // Links without Cut
		GlobalReductions.FindDoubleSum FindLinkswithCut = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount); // Links with cut

		int HistogramSize = ManxcatMDS.NumLinksHistogramBins;
		GlobalReductions.FindDoubleArraySum FindLinkHistogramBinCounts = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, HistogramSize);

		int[] DeletedPoint = new int[SALSAUtility.PointCount_Global]; // Set to zero if Point used = 1 if deleted
		for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
		{
			DeletedPoint[GlobalPointIndex] = 0;
		}

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
					// Unset Distances if either point on deleted list
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			FindLinkHistogramBinCounts.startthread(ThreadNo);
			for (int DistributedPointIndex = beginpoint; DistributedPointIndex < indexlen + beginpoint; DistributedPointIndex++)
			{
				int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
				int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
				if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] < SALSAUtility.SALSASHIFT)
				{
					continue;
				}
					int Countlinks = 0;
					for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < SALSAUtility.PointCount_Global; GlobalPointIndex2++)
					{
						int OriginalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex2];
						if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] < SALSAUtility.SALSASHIFT)
						{
							continue;
						}
							if (GlobalPointIndex1 == GlobalPointIndex2)
							{
								continue;
							}
								double tmp2 = SALSAParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
								if (tmp2 < -0.5)
								{
									FindMissingDistances.addapoint(ThreadNo, 1.0);
									continue;
								}
								if ((ManxcatMDS.PointStatus[GlobalPointIndex1] == -1) || (ManxcatMDS.PointStatus[GlobalPointIndex2] == -1))
								{
									FindMissingDistances.addapoint(ThreadNo, 1.0);
									SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2, -1.0);
									continue;
								}
								FindSystemMeanSigma.addapoint(ThreadNo, tmp2);
								FindSystemMax.addapoint(ThreadNo, tmp2);
								++Countlinks;
					}
					if (Countlinks <= SALSAUtility.LinkCut)
					{
						FindDisconnectedPoints.addapoint(ThreadNo, 1);
						DeletedPoint[GlobalPointIndex1] = 1;
					}
					if (Countlinks < HistogramSize)
					{
						FindLinkHistogramBinCounts.addapoint(ThreadNo, Countlinks);
					}
					FindLinkswithoutCut.addapoint(ThreadNo, (double) Countlinks);
					if (Countlinks > SALSAUtility.LinkCut)
					{
						FindLinkswithCut.addapoint(ThreadNo, (double) Countlinks);
					}
			}
		}
	   );

		FindMissingDistances.sumoverthreadsandmpi();
		FindDisconnectedPoints.sumoverthreadsandmpi();
		FindSystemMax.sumoverthreadsandmpi();
		FindSystemMeanSigma.sumoverthreadsandmpi();
		FindLinkswithCut.sumoverthreadsandmpi();
		FindLinkswithoutCut.sumoverthreadsandmpi();
		FindLinkHistogramBinCounts.sumoverthreadsandmpi();

		// Reconcile Deleted points across all processes
		if (SALSAUtility.MPI_Size > 1)
		{
			SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming);
			DeletedPoint = SALSAUtility.MPI_communicator.<Integer>Allreduce(DeletedPoint, Operation<Integer>.Add);
			SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming);
		}
		int countdeleted = 0;
		for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
		{
			if (DeletedPoint[GlobalPointIndex] == 0)
			{
				continue;
			}
			if (DeletedPoint[GlobalPointIndex] != 1)
			{
				SALSAUtility.printAndThrowRuntimeException("System Error Inconsistent Deleted point " + (new Integer(GlobalPointIndex)).toString() + " Status " + (new Integer(DeletedPoint[GlobalPointIndex])).toString());

			}
			ManxcatMDS.PointStatus[GlobalPointIndex] = -1;
			++countdeleted;
		}

		TotalDisconnectedPoints.argValue = FindDisconnectedPoints.TotalInt;
		TotalMissingDistances.argValue = FindMissingDistances.Total;

		double tmp = SALSAUtility.NumberVariedPoints;
		double pointsused = tmp * (tmp - 1) - TotalMissingDistances.argValue;
		double localtotalpoints = FindSystemMeanSigma.TotalNumberofPoints;
		if (Math.abs(localtotalpoints - pointsused) > 2.0)
		{
			SALSAUtility.printAndThrowRuntimeException("System Error and Must stop as Illegal Point Counts " + String.format("%.0f", pointsused) + " " + String.format("%.0f", localtotalpoints));

		}

		AverageDistance.argValue = FindSystemMeanSigma.Totalmean;
		SigmaDistance.argValue = FindSystemMeanSigma.Totalsigma;
		MaxDistance.argValue = FindSystemMax.TotalMax;



		// Output Link Statistics
		double TotalLinkswithoutCut = FindLinkswithoutCut.Total;
		double TotalLinkswithCut = FindLinkswithCut.Total;
		double TotalPointsbeforeLinkCut = FindLinkswithoutCut.TotalNumberofPoints;
		double TotalPointsafterLinkCut = FindLinkswithCut.TotalNumberofPoints;
		double AverageLinksperpoint = TotalLinkswithoutCut / TotalPointsbeforeLinkCut;
		double AverageLinksperpointaftercut = TotalLinkswithCut / TotalPointsafterLinkCut;

		if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
		{
			double fraction = TotalMissingDistances.argValue / (tmp * (tmp - 1.0));
			double EstimatedDimension = 2.0 * AverageDistance.argValue * AverageDistance.argValue / (SigmaDistance.argValue * SigmaDistance.argValue);
			SALSAUtility.SALSAPrint(1, "\nAs INPUT Disconnected Points " + (new Integer(TotalDisconnectedPoints.argValue)).toString() + " Missing Distances " + String.format("%.0f", TotalMissingDistances.argValue) + " Fraction " + String.format("%.4f", fraction));
			SALSAUtility.SALSAPrint(1, "As INPUT Max " + String.format("%.4E", MaxDistance.argValue) + " Average " + String.format("%.4E", AverageDistance.argValue) + " Sigma " + String.format("%.4E", SigmaDistance.argValue) + " Estimated Dimension " + String.format("%.2f", EstimatedDimension));
			SALSAUtility.SALSAPrint(1, " Overall Average Links " + String.format("%.2f", AverageLinksperpoint) + " from " + String.format("%.0f", TotalPointsbeforeLinkCut) + " Points and After cut " + String.format("%.2f", AverageLinksperpointaftercut) + " From Points " + String.format("%.0f", TotalPointsafterLinkCut));

			//  Output Link Histogram
			String histogramcounts = "";
			double total = 0.0;
			for (int binloop = 0; binloop < ManxcatMDS.NumLinksHistogramBins; binloop++)
			{
				histogramcounts += String.format("%.0f", FindLinkHistogramBinCounts.TotalSum[binloop]) + ", ";
				total += FindLinkHistogramBinCounts.TotalSum[binloop];
			}
			SALSAUtility.SALSAPrint(1, "\nLink Histogram Total " + String.format("%.0f", total) + " upto links " + (new Integer(ManxcatMDS.NumLinksHistogramBins)).toString() + " #Counts starting at zero\n" + histogramcounts);

			String deletedpoints = "\nNumber Deleted Points " + (new Integer(countdeleted)).toString() + " Deleted Points ";
			for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
			{
				if (ManxcatMDS.PointStatus[GlobalPointIndex] != -1)
				{
					continue;
				}
				deletedpoints += (new Integer(GlobalPointIndex)).toString() + ", ";
			}
			SALSAUtility.SALSAPrint(1, deletedpoints);
		}
		return;

	} //End IntialDistanceAnalysis

	// Second distance analysis designed to have same interface as IntialDistanceAnalysis
	public static void CleanandTransformDistances(tangible.RefObject<Double> AverageDistance, tangible.RefObject<Double> MaxDistce, tangible.RefObject<Double> SigmaDistance, tangible.RefObject<Integer> TotalDisconnectedPoints, tangible.RefObject<Double> TotalMissingDistances)
	{
		double initialAverageDistance = AverageDistance.argValue;
		double initialMaxDistance = MaxDistce.argValue;
		double initialSigmaDistance = SigmaDistance.argValue;

		double initialEstimatedDimension = 0.0;
		double initialScaleFactor = 0.0, initialTransformedMaxDistance = 0.0;

		if (SALSAUtility.TransformMethod == 8 || SALSAUtility.TransformMethod == 9)
		{
			// 4D or SQRT(4D) transformation specific statistics
			initialEstimatedDimension = 2.0 * initialAverageDistance * initialAverageDistance / (initialSigmaDistance * initialSigmaDistance);
			double initialIndividualSigma = Math.sqrt(initialAverageDistance / initialEstimatedDimension);
			initialScaleFactor = 2.0 * initialIndividualSigma * initialIndividualSigma;

			initialTransformedMaxDistance = initialScaleFactor * Transform4D(SpecialFunction.igamc(initialEstimatedDimension * 0.5, initialMaxDistance / initialScaleFactor));
			if (SALSAUtility.TransformMethod == 9)
			{
				initialTransformedMaxDistance = Math.sqrt(initialTransformedMaxDistance);
			}
		}

		if (SALSAUtility.NumberDistanceWeightCuts > 0)
		{
			for (int weightbin = 0; weightbin <= SALSAUtility.NumberDistanceWeightCuts; weightbin++)
			{
				SALSAUtility.ActualWeightCuts[weightbin] = 0.0;
			}
		}
		GlobalReductions.FindDoubleSum FindMissingDistances = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);
		GlobalReductions.FindIntSum FindDisconnectedPoints = new GlobalReductions.FindIntSum(SALSAUtility.ThreadCount);
		GlobalReductions.FindDoubleMax FindSystemMax = new GlobalReductions.FindDoubleMax(SALSAUtility.ThreadCount);
		GlobalReductions.FindMeanSigma FindSystemMeanSigma = new GlobalReductions.FindMeanSigma(SALSAUtility.ThreadCount);
		GlobalReductions.FindCorrelation FindNewoldDistancestatistics = new GlobalReductions.FindCorrelation(SALSAUtility.ThreadCount);
		GlobalReductions.FindDoubleArraySum FindDistanceWeightBinCounts = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, 1 + SALSAUtility.NumberDistanceWeightCuts);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
					// Unset Distances if either point on deleted list
						// Transform Distances if Requested
							// log(1-d) transformation
							// (x ** power with power given by TransformParameter
							// 4D Transformation. 
							// Note. Works only when no missing distances
							// SQRT(4D) Transformation. 
							// Note. Works only when no missing distances
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			FindDistanceWeightBinCounts.startthread(ThreadNo);
			for (int DistributedPointIndex = beginpoint; DistributedPointIndex < indexlen + beginpoint; DistributedPointIndex++)
			{
				int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
				int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
				if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] < SALSAUtility.SALSASHIFT)
				{
					continue;
				}
					int Countlinks = 0;
					for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < SALSAUtility.PointCount_Global; GlobalPointIndex2++)
					{
						int OriginalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex2];
						if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] < SALSAUtility.SALSASHIFT)
						{
							continue;
						}
							if (GlobalPointIndex1 == GlobalPointIndex2)
							{
								continue;
							}
								double tmp2 = SALSAParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
								if (tmp2 < -0.5)
								{
									FindMissingDistances.addapoint(ThreadNo, 1.0);
									continue;
								}
								if ((ManxcatMDS.PointStatus[GlobalPointIndex1] == -1) || (ManxcatMDS.PointStatus[GlobalPointIndex2] == -1))
								{
									FindMissingDistances.addapoint(ThreadNo, 1.0);
									SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2, -1.0);
								}
								else
								{
									double newvalue = tmp2;
									if (SALSAUtility.TransformMethod == 12)
									{
										tmp2 = Math.min(0.975, tmp2);
										newvalue = Math.log(1 - Math.pow(tmp2, SALSAUtility.TransformParameter)) / Math.log(1 - Math.pow(0.975, SALSAUtility.TransformParameter));
										SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2, newvalue);
									}
									else if (SALSAUtility.TransformMethod == 11)
									{
										tmp2 = Math.min(1.0, tmp2);
										newvalue = 1.0 - Math.pow(1.0 - tmp2, SALSAUtility.TransformParameter);
										SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2, newvalue);
									}
									else if (SALSAUtility.TransformMethod == 10)
									{
										tmp2 = Math.min(1.0, tmp2);
										newvalue = Math.pow(tmp2, SALSAUtility.TransformParameter);
										SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2, newvalue);
									}
									else if (SALSAUtility.TransformMethod == 8)
									{
										tmp2 = initialScaleFactor * Transform4D(SpecialFunction.igamc(initialEstimatedDimension * 0.5, tmp2 / initialScaleFactor));
										newvalue = tmp2 / initialTransformedMaxDistance;
										SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2, newvalue);
									}
									else if (SALSAUtility.TransformMethod == 9)
									{
										tmp2 = initialScaleFactor * Transform4D(SpecialFunction.igamc(initialEstimatedDimension * 0.5, tmp2 / initialScaleFactor));
										newvalue = Math.sqrt(tmp2) / initialTransformedMaxDistance;
										SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2, newvalue);
									}
									FindNewoldDistancestatistics.addapoint(ThreadNo, tmp2, newvalue);
									FindSystemMax.addapoint(ThreadNo, newvalue);
									FindSystemMeanSigma.addapoint(ThreadNo, newvalue);
									++Countlinks;
									if (SALSAUtility.NumberDistanceWeightCuts > 0)
									{
										int distpos = SALSAUtility.NumberDistanceWeightCuts;
										for (int weightbin = 0; weightbin < SALSAUtility.NumberDistanceWeightCuts; weightbin++)
										{
											if (SALSAUtility.DistanceWeightCuts[weightbin] > newvalue)
											{
												distpos = weightbin;
												break;
											}
										}
										FindDistanceWeightBinCounts.addapoint(ThreadNo, distpos);
									}
								}
					}
					if ((Countlinks <= 0) && (ManxcatMDS.PointStatus[GlobalPointIndex1] != -1))
					{
						FindDisconnectedPoints.addapoint(ThreadNo, 1);
					}
			}
		}
	   );

		// Accumulate Threads and invoke MPI
		FindNewoldDistancestatistics.sumoverthreadsandmpi();
		FindSystemMeanSigma.sumoverthreadsandmpi();
		FindSystemMax.sumoverthreadsandmpi();
		FindDisconnectedPoints.sumoverthreadsandmpi();
		FindMissingDistances.sumoverthreadsandmpi();
		FindDistanceWeightBinCounts.sumoverthreadsandmpi();

		AverageDistance.argValue = FindSystemMeanSigma.Totalmean;
		SigmaDistance.argValue = FindSystemMeanSigma.Totalsigma;
		MaxDistce.argValue = FindSystemMax.TotalMax;
		TotalDisconnectedPoints.argValue = FindDisconnectedPoints.TotalInt;
		TotalMissingDistances.argValue = FindMissingDistances.Total;
		double tmp = SALSAUtility.NumberVariedPoints;
		double pointsused = tmp * (tmp - 1) - TotalMissingDistances.argValue;
		double localtotalpoints = FindSystemMeanSigma.TotalNumberofPoints;
		if (Math.abs(localtotalpoints - pointsused) > 2.0)
		{
			SALSAUtility.printAndThrowRuntimeException("System Error and Must stop as Illegal Point Counts " + String.format("%.0f", pointsused) + " " + String.format("%.0f", localtotalpoints));

		}

		// Print correlations
		String label = "";
		if (SALSAUtility.TransformMethod != 0)
		{
			label = "\nTransform Method " + (new Integer(SALSAUtility.TransformMethod)).toString() + " Parameter " + String.format("%.4f", SALSAUtility.TransformParameter) + " ";
		}

		FindNewoldDistancestatistics.print(label + "\nOld and Transformed Distances\n", "%.4f");
		if (TotalDisconnectedPoints.argValue > 0)
		{
			SALSAUtility.printAndThrowRuntimeException("System Error and Must stop as Disconnected Points with ZERO links after first round " + (new Integer(TotalDisconnectedPoints.argValue)).toString());

		}

		// Set Distance Weightings
		if (SALSAUtility.NumberDistanceWeightCuts <= 0)
		{
			return;
		}
		String distanceresults = "";
		double fudge = 1.0 + (double) SALSAUtility.NumberDistanceWeightCuts;
		for (int weightbin = 0; weightbin <= SALSAUtility.NumberDistanceWeightCuts; weightbin++)
		{
			double weightforbin = 0.0;
			if (FindDistanceWeightBinCounts.TotalSum[weightbin] > 0.5)
			{
				weightforbin = pointsused / (fudge * FindDistanceWeightBinCounts.TotalSum[weightbin]);
			}
			String labelreason = " Rest";
			if (weightbin < SALSAUtility.NumberDistanceWeightCuts)
			{
				labelreason = " Distce < " + String.format("%.3f", SALSAUtility.DistanceWeightCuts[weightbin]);
			}
			distanceresults += String.format("%.0f", FindDistanceWeightBinCounts.TotalSum[weightbin]) + " Pts weight " + String.format("%.3f", weightforbin) + labelreason + " * ";
			SALSAUtility.ActualWeightCuts[weightbin] = Math.sqrt(weightforbin);
		}
		SALSAUtility.SALSAPrint(1, "\nDistance Counts and Cuts " + distanceresults);
		return;

	}

	private static double Transform4D(double higherDimensionInput)
	{
		/* Saliya - copying code directly from Adam Hughes project. 
		 * The comments except this one are from his code as well.*/

		double Einit = -Math.log(higherDimensionInput);
		//         double Einit = HigherDimInput;
		double E = Einit;
		double Eold = E;
		// double diff;

		for (int recurse = 0; recurse < 50; recurse++)
		{
			E = Einit + Math.log(1.0 + E);
			//      E = 1.0 + E - E / 2;
			/*      diff = E - Eold;
			      if (diff < 0)
			      {
			          diff = Eold - E;
			      }*/
			//                if (diff < 0.00001)
			if (Math.abs(E - Eold) < 0.00001)
			{
				return E;
			}
			Eold = E;
		}
		return E;
	}

//End CleanandTransformDistances

	// Find Center defined as Point with minimum mean distance to other points
	// Form Histogram of distance Bin counts
	public static void FindCenter(tangible.RefObject<Integer> Center, tangible.RefObject<Double> Radius, double Histmin, double Histmax, int NumberBins, tangible.RefObject<double[]> Bincounts)
	{
		GlobalReductions.FindMinorMaxValuewithIndex FindCenterCompute = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);

		int HistogramSize = 2 + NumberBins;
		double HistFudge = NumberBins / (Histmax - Histmin);
		GlobalReductions.FindDoubleArraySum FindDistanceHistogramBinCounts = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, HistogramSize);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
					// Form histogram
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			FindDistanceHistogramBinCounts.startthread(ThreadNo);
			for (int DistributedPointIndex = beginpoint; DistributedPointIndex < indexlen + beginpoint; DistributedPointIndex++)
			{
				double distcemeanperpoint = 0.0;
				int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
				if (ManxcatMDS.PointStatus[GlobalPointIndex1] == -1)
				{
					continue;
				}
					int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
					if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] < SALSAUtility.SALSASHIFT)
					{
						continue;
					}
						int Countlinks = 0;
						for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < SALSAUtility.PointCount_Global; GlobalPointIndex2++)
						{
							if (ManxcatMDS.PointStatus[GlobalPointIndex2] == -1)
							{
								continue;
							}
								int OriginalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex2];
								if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] < SALSAUtility.SALSASHIFT)
								{
									continue;
								}
									if (GlobalPointIndex1 == GlobalPointIndex2)
									{
										continue;
									}
										double tmp2 = SALSAParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
										if (tmp2 < -0.5)
										{
											continue;
										}
											distcemeanperpoint += tmp2;
											++Countlinks;
											double Histvalue = (tmp2 - Histmin) * HistFudge;
											int HistPosition = (int)Math.floor(Histvalue);
											if (HistPosition == NumberBins)
											{
												HistPosition = (int)Math.floor(Histvalue -0.00001);
											}
											if (HistPosition < 0)
											{
												HistPosition = -1;
											}
												if (HistPosition >= NumberBins)
												{
													HistPosition = NumberBins;
												}
													FindDistanceHistogramBinCounts.addapoint(ThreadNo, 1 + HistPosition);
						}
						if ((Countlinks >= ManxcatMDS.LinkCutforCenter) && (ManxcatMDS.PointStatus[GlobalPointIndex1] == 0))
						{
							distcemeanperpoint = distcemeanperpoint / Countlinks;
							FindCenterCompute.addapoint(ThreadNo, GlobalPointIndex1, distcemeanperpoint);
						}
			}
		}
	   );

		FindCenterCompute.sumoverthreadsandmpi();
		FindDistanceHistogramBinCounts.sumoverthreadsandmpi();
		Center.argValue = FindCenterCompute.TotalIndexValue;
		Radius.argValue = FindCenterCompute.TotalMaxOrMin;
		for (int binloop = 0; binloop < HistogramSize; binloop++)
		{
			Bincounts.argValue[binloop] = FindDistanceHistogramBinCounts.TotalSum[binloop];
		}
		return;

	} //End FindCenter

	/** 
	 Finds the minimum and maximum distances of original and Eculidean data
	 for both whole and selected set of clusters.
	 
	 @param xminWhole
	 @param xmaxWhole
	 @param yminWhole
	 @param ymaxWhole
	 @param xminSelected
	 @param xmaxSelected
	 @param yminSelected
	 @param ymaxSelected
	 @param ymaxSelectedInter
	 @param originalPnumToCnumTable
	 @param xminSelectedInter
	 @param xmaxSelectedInter
	 @param yminSelectedInter
	*/
	public static void FindMinMaxForDensity(tangible.RefObject<Double> xminWhole, tangible.RefObject<Double> xmaxWhole, tangible.RefObject<Double> yminWhole, tangible.RefObject<Double> ymaxWhole, tangible.RefObject<Double> xminSelected, tangible.RefObject<Double> xmaxSelected, tangible.RefObject<Double> yminSelected, tangible.RefObject<Double> ymaxSelected, tangible.RefObject<Double> xminSelectedInter, tangible.RefObject<Double> xmaxSelectedInter, tangible.RefObject<Double> yminSelectedInter, tangible.RefObject<Double> ymaxSelectedInter, java.util.Hashtable originalPnumToCnumTable)
	{
		/* MinorMax objects with index though index is not used here */
		GlobalReductions.FindMinorMaxValuewithIndex FindXminWhole = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
		GlobalReductions.FindMinorMaxValuewithIndex FindXminSelected = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
		GlobalReductions.FindMinorMaxValuewithIndex FindXminSelectedInter = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
		GlobalReductions.FindMinorMaxValuewithIndex FindYminWhole = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
		GlobalReductions.FindMinorMaxValuewithIndex FindYminSelected = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
		GlobalReductions.FindMinorMaxValuewithIndex FindYminSelectedInter = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);

		GlobalReductions.FindMinorMaxValuewithIndex FindXmaxWhole = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
		GlobalReductions.FindMinorMaxValuewithIndex FindXmaxSelected = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
		GlobalReductions.FindMinorMaxValuewithIndex FindXmaxSelectedInter = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
		GlobalReductions.FindMinorMaxValuewithIndex FindYmaxWhole = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
		GlobalReductions.FindMinorMaxValuewithIndex FindYmaxSelected = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
		GlobalReductions.FindMinorMaxValuewithIndex FindYmaxSelectedInter = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
					/* At this point the xval should be the transformed distances (if specified using TransformMethod & TransformParameter) */
					// Todo (html+density) - see if any distance cut needs to be considered. Also check if any pair is to be removed if not exist in pnumToCnum table
							// Intra cluster pairs (p1,p2) where both p1,p2 belong to one cluster.
							// So check if p1,p2 belong to the same cluster in our set of selected clusters
							// Inter cluster pairs (p1,p2) where both p1,p2 does NOT belong to one cluster.
							// So check if p1,p2 does NOT belong to the same cluster in our set of selected clusters
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int distributedPointIndex = beginpoint; distributedPointIndex < indexlen + beginpoint; distributedPointIndex++)
			{
				int globalPointIndex1 = distributedPointIndex + SALSAUtility.PointStart_Process;
				if (ManxcatMDS.PointStatus[globalPointIndex1] == -1)
				{
					continue;
				}
					int originalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex1];
					int cnum1 = SALSAUtility.IsClustersSelected ? ((int) originalPnumToCnumTable.get(originalPointIndex1)) : -1;
					if (SALSAUtility.OriginalPointDisposition[originalPointIndex1] < SALSAUtility.SALSASHIFT)
					{
						continue;
					}
						int usedPointIndex1 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex1];
						for (int globalPointIndex2 = 0; globalPointIndex2 < SALSAUtility.PointCount_Global; globalPointIndex2++)
						{
							if (ManxcatMDS.PointStatus[globalPointIndex2] == -1)
							{
								continue;
							}
								int originalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex2];
								if (SALSAUtility.OriginalPointDisposition[originalPointIndex2] < SALSAUtility.SALSASHIFT)
								{
									continue;
								}
									if (globalPointIndex1 == globalPointIndex2)
									{
										continue;
									}
										double xval = SALSAParallelism.getDistanceValue(globalPointIndex1, globalPointIndex2);
										if (xval < -0.5)
										{
											continue;
										}
											int usedPointIndex2 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex2];
											double yval = GetEuclideanDistance(Hotsun.GlobalParameter[usedPointIndex1], Hotsun.GlobalParameter[usedPointIndex2]);
											UpdateMinMax(xval, yval, FindXminWhole, FindXmaxWhole, FindYminWhole, FindYmaxWhole, ThreadNo);
											int cnum2 = SALSAUtility.IsClustersSelected ? ((int) originalPnumToCnumTable.get(originalPointIndex2)) : -1;
											if (cnum1 != -1 && cnum2 != -1 && SALSAUtility.SelectedClusters.contains(cnum1) && SALSAUtility.SelectedClusters.contains(cnum2))
											{
												if (cnum1 == cnum2)
												{
													UpdateMinMax(xval, yval, FindXminSelected, FindXmaxSelected, FindYminSelected, FindYmaxSelected, ThreadNo);
												}
												else
												{
													UpdateMinMax(xval, yval, FindXminSelectedInter, FindXmaxSelectedInter, FindYminSelectedInter, FindYmaxSelectedInter, ThreadNo);
												}
											}
						}
			}
		}
	   );

		FindXminWhole.sumoverthreadsandmpi();
		FindXminSelected.sumoverthreadsandmpi();
		FindXminSelectedInter.sumoverthreadsandmpi();
		FindYminWhole.sumoverthreadsandmpi();
		FindYminSelected.sumoverthreadsandmpi();
		FindYminSelectedInter.sumoverthreadsandmpi();

		FindXmaxWhole.sumoverthreadsandmpi();
		FindXmaxSelected.sumoverthreadsandmpi();
		FindXmaxSelectedInter.sumoverthreadsandmpi();
		FindYmaxWhole.sumoverthreadsandmpi();
		FindYmaxSelected.sumoverthreadsandmpi();
		FindYmaxSelectedInter.sumoverthreadsandmpi();

		xminWhole.argValue = FindXminWhole.TotalMaxOrMin;
		xminSelected.argValue = FindXminSelected.TotalMaxOrMin;
		xminSelectedInter.argValue = FindXminSelectedInter.TotalMaxOrMin;
		yminWhole.argValue = FindYminWhole.TotalMaxOrMin;
		yminSelected.argValue = FindYminSelected.TotalMaxOrMin;
		yminSelectedInter.argValue = FindYminSelectedInter.TotalMaxOrMin;

		xmaxWhole.argValue = FindXmaxWhole.TotalMaxOrMin;
		xmaxSelected.argValue = FindXmaxSelected.TotalMaxOrMin;
		xmaxSelectedInter.argValue = FindXmaxSelectedInter.TotalMaxOrMin;
		ymaxWhole.argValue = FindYmaxWhole.TotalMaxOrMin;
		ymaxSelected.argValue = FindYmaxSelected.TotalMaxOrMin;
		ymaxSelectedInter.argValue = FindYmaxSelectedInter.TotalMaxOrMin;
	}

	/** 
	 Given x and y values along their respective current minimum and maximum values, this will update
	 the minimum and maximum values.
	 
	 @param x
	 @param y
	 @param xmin
	 @param xmax
	 @param ymin
	 @param ymax
	 @param threadNo
	*/
	private static void UpdateMinMax(double x, double y, GlobalReductions.FindMinorMaxValuewithIndex xmin, GlobalReductions.FindMinorMaxValuewithIndex xmax, GlobalReductions.FindMinorMaxValuewithIndex ymin, GlobalReductions.FindMinorMaxValuewithIndex ymax, int threadNo)
	{
		// Ignore the index value
		xmin.addapoint(threadNo, 0, x);
		xmax.addapoint(threadNo, 0, x);
		ymin.addapoint(threadNo, 0, y);
		ymax.addapoint(threadNo, 0, y);
	}

	/** 
	 Generates density matrices of Euclidean Vs original distances for the whole
	 sample and the points in selected set of clusters. In selected clusters we
	 consider only the intra-cluster (within the cluster) pairs of points for distances.
	 
	 @param densityMatrixWhole 2D double array of size SALSAUtlity.Xres x SALSAUtility.Yres
	 @param yHistogramSelected
	 @param densityMatrixSelected 2D double array of size SALSAUtlity.Xres x SALSAUtility.Yres
	 @param yHistogramSelectedInter
	 @param xminWhole minimum original distance of a pair in the whole sample
	 @param xmaxWhole maximum original distance of a pair in the whole sample
	 @param yminWhole minimum Euclidean distance of a pair in the whole sample
	 @param ymaxWhole maximum Euclidean distance of a pair in the whole sample
	 @param xminSelected minimum original distance of a intra-cluster pair in the selected clusters
	 @param xmaxSelected maximum original distance of a intra-cluster pair in the selected clusters
	 @param yminSelected minimum Euclidean distance of a intra-cluster pair in the selected clusters
	 @param ymaxSelected maximum Euclidean distance of a intra-cluster pair in the selected clusters
	 @param countWhole
	 @param countSelected
	 @param originalPnumToCnumTable mapping from original point number to cluster number
	 @param xHistogramWhole
	 @param yHistogramWhole
	 @param xHistogramSelected
	 @param densityMatrixSelectedInter
	 @param xHistogramSelectedInter
	*/
	public static void GenerateDensityMatrix(tangible.RefObject<double[][]> densityMatrixWhole, tangible.RefObject<double[]> xHistogramWhole, tangible.RefObject<double[]> yHistogramWhole, tangible.RefObject<double[][]> densityMatrixSelected, tangible.RefObject<double[]> xHistogramSelected, tangible.RefObject<double[]> yHistogramSelected, tangible.RefObject<double[][]> densityMatrixSelectedInter, tangible.RefObject<double[]> xHistogramSelectedInter, tangible.RefObject<double[]> yHistogramSelectedInter, double xminWhole, double xmaxWhole, double yminWhole, double ymaxWhole, double xminSelected, double xmaxSelected, double yminSelected, double ymaxSelected, double xminSelectedInter, double xmaxSelectedInter, double yminSelectedInter, double ymaxSelectedInter, tangible.RefObject<Double> countWhole, tangible.RefObject<Double> countSelected, tangible.RefObject<Double> countSelectedInter, java.util.Hashtable originalPnumToCnumTable)
	{
		GlobalReductions.Find2DDoubleArraySum FindDensityMatrixWhole = new GlobalReductions.Find2DDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Yres, SALSAUtility.Xres);
		GlobalReductions.Find2DDoubleArraySum FindDensityMatrixSelected = new GlobalReductions.Find2DDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Yres, SALSAUtility.Xres);
		GlobalReductions.Find2DDoubleArraySum FindDensityMatrixSelectedInter = new GlobalReductions.Find2DDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Yres, SALSAUtility.Xres);

		GlobalReductions.FindDoubleArraySum FindXHistogramWhole = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Xres);
		GlobalReductions.FindDoubleArraySum FindXHistogramSelected = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Xres);
		GlobalReductions.FindDoubleArraySum FindXHistogramSelectedInter = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Xres);
		GlobalReductions.FindDoubleArraySum FindYHistogramWhole = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Yres);
		GlobalReductions.FindDoubleArraySum FindYHistogramSelected = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Yres);
		GlobalReductions.FindDoubleArraySum FindYHistogramSelectedInter = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount, SALSAUtility.Yres);

		double deltaxWhole = (xmaxWhole - xminWhole) / SALSAUtility.Xres;
		double deltayWhole = (ymaxWhole - yminWhole) / SALSAUtility.Yres;

		double deltaxSelected = (xmaxSelected - xminSelected) / SALSAUtility.Xres;
		double deltaySelected = (ymaxSelected - yminSelected) / SALSAUtility.Yres;

		double deltaxSelectedInter = (xmaxSelectedInter - xminSelectedInter) / SALSAUtility.Xres;
		double deltaySelectedInter = (ymaxSelectedInter - yminSelectedInter) / SALSAUtility.Yres;

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
					/* At this point the xval should be the transformed distances (if specified using TransformMethod & TransformParameter) */
					// Todo (html+density) - see if any distance cut needs to be considered. Also check if any pair is to be removed if not exist in pnumToCnum table
							// Intra cluster pairs (p1,p2) where both p1,p2 belong to one cluster.
							// So check if p1,p2 belong to the same cluster in our set of selected clusters
							// Inter cluster pairs (p1,p2) where both p1,p2 does NOT belong to one cluster.
							// So check if p1,p2 does NOT belong to the same cluster in our set of selected clusters
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			FindDensityMatrixWhole.startthread(ThreadNo);
			FindDensityMatrixSelected.startthread(ThreadNo);
			FindDensityMatrixSelectedInter.startthread(ThreadNo);
			FindXHistogramWhole.startthread(ThreadNo);
			FindXHistogramSelected.startthread(ThreadNo);
			FindXHistogramSelectedInter.startthread(ThreadNo);
			FindYHistogramWhole.startthread(ThreadNo);
			FindYHistogramSelected.startthread(ThreadNo);
			FindYHistogramSelectedInter.startthread(ThreadNo);
			for (int distributedPointIndex = beginpoint; distributedPointIndex < indexlen + beginpoint; distributedPointIndex++)
			{
				int globalPointIndex1 = distributedPointIndex + SALSAUtility.PointStart_Process;
				if (ManxcatMDS.PointStatus[globalPointIndex1] == -1)
				{
					continue;
				}
					int originalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex1];
					int cnum1 = SALSAUtility.IsClustersSelected ? ((int) originalPnumToCnumTable.get(originalPointIndex1)) : -1;
					if (SALSAUtility.OriginalPointDisposition[originalPointIndex1] < SALSAUtility.SALSASHIFT)
					{
						continue;
					}
						int usedPointIndex1 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex1];
						for (int globalPointIndex2 = 0; globalPointIndex2 < SALSAUtility.PointCount_Global; globalPointIndex2++)
						{
							if (ManxcatMDS.PointStatus[globalPointIndex2] == -1)
							{
								continue;
							}
								int originalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex2];
								if (SALSAUtility.OriginalPointDisposition[originalPointIndex2] < SALSAUtility.SALSASHIFT)
								{
									continue;
								}
									if (globalPointIndex1 == globalPointIndex2)
									{
										continue;
									}
										double xval = SALSAParallelism.getDistanceValue(globalPointIndex1, globalPointIndex2);
										if (xval < -0.5)
										{
											continue;
										}
											int usedPointIndex2 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex2];
											double yval = GetEuclideanDistance(Hotsun.GlobalParameter[usedPointIndex1], Hotsun.GlobalParameter[usedPointIndex2]);
											UpdateCells(xval, yval, xmaxWhole, xminWhole, ymaxWhole, yminWhole, deltaxWhole, deltayWhole, FindDensityMatrixWhole, FindXHistogramWhole, FindYHistogramWhole, ThreadNo);
											int cnum2 = SALSAUtility.IsClustersSelected ? ((int) originalPnumToCnumTable.get(originalPointIndex2)) : -1;
											if (cnum1 != -1 && cnum2 != -1 && SALSAUtility.SelectedClusters.contains(cnum1) && SALSAUtility.SelectedClusters.contains(cnum2))
											{
												if (cnum1 == cnum2)
												{
													UpdateCells(xval, yval, xmaxSelected, xminSelected, ymaxSelected, yminSelected, deltaxSelected, deltaySelected, FindDensityMatrixSelected, FindXHistogramSelected, FindYHistogramSelected, ThreadNo);
												}
												else
												{
													UpdateCells(xval, yval, xmaxSelectedInter, xminSelectedInter, ymaxSelectedInter, yminSelectedInter, deltaxSelectedInter, deltaySelectedInter, FindDensityMatrixSelectedInter, FindXHistogramSelectedInter, FindYHistogramSelectedInter, ThreadNo);
												}
											}
						}
			}
		}
	   );

		FindDensityMatrixWhole.sumoverthreadsandmpi();
		FindDensityMatrixSelected.sumoverthreadsandmpi();
		FindDensityMatrixSelectedInter.sumoverthreadsandmpi();

		FindXHistogramWhole.sumoverthreadsandmpi();
		FindXHistogramSelected.sumoverthreadsandmpi();
		FindXHistogramSelectedInter.sumoverthreadsandmpi();
		FindYHistogramWhole.sumoverthreadsandmpi();
		FindYHistogramSelected.sumoverthreadsandmpi();
		FindYHistogramSelectedInter.sumoverthreadsandmpi();

		densityMatrixWhole.argValue = FindDensityMatrixWhole.TotalSum;
		densityMatrixSelected.argValue = FindDensityMatrixSelected.TotalSum;
		densityMatrixSelectedInter.argValue = FindDensityMatrixSelectedInter.TotalSum;

		xHistogramWhole.argValue = FindXHistogramWhole.TotalSum;
		xHistogramSelected.argValue = FindXHistogramSelected.TotalSum;
		xHistogramSelectedInter.argValue = FindXHistogramSelectedInter.TotalSum;
		yHistogramWhole.argValue = FindYHistogramWhole.TotalSum;
		yHistogramSelected.argValue = FindYHistogramSelected.TotalSum;
		yHistogramSelectedInter.argValue = FindYHistogramSelectedInter.TotalSum;

		countWhole.argValue = FindDensityMatrixWhole.TotalNumberofPoints;
		countSelected.argValue = FindDensityMatrixSelected.TotalNumberofPoints;
		countSelectedInter.argValue = FindDensityMatrixSelectedInter.TotalNumberofPoints;
	}

	private static void UpdateCells(double x, double y, double xmax, double xmin, double ymax, double ymin, double deltax, double deltay, GlobalReductions.Find2DDoubleArraySum FindDensityMatrix, GlobalReductions.FindDoubleArraySum FindXHistogram, GlobalReductions.FindDoubleArraySum FindYHistogram, int threadNo)
	{
		// cell number based on zero index from bottom left corner
		// if x is equal to xmax then it's placed in the last cell, which is xres-1 in zero based index
		// same is done for y when y == ymax
		int cellx = Math.abs(x - xmax) < Double.MIN_VALUE ? SALSAUtility.Xres - 1 : (int)Math.floor((x - xmin) / deltax);
		int celly = Math.abs(y - ymax) < Double.MIN_VALUE ? SALSAUtility.Yres - 1 : (int)Math.floor((y - ymin) / deltay);

		if (x > xmax || y > ymax || x < xmin || y < ymin)
		{
			// now this should never be reached
			throw new RuntimeException("bad(1)-> x: " + x + " y: " + y + " xmax: " + xmax + " xmin: " + xmin + " ymax: " + ymax + " ymin: " + ymin);
		}

		if (cellx >= SALSAUtility.Xres || celly >= SALSAUtility.Yres)
		{
			// now this should never be reached
			throw new RuntimeException("bad(2)-> x: " + x + " y:" + y + " xmax: " + xmax + " xmin: " + xmin + " ymax: " + ymax + " ymin: " + ymin + " cellx: " + cellx + " celly: " + celly);
		}

		FindDensityMatrix.addapoint(threadNo, celly, cellx);
		FindXHistogram.addapoint(threadNo, cellx);
		FindYHistogram.addapoint(threadNo, celly);
	}




	public static double GetEuclideanDistance(double[] vector1, double[] vector2)
	{
		if (vector1.length != vector2.length)
		{
			throw new RuntimeException(ManxcatErrorMessages.UnequalVectorLengths);
		}

		int dimension = vector1.length;

		double d = 0.0;
		for (int i = 0; i < dimension; ++i)
		{
			d += Math.pow((vector2[i] - vector1[i]), 2);
		}
		return Math.sqrt(d);
	}


	// Find x Axis position and number of nearby distances and number of points near another
	public static void FindxAxis(int Center, tangible.RefObject<Integer> xAxis, tangible.RefObject<Double> MaxDistceGlobal, tangible.RefObject<Double> DistancesNearEachOther, tangible.RefObject<Integer> NotLonelyPoints)
	{

		GlobalReductions.FindDoubleSum FindNearbyPairs = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);
		GlobalReductions.FindIntSum FindCozyPoints = new GlobalReductions.FindIntSum(SALSAUtility.ThreadCount);
		GlobalReductions.FindMinorMaxValuewithIndex FindxAxisCompute = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int DistributedPointIndex = beginpoint; DistributedPointIndex < indexlen + beginpoint; DistributedPointIndex++)
			{
				int notlonely = 0;
				int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
				if (ManxcatMDS.PointStatus[GlobalPointIndex1] != 0)
				{
					continue;
				}
					int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
					if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] < SALSAUtility.SALSASHIFT)
					{
						continue;
					}
						for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < SALSAUtility.PointCount_Global; GlobalPointIndex2++)
						{
							if (GlobalPointIndex2 == GlobalPointIndex1)
							{
								continue;
							}
								int OriginalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex2];
								if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] < SALSAUtility.SALSASHIFT)
								{
									continue;
								}
									double tmp = SALSAParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
									if (tmp < -0.5)
									{
										continue;
									}
										if (tmp < ManxcatMDS.MinimumDistance)
										{
											++notlonely;
											FindNearbyPairs.addapoint(ThreadNo, 1.0);
										}
						}
						if (notlonely > 0)
						{
							FindCozyPoints.addapoint(ThreadNo, 1);
						}
							double MaxDistcePoint = SALSAParallelism.getDistanceValue(GlobalPointIndex1, Center);
							if (MaxDistcePoint < -0.5)
							{
								continue;
							}
								FindxAxisCompute.addapoint(ThreadNo, GlobalPointIndex1, MaxDistcePoint);
			}
		}
	   );

		FindxAxisCompute.sumoverthreadsandmpi();
		FindCozyPoints.sumoverthreadsandmpi();
		FindNearbyPairs.sumoverthreadsandmpi();

		xAxis.argValue = FindxAxisCompute.TotalIndexValue;
		MaxDistceGlobal.argValue = FindxAxisCompute.TotalMaxOrMin;
		NotLonelyPoints.argValue = FindCozyPoints.TotalInt;
		DistancesNearEachOther.argValue = FindNearbyPairs.Total;
		return;

	} //End FindxAxis

	//  Find point defining xy Plane
	public static void FindxyPlane(int Center, int xAxis, tangible.RefObject<Integer> xyPlane, tangible.RefObject<Double> MaxDistceGlobal)
	{
		GlobalReductions.FindMinorMaxValuewithIndex FindxyPlaneCompute = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);

		Parallel.For(0, SALSAUtility.getParallelOptions().MaxDegreeOfParallelism, SALSAUtility.getParallelOptions(), (ThreadNo) => // End loop over Point dependent quantities
		{
			int indexlen = SALSAUtility.PointsperThread[ThreadNo];
			int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] - SALSAUtility.PointStart_Process;
			for (int DistributedPointIndex = beginpoint; DistributedPointIndex < indexlen + beginpoint; DistributedPointIndex++)
			{
				int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
				int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
				if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] < SALSAUtility.SALSASHIFT)
				{
					continue;
				}
					if (ManxcatMDS.PointStatus[GlobalPointIndex1] != 0)
					{
						continue;
					}
						double tmp1 = SALSAParallelism.getDistanceValue(GlobalPointIndex1, Center);
						if (tmp1 < -0.5)
						{
							continue;
						}
							double tmp2 = SALSAParallelism.getDistanceValue(GlobalPointIndex1, xAxis);
							if (tmp2 < -0.5)
							{
								continue;
							}
								double MaxDistcePoint = tmp1 * tmp2;
								FindxyPlaneCompute.addapoint(ThreadNo, GlobalPointIndex1, MaxDistcePoint);
			}
		}
	   );

		FindxyPlaneCompute.sumoverthreadsandmpi();
		xyPlane.argValue = FindxyPlaneCompute.TotalIndexValue;
		MaxDistceGlobal.argValue = FindxyPlaneCompute.TotalMaxOrMin;
		return;

	} //End FindxyPlane

	public static void SetUpHistogramRange(int PointsinHistogram, tangible.RefObject<Double> Histmin, tangible.RefObject<Double> Histmax)
	{ // Choose good range to get rounded labels

		if (Histmax.argValue <= Histmin.argValue)
		{
			return;
		}
		if (Math.abs(Histmin.argValue) < 0.000000001)
		{
			Histmin.argValue = 0.0;
		}
		else
		{
			double newminvalue = NewBinSize(Math.abs(Histmin.argValue));
			if (Histmin.argValue < 0.0)
			{
				Histmin.argValue = -newminvalue;
			}
			else
			{
				Histmin.argValue = newminvalue;
			}
		}
		double binsize = (Histmax.argValue - Histmin.argValue) / PointsinHistogram;
		Histmax.argValue = Histmin.argValue + PointsinHistogram * NewBinSize(binsize);
		return;

	} // End SetUpHistogramRange

	public static double NewBinSize(double oldBinSize)
	{ // Round Bin size up to a pretty value
		if (oldBinSize <= 0.000000001 || oldBinSize > 10000000000.0)
		{
			return oldBinSize;
		}

		double logvalue = Math.log10(oldBinSize);
		int intlogvalue = (int)Math.floor(logvalue);
		double fudgepower = 1.0 - (double)intlogvalue;
		double fudge = Math.pow(10.0, fudgepower);
		double scaled = fudge * oldBinSize;
		scaled = Math.min(scaled, 100.0);
		double Intversionofnewbinsize = Math.ceil(scaled) / fudge;
//            SALSAUtility.SALSAPrint(1, "Hist " + oldBinSize.ToString("F4") + " " + intlogvalue + " " + Intversionofnewbinsize.ToString("F4") + " scaled " + scaled.ToString("F2"));
		return Intversionofnewbinsize;
	}





}