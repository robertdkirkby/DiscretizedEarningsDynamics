# DiscretizedEarningsDynamics
Codes for paper introducing discretization method for life-cycle AR(1) with gaussian-mixture innovations.

The actual discretization commands are in VFI Toolkit: https://github.com/vfitoolkit/VFIToolkit-matlab/tree/master/DiscretizationMethods

Following is about generating all the results in the paper.

First run 'JustDiscretizeAndEvaluate.m'
This creates all the discretizations, including a variety of grid sizes (n_z) and other settings.
It also produces lots of output that evaluates the different discretizations.
For each of the six models it is calling 'DiscretizeEarningsDynamicsModel.m' and then 'EvaluateDiscretization.m'.

Then run EarningsDynamics.m six times, changing 'useModel' each time (set useModel=1,2,3,4,5,6). (line 10)
I do this with LoadDiscretization=1 (so it just uses the discretization created earlier)
This calibrates each model, then produces lots of model outputs and saves them in .mat files.
The EarningsDynamics.m can be run in pieces using 'preCalib','doCalib','SolveVandStationaryDist' and 'CalculateStatistics', which can be set to zero to just reload what was done on a previous run.

Then run 'EarningsDynamics_CreateResults.m' which loads things created previously and uses them to create some graphs and tables.
(Many graphs and tables are also created by the earlier codes)

'CompareAccuracyOfGrids.m' creates the Appendix of the paper about accuracy of the discretization for different grid sizes.

Note: For reasons unclear 'JustDiscretizeAndEvaluate.m' runs fine on my laptop but fails on my desktop. Matlab errors are unclear on why this might be?

'DiscretizedEarningDynamics_ExtraAboutz_nz51.pdf' contains a large number of further results on the performance of the discretization method for life-cycle AR(1) with gaussian-mixture shocks (here with 51 grid points).

The codes assume the existence of a bunch of subfolders, the following is a list (hopefully complete) of them:
./SavedOutput/BasicDiscretization/
./SavedOutput/CompareGridSize/
./SavedOutput/EvaluateDiscretization/
./SavedOutput/Main/
./SavedOutput/LatexInputs/
./SavedOutput/Graphs/

