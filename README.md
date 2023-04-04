# DiscretizedEarningsDynamics
Codes for paper introducing discretization method for life-cycle AR(1) with gaussian-mixture innovations
The actual discretization commands are in VFI Toolkit: https://github.com/vfitoolkit/VFIToolkit-matlab/tree/master/DiscretizationMethods

Following is about generating all the results in the paper.

First run 'JustDiscretizeAndEvaluate.m'
This creates all the discretizations, including a variety of grid sizes (n_z) and other settings.
It also produces lots of output that evaluates the different discretizations.
For each of the six models it is running 'DiscretizeEarningsDynamicsModel.m' and then 'EvaluateDiscretization.m'.

Then run EarningsDynamics.m six times, changing 'useModel' each time (1,2,3,4,5,6). (line 5)
I do this with LoadDiscretization=1 (so it just uses the discretization created earlier)
I run all six with 'Setmedianwage=0' (line 13)
I then run all six again with 'Setmedianwage=1' (line 13) [This is to calibrate the parameter w]

Then run 'EarningsDynamics2_CreateResults.m' which loads things created previously and uses them to create some graphs and tables.
(Many graphs and tables are also created by the earlier codes)

Note: For reasons unclear 'JustDiscretizeAndEvaluate.m' runs fine on my laptop but fails on my desktop. Matlab errors are unclear on why this might be?
