# Farahani2021

MATLAB Scripts Associated with Farahani et al Cell Reports 2021

1. Steady-state Erk dynamics — scripts for analyzing KTR-reported Erk activity under steady-state conditions in growth medium. These scripts require raw KTR trajectories and background fluorescence intensity to be entered into an Excel sheet (1_KTR_traces.xlsx), in a folder labeled “KTR measurements” (sample data from Figure 1 included here). The peakfinder plugin (https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate) should be saved in a folder named “utils.”

2. EGF dose response — scripts for analyzing KTR-reported Erk activity for cells treated with EGF. Raw KTR trajectories and background fluorescence intensities should be entered into spreadsheets in a folder labeled “KTR measurements” (sample data from Figure 2 included here). 

3. EGF internalization — scripts for analyzing EGF-488 internalization from z-stack images. Raw .tif files for EGF-488 and EGFR should be saved such that files end with “_EGF.tif” and “_EGFR.tif”, respectively. Script saves a binary z-stack images of segmented samples for validating size and intensity thresholding of objects.

4. EGF membrane binding — scripts for analyzing EGF-488 membrane binding from max intensity projections of z-stack images. Raw z-stack images must first be subtracted of background and gaussian filtered using “1_EGF_Membrane_BGsubtract_Blur” macro in ImageJ. “2_EGF_Membrane_Analysis.m” requires 1) .txt files of cell outlines be fed to “EGF_Membrane_Outlines.m”, which converts .txt of outlines (provided by CellPose) into tiffs; 2) EGF-488 and EGFR tiffs to be saved with “EGF” and “EGFR” in their names, respectively. “3_EGF_Membrane_Compile.m” was used to pool together analyzed data from all experiments. 
