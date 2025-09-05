# SXI_L3
Code to read Steve's new L3 SXI Product and rebin it in time and spatially. 

Author: S. J. Wharton 
Date: 1st August 2025 

You will need astropy in order to read the FITS files. Try: 
pip install astropy 

The code contains the following scripts in the SXI_L3 folder: 
paths.py                Functions to hardcode certain paths to the data. 
rebin_funcs.py          Functions to rebin data spatially. 
read_cmap.py            Function to get the lundi colourmap. 
make_image_axes.py      Function to plot an extension onto an axes and add a colourbar. 
read_bkg_file.py        Function to read in the background file. 
read_tot_file.py        Function to read in the total counts file. 
read_vcy_file.py        Function to read in the vignetting file. 
rebin_file_spatially.py Function to read in the 3 file types and spatially rebin and calc CXFOV.
rebin_files.py          Function to rebin a set of files both spatially and temporally. 
read_rebinned_file.py   Function to read in the rebinned file. 

TO RUN THE EXAMPLE SCRIPT IN THE UNIX TERMINAL, TYPE: 
python3 examples.py 

REDUNDANT FILES. 
read_mapL3_OLD.py
read_cxfov_OLD.py
read_cxfov_time_OLD.py
examples_OLD.py 


INPUT/OUTPUT. 
The paths for locating the L3 files and putting the output files can be changed in paths.py. They are currently set to: 

datapath = '/data/smile/shared/sims/GAC_OUTPUT/' 
outpath = '/data/smile/shared/sims/GAC_OUTPUT/binned_examples/' 
