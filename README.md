# SXI_L3
Code to read Steve's new L3 SXI Product and rebin it in time and spatially. 

Author: S. J. Wharton 
Date: 1st August 2025 

You will need astropy in order to read the FITS files. Try: 
pip install astropy 

The code contains the following scripts in the SXI_L3 folder: 
paths           Functions to hardcode certain paths to the data. 
read_mapL3      Class to read in the original L3 SXI file. 
read_cxfov      Class to read to rebin a single L3 SXI file spatially and produce a 
                new fits file. 
read_cxfov_time Class to read in a sequence of consecutive files and combine the data 
                temporally. Then rebins spatially and produces a new fits file. 
rebin_funcs     Functions to rebin data spatially. 
read_cmap       Function to get the lundi colourmap. 
make_image_axes Function to plot an extension onto an axes and add a colourbar. 

TO RUN THE EXAMPLE SCRIPT IN THE UNIX TERMINAL, TYPE: 
python3 examples.py 


INPUT/OUTPUT. 
The paths for locating the L3 files and putting the output files can be changed in paths.py. They are currently set to: 

datapath = '/data/smile/shared/sims/GAC_OUTPUT/' 
outpath = '/data/smile/shared/sims/GAC_OUTPUT/binned_examples/' 
