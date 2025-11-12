#This is an example script for running it. 

import SXI_L3

#CREATE A REBINNED FILE. 
plot = False 


#Enter the start and end times as strings. 
start = '20260317T0240'
end = '20260317T0245'

#Enter the spatial resolution of the binned image in degrees. 
#Fov is 16x27, so 16/xres should give an integer and 27/yres should give an integer. 
xres = 1
yres = 1 

#The code assumes that the raw L3 files are stored in '/data/smile/shared/sims/GAC_V2/' 
#You can change this in the file SXI_L3/paths.py. 

#Put the full name of the output file location. RECOMMEND CHANGING THIS! 
fitspath = './' 

#Run the rebinning code. 
rebin = SXI_L3.rebin_files.rebin_files(stime=start, etime=end, xres=xres, yres=yres, fitspath=fitspath) 

#Create a new fits file of the rebinned file. 
rebin.create_fits_file() 

#PLOT THE CONTENTS OF A REBINNED FILE. 
#By default, it performs all plotting functions. 
#All plots will save to fitspath if save=True. 
#Use vmin and vmax to scale the colourbar. 

if plot: 
    
    #Show the key extensions. 
    rebin.plot_key_extensions(cmap='lundi', vmin=0, vmax=20, save=False)
    
    #Show all the extensions. 
    rebin.plot_all_extensions(cmap='lundi', vmin=0, vmax=20, save=False) 
    
