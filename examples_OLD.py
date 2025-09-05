#These are standard examples of how you might run the code. 

import SXI_L3 

#EXAMPLE 1
#To read in a single file and plot one of the extensions. 
l3_file = SXI_L3.read_mapL3.read_mapL3(folder='sim_0600/', filename='mapL3_0600.fits') 
l3_file.plot_raw_extension(ext='CTSMAP', cmap='lundi', vmax=20) 



#EXAMPLE 2
#To read in a single file, rebin it spatially and calculate CXFOV. 
#xres and yres are the angular resolution in degrees. Recommend 0.25, 0.5 or 1. 
cxfov = SXI_L3.read_cxfov.read_cxfov(folder='sim_0600/', filename='mapL3_0600.fits', xres=1, yres=1) 

#To plot the original vs the binned data. CXFOV doesn't have an original. 
cxfov.plot_raw_extension_vs_binned(ext='CTSMAP', cmap='lundi', vmax=20)
cxfov.plot_raw_extension_vs_binned(ext='CXFOV', cmap='lundi', vmax=20) 

#To create an output fits file. 
cxfov.create_fits_file() 


#EXAMPLE 3
#Reads in a set of files between the first and last stated. 
#Combines the data in time, then rebins spatially as before. 
cxfovt = SXI_L3.read_cxfov_time.read_cxfov_time(start_folder='sim_0600', end_folder='sim_0604', xres=1, yres=1) 

#To plot one of the extensions of the combined data. 
cxfovt.plot_binned_extension(ext='CXFOV', cmap='lundi', vmax=20) 

#To create an output fits file (same format as example 2). 
cxfovt.create_fits_file() 
