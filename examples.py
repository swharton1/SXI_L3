#This shows you how to create a rebinned file. 

import SXI_L3 
import matplotlib.pyplot as plt

#To rebin a file. 
rebin_file = SXI_L3.rebin_files.rebin_files(start_folder='sim_0600', end_folder='sim_0603', xres=1, yres=1)

#To save it as a new FITS file. 
rebin_file.create_fits_file() 

#To plot the key extensions of this file. 
rebin_file.plot_key_extensions() 

#To plot all extensions of this file. 
rebin_file.plot_all_extensions() 

#To read this file in... 
read_file = SXI_L3.read_rebinned_file.read_rebinned_file() 

plt.show()
