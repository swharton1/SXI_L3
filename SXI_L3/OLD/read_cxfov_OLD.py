#This file will read in the CXFOV information, doing all the binning and calculation, from the original mapL3 file. This only bins spatially. It does not combine files together to bin in time as well. 

import os 
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits 

from . import paths 
from . import read_mapL3_OLD 
from . import rebin_funcs 
from . import read_cmap 
from . import make_image_axes


class read_cxfov():
    '''This class will read in the mapL3 file, rebin the data spatially to the desired size and calculate the CXFOV data.''' 
    
    def __init__(self, folder='sim_0600/', filename='mapL3_0600.fits', xres=1, yres=1):
        '''This reads in the file and constructs the filename. 
        
        Parameters
        ----------
        folder - name of folder. Currently in format 'sim_HHMM' 
        filename - Currently in format 'mapL3_HHMM' 
        xres - x angular resolution to bin too. def = 1. Must be greater than that in the original file. 
        yres - y angular resolution to bin too. def = 1. Must be greater than that in the original file. 
        
        '''
        
        self.filename = filename 
        self.folder = folder 
              
        #Read in the file using read_mapL3.py 
        self.mapL3 = read_mapL3_OLD.read_mapL3(folder=folder, filename=filename)        
        
        #Get the binning resolution data and do basic checks. 
        self.xres = xres 
        self.yres = yres 
        
        #Rebin the arrays. 
        self.rebin_arrays() 
        
        #Calculate CXFOV. 
        self.calc_cxfov() 

    
    def calc_cxfov(self):
        '''This will calculate the CXFOV with the rebinned arrays.
        Must execute CXFOV array first.
        
        THIS IS CURRENTLY WRONG. ''' 
        
        #Add all the background contributions. 
        background = self.data['XBMAP']+self.data['PSMAP']+self.data['PBMAP']+self.data['CLMAP']+self.data['SPMAP']
        
        #Calculate CXFOV.
        self.data['CXFOV'] = (self.data['CTSMAP'] - background)/self.data['VIGMAP']  
        
        
        
    
    
    def rebin_arrays(self):
        '''This will rebin all the arrays and add them to a new dictionary.''' 
        
        print ('Rebinning data...') 
        
        #Check xres and yres are appropriate choices. 
        assert self.xres >= self.mapL3.xdeg_sep, f'xres must be greater than {self.mapL3.xdeg_sep}'
        assert self.yres >= self.mapL3.ydeg_sep, f'yres must be greater than {self.mapL3.ydeg_sep}'
        
        #Get new shape for the arrays. 
        self.m_pixels = self.mapL3.phi_fov/self.xres 
        self.n_pixels = self.mapL3.theta_fov/self.yres 
        
        #Check this gives an integer shape. 
        assert self.m_pixels.is_integer(), f'xres = {self.xres} did not give an integer number of pixels for m_pixels. Try something else.' 
        assert self.n_pixels.is_integer(), f'yres = {self.yres} did not give an integer number of pixels for n_pixels. Try something else.' 
        self.m_pixels = int(self.m_pixels)
        self.n_pixels = int(self.n_pixels) 
        
        #Now get the new shape. 
        self.new_shape = (self.n_pixels, self.m_pixels)  
        
        #Now start rebinning all the data products.  
        self.data = {}  
        self.data['CTSMAP'] = rebin_funcs.rebin_sum(self.mapL3.data['CTSMAP'], self.new_shape) 
        self.data['XBMAP'] = rebin_funcs.rebin_sum(self.mapL3.data['XBMAP'], self.new_shape) 
        self.data['PSMAP'] = rebin_funcs.rebin_sum(self.mapL3.data['PSMAP'], self.new_shape) 
        self.data['PBMAP'] = rebin_funcs.rebin_sum(self.mapL3.data['PBMAP'], self.new_shape) 
        self.data['CLMAP'] = rebin_funcs.rebin_sum(self.mapL3.data['CLMAP'], self.new_shape) 
        self.data['SPMAP'] = rebin_funcs.rebin_sum(self.mapL3.data['SPMAP'], self.new_shape) 
        self.data['VIGMAP'] = rebin_funcs.rebin_mean(self.mapL3.data['VIGMAP'], self.new_shape) 
        
                

    #PLOTTING FUNCTIONS. 
    ###################
        
    def plot_raw_extension_vs_binned(self, ext = 'CTSMAP', cmap='lundi', vmin=0, vmax=10):
        '''This will plot the original extension against the spatially rebinned extension. Use only for when spatial rebinning has taken place. 
        
        ''' 
        plt.close("all")
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   
        
        #Use a scale from 0-1 for the vignetting map. Overwrite. 
        if ext.upper() == 'VIGMAP': 
            vmin = 0
            vmax = 1 
            
        #Create the figure. 
        fig = plt.figure(figsize=(8,6))
        fig.subplots_adjust(top=0.8, left=0.10, right=0.85, wspace=0.4)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122) 
        
        #Plot the original on the left. 
        #Make the axis. 
        if ext.upper() in ['CTSMAP', 'XBMAP', 'PSMAP', 'PBMAP', 'CLMAP', 'SPMAP', 'VIGMAP']: 
            make_image_axes.make_image_axes(ax1, self.mapL3.data[ext.upper()], self.mapL3.xdeg_min, self.mapL3.ydeg_min, self.mapL3.n_pixels, self.mapL3.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title='Counts/pixel')
            ax1.set_title(ext.upper()+'\n\n')
        else:
            ax1.text(0.5, 0.5, 'None', ha='center', va='center', transform=ax1.transAxes) 
        #Plot the original on the left. 
        #Make the axis. 
        make_image_axes.make_image_axes(ax2, self.data[ext.upper()], self.mapL3.xdeg_min, self.mapL3.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title='Counts/pixel')
        ax2.set_title(ext.upper()+'-Rebinned\n\n')
        
        
        #Add figure title with key meta information. 
        time_title = f'{self.mapL3.date_obs}'
        exp_title = f'Expos: {self.mapL3.expos}s'
        pos_title = 'SMILE: ({:.2f},{:.2f},{:.2f})'.format(*self.mapL3.pos) 
        aim_title = 'AIM: ({:.2f},{:.2f},{:.2f})'.format(*self.mapL3.aim) 
        metatitle = time_title+'\n'+exp_title+'\n'+pos_title+'\n'+aim_title 
        
        fig.text(0.15, 0.95, metatitle, ha='left', va='top') 
        
    def create_fits_file(self):
        '''This will create a FITS file with similar extensions to that produced by the SXI simulator. It should be readable by my read_fits_image.py script.''' 
        
        #Create a filename. 
        self.fitspath = paths.get_fits_path()
        filename = f'cts_{self.folder[:-1]}_{self.folder[:-1]}_{self.xres}x{self.yres}.fits'
        self.outname=self.fitspath+filename 
        
        #Create a new HDU object. 
        self.hdu = fits.PrimaryHDU()
        self.hdu.header
        
        #Add primary header info here.  
        self.hdu.header['POS_X'] = self.mapL3.pos[0]
        self.hdu.header['POS_Y'] = self.mapL3.pos[1]
        self.hdu.header['POS_Z'] = self.mapL3.pos[2]
        self.hdu.header['AIM_X'] = self.mapL3.aim[0]
        self.hdu.header['AIM_Y'] = self.mapL3.aim[1]
        self.hdu.header['AIM_Z'] = self.mapL3.aim[2]
        
        #Add the other extensions. Start with CTSMAP and adapt its header.  
        self.hdu1 = fits.ImageHDU(data=self.data['CTSMAP'], name='CTSMAP')
        self.hdu1.header['CTYPE1'] = self.mapL3.x_unit
        self.hdu1.header['CRVAL1'] = self.mapL3.xdeg_min 
        self.hdu1.header['CDELT1'] = self.xres 
        self.hdu1.header['CTYPE2'] = self.mapL3.y_unit
        self.hdu1.header['CRVAL2'] = self.mapL3.ydeg_min 
        self.hdu1.header['CDELT2'] = self.yres 
        self.hdu1.header['EXPOS'] = self.mapL3.expos
        self.hdu1.header['EMIN'] = self.mapL3.emin
        self.hdu1.header['EMAX'] = self.mapL3.emax
        self.hdu1.header['RA'] = self.mapL3.ra
        self.hdu1.header['DEC'] = self.mapL3.dec 
        self.hdu1.header['DATE-OBS'] = self.mapL3.date_obs
        self.hdu1.header['DATE-END'] = self.mapL3.date_end 
        header = self.hdu1.header 
        
        #Add the other extensions. 
        self.hdu2 = fits.ImageHDU(data=self.data['XBMAP'], name='XBMAP', header=header) 
        self.hdu3 = fits.ImageHDU(data=self.data['PSMAP'], name='PSMAP', header=header) 
        self.hdu4 = fits.ImageHDU(data=self.data['PBMAP'], name='PBMAP', header=header) 
        self.hdu5 = fits.ImageHDU(data=self.data['CLMAP'], name='CLMAP', header=header) 
        self.hdu6 = fits.ImageHDU(data=self.data['SPMAP'], name='SPMAP', header=header) 
        self.hdu7 = fits.ImageHDU(data=self.data['VIGMAP'], name='VIGMAP', header=header) 

        #Add Original Comments. 
        self.hdu1.header['COMMENT'] = 'Observed Total Counts Map' 
        self.hdu2.header['COMMENT'] = 'Predicted diffuse X-ray Map' 
        self.hdu3.header['COMMENT'] = 'Predicted Point Source Map' 
        self.hdu4.header['COMMENT'] = 'Predicted Particle Background Map' 
        self.hdu5.header['COMMENT'] = 'Predicted Calibration Source Map' 
        self.hdu6.header['COMMENT'] = 'Predicted Soft Proton Map' 
        self.hdu7.header['COMMENT'] = 'Predicted Relative Exposure Map' 
        
        #Add comment to primary file stating which files were used to create this one. 
        self.hdu.header['COMMENT'] = f'Derived from {self.folder}'
        self.hdu.header['COMMENT'] = 'Made by read_cxfov.py' 
        #Make the HDU list.  
        self.hdul = fits.HDUList([self.hdu, self.hdu1, self.hdu2, self.hdu3, self.hdu4, self.hdu5, self.hdu6, self.hdu7]) 
        
        
        
        #Write the fits file. 
        self.hdul.writeto(self.outname, overwrite=True) 
        print ('Created: {}'.format(self.outname))     
                
        
        
        
              
