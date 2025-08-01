#This file will read in all files within a given timeframe and bin both temporally and spatially. 

#How it gathers the files depends on their location and naming format, so this code may have to change. 

import os 
import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits 

from . import paths 
from . import read_mapL3 
from . import rebin_funcs 
from . import read_cmap 
from . import make_image_axes


class read_cxfov_time():
    '''This will read in a set number of files and combine them together.''' 
    
    def __init__(self, start_folder='sim_0600', end_folder='sim_0601', xres=1, yres=1):
        '''This takes in the start and end folders, gets the filenames and then combines data temporally before combining it spatially. 
        
        Parameters
        ----------
        start_folder - First folder to include. 
        end_folder - Last fodler to include. 
        xres - x angular resolution to bin too. def = 1. Must be greater than that in the original file. 
        yres - y angular resolution to bin too. def = 1. Must be greater than that in the original file. 
        
        '''
        
        #Get lists of folders and filenames. 
        self.start_folder = start_folder
        self.end_folder = end_folder 
        self.get_folders_and_filenames() 
        
        #Now read in all of the mapL3 files. 
        self.mapL3 = [] 
        for n in range(len(self.folders)):
            self.mapL3.append(read_mapL3.read_mapL3(folder=self.folders[n], filename=self.filenames[n]))    

        #Bin temporally. 
        self.rebin_temporally() 

        #Bin spatially. 
        self.xres = xres
        self.yres = yres 
        self.rebin_arrays() 

        #Calculate CXFOV. 
        self.calc_cxfov()
        
        #Pull out any other meta information you might need. 
        self.xdeg_min = self.mapL3[0].xdeg_min 
        self.ydeg_min = self.mapL3[0].ydeg_min 
        self.date_obs = self.mapL3[0].date_obs
        self.date_end = self.mapL3[-1].date_end 
        self.pos = self.mapL3[0].pos
        self.aim = self.mapL3[0].aim 
        
        #Calculate the combined exposure time. 
        self.expos = 0 
        for m in self.mapL3:
            self.expos += m.expos 


    def get_folders_and_filenames(self):
        '''This will work out the names of all the folders and filenames. It depends on the format of your folders and filenames, so this could change. ''' 
        
        #Get start and end numbers. 
        self.start_num = int(self.start_folder[-4:])
        self.end_num = int(self.end_folder[-4:]) 
        
        #Make list. 
        self.file_list_nums = [] 
        num = self.start_num 
        while num <= self.end_num:
            self.file_list_nums.append(num)
            num += 1 
        
        self.folders = [f'sim_{n:0>4}/' for n in self.file_list_nums]  
        self.filenames = [f'mapL3_{n:0>4}.fits' for n in self.file_list_nums]    
    
    def calc_cxfov(self):
        '''This will calculate the CXFOV with the rebinned arrays.
        Must execute CXFOV array first.
        
        THIS IS CURRENTLY WRONG. ''' 
        
        #Add all the background contributions. 
        background = self.data['XBMAP']+self.data['PSMAP']+self.data['PBMAP']+self.data['CLMAP']+self.data['SPMAP']
        
        #Calculate CXFOV.
        self.data['CXFOV'] = (self.data['CTSMAP'] - background)/self.data['VIGMAP']  
            
    #FUNCTIONS FOR REBINNING IN TIME AND SPATIALLY. 
    ###############################################
    
    def rebin_temporally(self):
        '''This will add the data up from each file to get total counts over the full time period. It currently assumes that the vignetting map doesn't change over the time period and uses it from the start. '''
        
        print ('Rebinning data temporally...') 
        
        #Now you need to combine data from each file together (except VIGMAP). 
        self.data_t = {} 
        self.data_t['CTSMAP'] = np.zeros((self.mapL3[0].data['CTSMAP'].shape))
        self.data_t['XBMAP'] = np.zeros((self.mapL3[0].data['XBMAP'].shape))
        self.data_t['PSMAP'] = np.zeros((self.mapL3[0].data['PSMAP'].shape))
        self.data_t['PBMAP'] = np.zeros((self.mapL3[0].data['PBMAP'].shape))
        self.data_t['CLMAP'] = np.zeros((self.mapL3[0].data['CLMAP'].shape))
        self.data_t['SPMAP'] = np.zeros((self.mapL3[0].data['SPMAP'].shape))
        
        for n in range(len(self.mapL3)):
            self.data_t['CTSMAP'] += self.mapL3[n].data['CTSMAP']
            self.data_t['XBMAP'] += self.mapL3[n].data['XBMAP']
            self.data_t['PSMAP'] += self.mapL3[n].data['PSMAP']
            self.data_t['PBMAP'] += self.mapL3[n].data['PBMAP']
            self.data_t['CLMAP'] += self.mapL3[n].data['CLMAP']
            self.data_t['SPMAP'] += self.mapL3[n].data['SPMAP']
        
        #Get the vignetting map on the assumption it is the same for all files. 
        self.data_t['VIGMAP'] = self.mapL3[0].data['VIGMAP']     
        
    def rebin_arrays(self):
        '''This will rebin all the arrays spatially and add them to a new dictionary.''' 
        
        print ('Rebinning data spatially...') 
        
        #Check xres and yres are appropriate choices. 
        assert self.xres >= self.mapL3[0].xdeg_sep, f'xres must be greater than {self.mapL3.xdeg_sep}'
        assert self.yres >= self.mapL3[0].ydeg_sep, f'yres must be greater than {self.mapL3.ydeg_sep}'
        
        #Get new shape for the arrays. 
        self.m_pixels = self.mapL3[0].phi_fov/self.xres 
        self.n_pixels = self.mapL3[0].theta_fov/self.yres 
        
        #Check this gives an integer shape. 
        assert self.m_pixels.is_integer(), f'xres = {self.xres} did not give an integer number of pixels for m_pixels. Try something else.' 
        assert self.n_pixels.is_integer(), f'yres = {self.yres} did not give an integer number of pixels for n_pixels. Try something else.' 
        self.m_pixels = int(self.m_pixels)
        self.n_pixels = int(self.n_pixels) 
        
        #Now get the new shape. 
        self.new_shape = (self.n_pixels, self.m_pixels)  
        
        #Now start rebinning all the data products.  
        self.data = {}  
        self.data['CTSMAP'] = rebin_funcs.rebin_sum(self.data_t['CTSMAP'], self.new_shape) 
        self.data['XBMAP'] = rebin_funcs.rebin_sum(self.data_t['XBMAP'], self.new_shape) 
        self.data['PSMAP'] = rebin_funcs.rebin_sum(self.data_t['PSMAP'], self.new_shape) 
        self.data['PBMAP'] = rebin_funcs.rebin_sum(self.data_t['PBMAP'], self.new_shape) 
        self.data['CLMAP'] = rebin_funcs.rebin_sum(self.data_t['CLMAP'], self.new_shape) 
        self.data['SPMAP'] = rebin_funcs.rebin_sum(self.data_t['SPMAP'], self.new_shape) 
        self.data['VIGMAP'] = rebin_funcs.rebin_mean(self.data_t['VIGMAP'], self.new_shape) 
        
        
        

    #PLOTTING FUNCTIONS. 
    ###################
        
    def plot_binned_extension(self, ext = 'CTSMAP', cmap='lundi', vmin=0, vmax=10):
        '''This will plot one of the extensions for you.''' 
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   

        #Use a scale from 0-1 for the vignetting map. Overwrite. 
        if ext.upper() == 'VIGMAP': 
            vmin = 0
            vmax = 1 
                        
        #Create the figure. 
        fig = plt.figure(figsize=(4,6))
        fig.subplots_adjust(top=0.8, left=0.20, right=0.85)
        ax = fig.add_subplot(111)
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax, self.data[ext.upper()], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title='Counts/pixel')
        ax.set_title(ext.upper()+'\n\n')
        
        #Add figure title with key meta information. 
        time_title = f'{self.date_obs}'
        exp_title = f'Expos: {self.expos}s'
        pos_title = 'SMILE: ({:.2f},{:.2f},{:.2f})'.format(*self.pos) 
        aim_title = 'AIM: ({:.2f},{:.2f},{:.2f})'.format(*self.aim) 
        metatitle = time_title+'\n'+exp_title+'\n'+pos_title+'\n'+aim_title 
        
        fig.text(0.15, 0.95, metatitle, ha='left', va='top')         


    def create_fits_file(self):
        '''This will create a FITS file with similar extensions to that produced by the SXI simulator. It should be readable by my read_fits_image.py script.''' 
        
        #Create a filename. 
        self.fitspath = paths.get_fits_path()
        filename = f'cts_{self.folders[0][:-1]}_{self.folders[-1][:-1]}_{self.xres}x{self.yres}.fits'
        self.outname=self.fitspath+filename 
        
        #Create a new HDU object. 
        self.hdu = fits.PrimaryHDU()
        self.hdu.header
        
        #Add primary header info here.  
        self.hdu.header['POS_X'] = self.pos[0]
        self.hdu.header['POS_Y'] = self.pos[1]
        self.hdu.header['POS_Z'] = self.pos[2]
        self.hdu.header['AIM_X'] = self.aim[0]
        self.hdu.header['AIM_Y'] = self.aim[1]
        self.hdu.header['AIM_Z'] = self.aim[2]
        
        #Add the other extensions. Start with CTSMAP and adapt its header.  
        self.hdu1 = fits.ImageHDU(data=self.data['CTSMAP'], name='CTSMAP')
        self.hdu1.header['CTYPE1'] = self.mapL3[0].x_unit
        self.hdu1.header['CRVAL1'] = self.xdeg_min 
        self.hdu1.header['CDELT1'] = self.xres 
        self.hdu1.header['CTYPE2'] = self.mapL3[0].y_unit
        self.hdu1.header['CRVAL2'] = self.ydeg_min 
        self.hdu1.header['CDELT2'] = self.yres 
        self.hdu1.header['EXPOS'] = self.expos
        self.hdu1.header['EMIN'] = self.mapL3[0].emin
        self.hdu1.header['EMAX'] = self.mapL3[0].emax
        self.hdu1.header['RA'] = self.mapL3[0].ra
        self.hdu1.header['DEC'] = self.mapL3[0].dec 
        self.hdu1.header['DATE-OBS'] = self.date_obs
        self.hdu1.header['DATE-END'] = self.date_end 
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
        for f in self.folders:
            self.hdu.header['COMMENT'] = f'Derived from {f}'
        self.hdu.header['COMMENT'] = 'Made by read_cxfov_time.py' 
        #Make the HDU list.  
        self.hdul = fits.HDUList([self.hdu, self.hdu1, self.hdu2, self.hdu3, self.hdu4, self.hdu5, self.hdu6, self.hdu7]) 
        
        
        
        #Write the fits file. 
        self.hdul.writeto(self.outname, overwrite=True) 
        print ('Created: {}'.format(self.outname))     
