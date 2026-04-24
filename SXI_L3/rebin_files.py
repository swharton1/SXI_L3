#This will use the new file format and rebin both spatially and temporally. 
#This is the file you should use to create rebinned files for fitting going forwards. 

import os 
import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits 
import datetime as dt 

from . import paths 
from . import read_bkg_file
from . import read_tot_file
from . import read_vcy_file  
from . import rebin_funcs 
from .SXI_Core import make_image_axes
from . import rebin_file_spatially 
from .SXI_Core import read_cmap 
from . import quality_flag as qf 


class rebin_files(): 
    '''This will read in a set number of files and combine them in space and time.''' 
    
    def __init__(self, stime='20260317T0240', etime='20260317T0244', xres=1, yres=1, fitspath=None):
        '''This takes in the start and end folders, gets the filenames and then combines data temporally before combining it spatially. 
        
        Parameters
        ----------
        stime - start time in string format 'YYYYmmddTHHMM'
        etime - end time in string format 'YYYYmmddTHHMM'  
        xres - x angular resolution to bin too. def = 1. Must be greater than that in the original file. 
        yres - y angular resolution to bin too. def = 1. Must be greater than that in the original file.
        fitspath - absolute path for the output file.  
        
        '''
        
        #Get lists of folders and filenames. 
        self.stime = stime
        self.etime = etime 
        self.datapath = paths.get_data_path()
        
        #Use default output location of no fitspath is given. 
        if fitspath is None: 
            self.fitspath = paths.get_fits_path()
        else: 
            self.fitspath = fitspath 
            
        self.get_folders_and_filenames() 
        
        self.xres = xres
        self.yres = yres 
        
        #Now loop through each folder and create a combined object with the spatial rebinning, devignetting and CXFOV calculation done. 
        self.rebin_files = [] 
        for f, fval in enumerate(self.folders): 
            rebin_file = rebin_file_spatially.rebin_file(folder=fval, bkg_file=self.bkg_files[f], tot_file=self.tot_files[f], vcy_file=self.vcy_files[f], xres=xres, yres=yres) 
            self.rebin_files.append(rebin_file) 
        
        #Now combine files temporally. 
        self.rebin_temporally()
        
        #Now calculate CXFOV. 
        print ('Calculating CXFOV...') 
        self.rebin_final['CXFOV'] = (self.rebin_final['CTSMAP'] - self.rebin_final['BKGMAP'])/self.rebin_final['VIGMAP'] 
        
        #Calculate the error in the FOV. 
        #FOR NOW, THIS JUST COMES FROM THE ERROR IN CTSMAP. 
        print ('Estimating Error...') 
        self.rebin_final['ERRFOV'] = np.sqrt(self.rebin_final['CTSMAP'])
        
        #Extract some key information. 
        self.aim = self.rebin_files[0].bkg.aim
        self.pos = self.rebin_files[0].bkg.pos 
        self.n_pixels = self.rebin_files[0].n_pixels
        self.m_pixels = self.rebin_files[0].m_pixels 
        self.xdeg_min = self.rebin_files[0].xdeg_min
        self.ydeg_min = self.rebin_files[0].ydeg_min 
        self.x_unit = self.rebin_files[0].x_unit
        self.y_unit = self.rebin_files[0].y_unit 
        self.ra = self.rebin_files[0].ra
        self.dec = self.rebin_files[0].dec 
        self.emin = self.rebin_files[0].emin
        self.emax = self.rebin_files[0].emax
        self.date_obs = self.rebin_files[0].date_obs
        self.date_end = self.rebin_files[-1].date_end 
        
        self.expos = self.rebin_files[0].expos*len(self.folders) 
        
        #Make datetime strings in SMILE formating convention. 
        date_obs_obj = dt.datetime.strptime(self.date_obs, '%Y-%m-%dT%H:%M:%S.%f') 
        date_end_obj = dt.datetime.strptime(self.date_end, '%Y-%m-%dT%H:%M:%S.%f') 
        self.date_obs_str = dt.datetime.strftime(date_obs_obj, '%Y%m%dT%H%M') 
        self.date_end_str = dt.datetime.strftime(date_end_obj, '%Y%m%dT%H%M') 
    
        #Calculate quality flag. 
        self.qf = qf.calc_quality_flag(self.aim, self.pos, self.expos, self.rebin_final['CXFOV']) 
    

        
    def get_folders_and_filenames(self):
        '''This is the updated way of getting the folder names using the time.'''
        
        print ('Get folder names...') 
        
        stime_dt = dt.datetime.strptime(self.stime, '%Y%m%dT%H%M')
        etime_dt = dt.datetime.strptime(self.etime, '%Y%m%dT%H%M')  
        minute = dt.timedelta(hours=1/60)
        n = int((etime_dt-stime_dt).seconds/60)
        
        #Create list of datetimes. 
        self.times = [] 
        for i in range(n+1):
            self.times.append(stime_dt + i*minute)
            
        #Now loop through these to make folder names. 
        self.folders = []
        self.bkg_files = []
        self.tot_files = []
        self.vcy_files = [] 
        for i in range(n):
            
            start = dt.datetime.strftime(self.times[i], '%Y%m%dT%H%M')
            end = dt.datetime.strftime(self.times[i+1], '%Y%m%dT%H%M')
            
            #Folder names. 
            folder_name = f'L3_{start}-{end}/'
            self.folders.append(folder_name)
            
            #Filenames. 
            self.bkg_files.append(f'SMILE_SXI_L3_SCIM15-SCI-BKG_{start}-{end}_V01.fits') 
            self.tot_files.append(f'SMILE_SXI_L3_SCIM15-SCI-TOT_{start}-{end}_V01.fits') 
            self.vcy_files.append(f'SMILE_SXI_L3_SCIM15-SCI-VCY_{start}-{end}_V01.fits') 
            
           
#    def get_folders_and_filenames(self):
#        '''This will work out the names of all the folders and filenames. It depends on the format of your folders and filenames, so this could change. ''' 
        
#        print ('Get folder and filenames...') 
#        #Get start and end numbers. 
#        self.start_num = int(self.start_folder[-4:])
#        self.end_num = int(self.end_folder[-4:]) 
        
        #Make list. 
#        self.folder_nums = np.arange(self.start_num, self.end_num+1)
#        self.folders = [f'sim_{n:0>4}/' for n in self.folder_nums]  
        
        #Prepare empty lists of filenames. 
#        self.bkg_files = []
#        self.tot_files = []
#        self.vcy_files = []  
        
        #Go through each folder and pull out the name. 
#        for f in self.folders:
            
#            self.bkg_files.append(glob.glob1(self.datapath+f, '*-BKG_*')[0])
#            self.tot_files.append(glob.glob1(self.datapath+f, '*-TOT_*')[0])
#            self.vcy_files.append(glob.glob1(self.datapath+f, '*-VCY_*') [0])


    def rebin_temporally(self):
        '''This will add the data up from each file to get total counts over the full time period. It currently assumes that the vignetting map doesn't change over the time period and uses it from the start. '''
        
        print ('Rebinning data temporally...') 
        
        #Now you need to combine data from each file together (except VIGMAP). 
        self.rebin_final = {} 
        self.rebin_final['CTSMAP'] = np.zeros((self.rebin_files[0].rebin_data['CTSMAP'].shape))
        self.rebin_final['BKGMAP'] = np.zeros((self.rebin_files[0].rebin_data['BKGMAP'].shape))
        self.rebin_final['XBMAP'] = np.zeros((self.rebin_files[0].rebin_data['XBMAP'].shape))
        self.rebin_final['PSMAP'] = np.zeros((self.rebin_files[0].rebin_data['PSMAP'].shape))
        self.rebin_final['PBMAP'] = np.zeros((self.rebin_files[0].rebin_data['PBMAP'].shape))
        self.rebin_final['CLMAP'] = np.zeros((self.rebin_files[0].rebin_data['CLMAP'].shape))
        self.rebin_final['SPMAP'] = np.zeros((self.rebin_files[0].rebin_data['SPMAP'].shape))
        
        for n in range(len(self.rebin_files)):
            self.rebin_final['CTSMAP'] += self.rebin_files[n].rebin_data['CTSMAP']
            self.rebin_final['BKGMAP'] += self.rebin_files[n].rebin_data['BKGMAP'] 
            self.rebin_final['XBMAP'] += self.rebin_files[n].rebin_data['XBMAP']
            self.rebin_final['PSMAP'] += self.rebin_files[n].rebin_data['PSMAP']
            self.rebin_final['PBMAP'] += self.rebin_files[n].rebin_data['PBMAP']
            self.rebin_final['CLMAP'] += self.rebin_files[n].rebin_data['CLMAP']
            self.rebin_final['SPMAP'] += self.rebin_files[n].rebin_data['SPMAP']
            
        
        #Get the vignetting map on the assumption it is the same for all files. 
        #Choose the first one. 
        self.rebin_final['VIGMAP'] = self.rebin_files[0].rebin_data['VIGMAP']         

    def plot_ctsmap_all(self, cmap='lundi', vmin=0, vmax=80, save=False):
        '''This will make a plot showing how all the CTSMAP extensions are added together to produce one with a single, larger exposure time.''' 

        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   
  
          
        fig = plt.figure(figsize=(8,3))
        fig.subplots_adjust(left=0.05, wspace=0.5)
        
        #How many axes do you need.
        n = len(self.folders)
        
        if n > 4: 
            print ("Too many axes to plot. Only showing first four.") 
            n = 4 
        
        #Plot individual frames. 
        for i in range(n):
            ax = fig.add_subplot(1, n+1, i+1)
                
            #Make the axis. 
            make_image_axes.make_image_axes(ax, self.rebin_files[i].rebin_data['CTSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=False)   
            
        
        axt = fig.add_subplot(1, n+1, n+1)
        
        #Make the axis. 
        make_image_axes.make_image_axes(axt, self.rebin_final['CTSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='Counts/pixel', ylabel=False, add_cbar=False)
        axt.set_title('Combined\n', fontsize=10) 
        
        fig.text(0.5, 0.9, 'Integrated CTSMAP', ha='center')
        
        if save: 
            filename = f'SMILE_SXI_L3_SCIM{self.xres*60}x{self.yres*60}-SCI-CTS_{self.date_obs_str}-{self.date_end_str}_V01_integrated.png'
            print ('Saving: ', self.fitspath+filename)
            fig.savefig(self.fitspath+filename)
            
    def plot_key_extensions(self, cmap='lundi', vmin=0, vmax=20, save=False):
        '''This will plot the final most important extensions, CTSMAP, BKGMAP and CXFOV.'''
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   
  
          
        fig = plt.figure(figsize=(8,4))
        fig.subplots_adjust(left=0.1, wspace=0.5)
        
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133) 
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax1, self.rebin_final['CTSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=True, add_cbar=True)
        ax1.set_title('CTSMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax2, self.rebin_final['BKGMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True)
        ax2.set_title('BKGMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax3, self.rebin_final['CXFOV'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='Counts/pixel', ylabel=False, add_cbar=True)
        ax3.set_title('CXFOV\n', fontsize=10) 

        if save: 
            filename = f'SMILE_SXI_L3_SCIM{self.xres*60}x{self.yres*60}-SCI-CXF_{self.date_obs_str}-{self.date_end_str}_V01_key_ext.png'
            print ('Saving: ', self.fitspath+filename)
            fig.savefig(self.fitspath+filename)

    def plot_all_extensions(self, cmap='lundi', vmin=0, vmax=20, save=False):
        '''This will plot the final most important extensions, CTSMAP, BKGMAP and CXFOV.'''
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   
  
          
        fig = plt.figure(figsize=(8,4))
        fig.subplots_adjust(left=0.1, wspace=0.7, hspace=0.4)
        
        ax1 = fig.add_subplot(251)
        ax2 = fig.add_subplot(252)
        ax3 = fig.add_subplot(253) 
        ax4 = fig.add_subplot(254)
        ax5 = fig.add_subplot(255)
        ax6 = fig.add_subplot(256) 
        ax7 = fig.add_subplot(257)
        ax8 = fig.add_subplot(258)
        ax9 = fig.add_subplot(259)
        ax10 = fig.add_subplot(2,5,10)
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax1, self.rebin_final['CTSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=True, add_cbar=True, xlabel=False)
        ax1.set_title('CTSMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax2, self.rebin_final['BKGMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True, xlabel=False)
        ax2.set_title('BKGMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax3, self.rebin_final['XBMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True, xlabel=False)
        ax3.set_title('XBMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax4, self.rebin_final['PSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True, xlabel=False)
        ax4.set_title('PSMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax5, self.rebin_final['PBMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='Counts/pixel', ylabel=False, add_cbar=True, xlabel=False)
        ax5.set_title('PBMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax6, self.rebin_final['CLMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=True, add_cbar=True, xlabel=True)
        ax6.set_title('CLMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax7, self.rebin_final['SPMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True)
        ax7.set_title('SPMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax8, self.rebin_final['VIGMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=1, cbar_title='', ylabel=False, add_cbar=True)
        ax8.set_title('VIGMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax9, self.rebin_final['ERRFOV'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True)
        ax9.set_title('ERRFOV\n', fontsize=10)       
        
        make_image_axes.make_image_axes(ax10, self.rebin_final['CXFOV'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='Counts/pixel', ylabel=False, add_cbar=True)
        ax10.set_title('CXFOV\n', fontsize=10)    
        
        if save: 
            filename = f'SMILE_SXI_L3_SCIM{self.xres*60}x{self.yres*60}-SCI-CXF_{self.date_obs_str}-{self.date_end_str}_V01_all_ext.png'
            print ('Saving: ', self.fitspath+filename)
            fig.savefig(self.fitspath+filename)
            
            
    def create_fits_file(self):
        '''This will create a FITS file with similar extensions to that produced by the SXI simulator. It should be readable by my read_fits_image.py script.''' 
        
        #Create a filename. 
        filename = f'SMILE_SXI_L3_SCIM{self.xres*60}x{self.yres*60}-SCI-CXF_{self.date_obs_str}-{self.date_end_str}_V01.fits'
        self.outname=self.fitspath+filename 
        
        #Create a new HDU object. 
        self.hdu = fits.PrimaryHDU(data=self.rebin_final['CXFOV'])
        self.hdu.header
        
        #Add primary header info here.  
        self.hdu.header['POS_X'] = self.pos[0]
        self.hdu.header['POS_Y'] = self.pos[1]
        self.hdu.header['POS_Z'] = self.pos[2]
        self.hdu.header['AIM_X'] = self.aim[0]
        self.hdu.header['AIM_Y'] = self.aim[1]
        self.hdu.header['AIM_Z'] = self.aim[2]
        
        #Add the other extensions. Start with CTSMAP and adapt its header.  
        #self.hdu1 = fits.ImageHDU(data=self.data['CXFOV'], name='CXFOV')
        self.hdu.header['CTYPE1'] = self.x_unit
        self.hdu.header['CRVAL1'] = self.xdeg_min 
        self.hdu.header['CDELT1'] = self.xres 
        self.hdu.header['CTYPE2'] = self.y_unit
        self.hdu.header['CRVAL2'] = self.ydeg_min 
        self.hdu.header['CDELT2'] = self.yres 
        self.hdu.header['EXPOS'] = self.expos
        self.hdu.header['EMIN'] = self.emin
        self.hdu.header['EMAX'] = self.emax
        self.hdu.header['RA'] = self.ra
        self.hdu.header['DEC'] = self.dec 
        self.hdu.header['DATE-OBS'] = self.date_obs
        self.hdu.header['DATE-END'] = self.date_end 
        self.hdu.header['COMMENT'] = 'Made by rebin_files.py' 
        self.hdu.header['QF'] = self.qf 
        
        #Add comment to primary file stating which files were used to create this one. 
        for f in self.folders:
            self.hdu.header['COMMENT'] = f'Derived from {f}'
            
        header = self.hdu.header 
        
        #Add the other extensions. 
        self.hdu1 = fits.ImageHDU(data=self.rebin_final['CTSMAP'], name='CTSMAP', header=header)
        self.hdu2 = fits.ImageHDU(data=self.rebin_final['BKGMAP'], name='BKGMAP', header=header)
        self.hdu3 = fits.ImageHDU(data=self.rebin_final['XBMAP'], name='XBMAP', header=header) 
        self.hdu4 = fits.ImageHDU(data=self.rebin_final['PSMAP'], name='PSMAP', header=header) 
        self.hdu5 = fits.ImageHDU(data=self.rebin_final['PBMAP'], name='PBMAP', header=header) 
        self.hdu6 = fits.ImageHDU(data=self.rebin_final['CLMAP'], name='CLMAP', header=header) 
        self.hdu7 = fits.ImageHDU(data=self.rebin_final['SPMAP'], name='SPMAP', header=header) 
        self.hdu8 = fits.ImageHDU(data=self.rebin_final['VIGMAP'], name='VIGMAP', header=header) 
        self.hdu9 = fits.ImageHDU(data=self.rebin_final['CXFOV'], name='CXFOV', header=header)
        self.hdu10 = fits.ImageHDU(data=self.rebin_final['ERRFOV'], name='ERRFOV', header=header) 
        
        #Add Original Comments. 
        self.hdu1.header['COMMENT'] = 'Observed Total Counts Map' 
        self.hdu2.header['COMMENT'] = 'Predicted Total Background Counts Map'
        self.hdu3.header['COMMENT'] = 'Predicted diffuse X-ray Map' 
        self.hdu4.header['COMMENT'] = 'Predicted Point Source Map' 
        self.hdu5.header['COMMENT'] = 'Predicted Particle Background Map' 
        self.hdu6.header['COMMENT'] = 'Predicted Calibration Source Map' 
        self.hdu7.header['COMMENT'] = 'Predicted Soft Proton Map' 
        self.hdu8.header['COMMENT'] = 'Predicted Relative Exposure Map (Vignetting)' 
        self.hdu9.header['COMMENT'] = 'Foreground, Devignetted Count Map (CXFOV)'
        self.hdu10.header['COMMENT'] = 'Error in counts. Assumed square root of CTSMAP'
        
        #Make the HDU list.  
        self.hdul = fits.HDUList([self.hdu, self.hdu1, self.hdu2, self.hdu3, self.hdu4, self.hdu5, self.hdu6, self.hdu7, self.hdu8, self.hdu9, self.hdu10]) 
        
        
        #Write the fits file. 
        self.hdul.writeto(self.outname, overwrite=True) 
        print ('Created: {}'.format(self.outname))      
