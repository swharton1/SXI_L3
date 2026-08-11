#This will read in the rebinned file created by rebin_files.py. 
#It should do the same job as read_fits_image.py in CMEM_Ops. 

import os 
import numpy as np 
from astropy.io import fits 
import matplotlib.pyplot as plt 
import datetime as dt

#from . import paths 
from . import read_config
from .SXI_Core import read_cmap 
from .SXI_Core import make_image_axes 

class read_rebinned_file():
    '''This reads in the new rebinned file.''' 
    
    def __init__(self, folder='/home/s/sw682/Code/SXI_L3/SXI_L3/binned_examples/', filename='SMILE_SXI_L3_SCIM60x60-SCI-CXF_20260317T0240-20260317T0244_V01.fits', fitspath=None):
        
        '''This reads in a rebinned file.
        
        Parameters
        ----------
        folder - absolute path to the file. 
        filename - filename. 
        fitspath - absolute path to output files. None - uses default. 
        
        '''
        
        #Get path to the data. 
        self.filename = filename 
        if folder is None: 
            #self.datapath = paths.get_fits_path()
            self.datapath = read_config.read_config(path_type='data_path') 
        else: 
            self.datapath = folder 
        self.fullname = os.path.join(self.datapath, filename) 
        
        if fitspath is None: 
            #self.fitspath = paths.get_fits_path()
            self.fitspath = read_config.read_config(path_type='fits_path')
        else:
            self.fitspath = fitspath 
            
        
        
        #Check the file exists. 
        assert os.path.isfile(self.fullname), f'{self.fullname} does not exist' 
        
        #Now you know it exists, open it. 
        print (f'Read {self.fullname}...')
        with fits.open(self.fullname) as hdul: 
            self.hdul = hdul 
            
            #Read out the entire file. Headers then data. 
            self.primary_header = self.hdul['PRIMARY'].header 
            
            #Works out if it has read a constrained file here. 
            if 'CALC_CX' in self.primary_header:
                self.calc_cxfov = self.primary_header['CALC_CX'] 
            else:
                self.calc_cxfov = 'normal' 
                
            #Now data. 
            self.data = {} 
            self.data['CTSMAP'] = self.hdul['CTSMAP'].data
            self.data['BKGMAP'] = self.hdul['BKGMAP'].data 
            self.data['XBMAP'] = self.hdul['XBMAP'].data 
            self.data['PSMAP'] = self.hdul['PSMAP'].data 
            self.data['PBMAP'] = self.hdul['PBMAP'].data 
            self.data['CLMAP'] = self.hdul['CLMAP'].data 
            self.data['SPMAP'] = self.hdul['SPMAP'].data 
            self.data['VIGMAP'] = self.hdul['VIGMAP'].data 
            self.data['CXFOV'] = self.hdul['CXFOV'].data
            self.data['ERRFOV'] = self.hdul['ERRFOV'].data
            
            if self.calc_cxfov == 'constrained':
                self.data['BKGCON'] = self.hdul['BKGCON'].data 
                self.data['CXFOV_CON'] = self.hdul['CXFOV_CON'].data
                
        #Further specific extractions of data. 
        self.get_orbit_info() 
        self.get_camera_info() 
        
        
    #FUNCTIONS TO EXTRACT KEY HEADER INFO FROM THE FILE, INCLUDING DIPOLE ANGLE.
    ########################################################################
           
    def get_orbit_info(self):
        '''This will extract the spacecraft position, aim point and time.''' 
        
        #Smile location 
        self.pos = np.array([self.primary_header['POS_X'], self.primary_header['POS_Y'], self.primary_header['POS_Z']]) 
        
        #SXI Aim point 
        self.aim = np.array([self.primary_header['AIM_X'], self.primary_header['AIM_Y'], self.primary_header['AIM_Z']])    
        
        #Time  
        self.date_obs = self.primary_header['DATE-OBS'] 
        self.date_end = self.primary_header['DATE-END']
        
        
        #Get datetime object for the start. 
        self.stime = dt.datetime.strptime(self.date_obs, '%Y-%m-%dT%H:%M:%S.%f') 
        self.etime = dt.datetime.strptime(self.date_end, '%Y-%m-%dT%H:%M:%S.%f') 
        
        self.date_obs_str = dt.datetime.strftime(self.stime, '%Y%m%dT%H%M') 
        self.date_end_str = dt.datetime.strftime(self.etime, '%Y%m%dT%H%M') 
        
        
        #Get Energy Bands. 
        self.emin = self.primary_header['EMIN']
        self.emax = self.primary_header['EMAX']
        
        #Get Pointing in Sky coords. 
        self.ra = self.primary_header['RA']
        self.dec = self.primary_header['DEC'] 

    def get_camera_info(self):
        '''This will get information about the camera and resolution, including plotting arrays.''' 
            
        #Get information about the camera out. 
        #Number of pixels. 
        self.m_pixels = self.primary_header['NAXIS1']
        self.n_pixels = self.primary_header['NAXIS2'] 
        
        #Pixel widths. 
        self.xres = self.primary_header['CDELT1']
        self.yres = self.primary_header['CDELT2']
        
        #Lower bounds. 
        self.xdeg_min = self.primary_header['CRVAL1']
        self.ydeg_min = self.primary_header['CRVAL2']
        
        #Units. 
        self.x_unit = self.primary_header['CTYPE1']
        self.y_unit = self.primary_header['CTYPE2'] 
        
        #Calculate FOV. 
        self.phi_fov = -2*self.xdeg_min 
        self.theta_fov = -2*self.ydeg_min 
        
        #Get 1D pixel arrays for plotting. The edges of the pixels.  
        self.xarray = np.linspace(self.xdeg_min, -self.xdeg_min, self.m_pixels+1)
        self.yarray = np.linspace(self.ydeg_min, -self.ydeg_min, self.n_pixels+1)
        
        #Make 2D arrays for x and y. 
        self.X, self.Y = np.meshgrid(self.xarray, self.yarray)
        
        #Exposure 
        self.expos = self.primary_header['EXPOS']   
        
    def plot_key_extensions(self, cmap='lundi', vmin=0, vmax=20, save=False, close=False, per_pixel=False):
        '''This will plot the final most important extensions, CTSMAP, BKGMAP and CXFOV.'''
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   
  
        if close: 
            plt.ioff() 
        else: 
            plt.ion() 
              
        fig = plt.figure(figsize=(8,4))
        fig.subplots_adjust(left=0.1, wspace=0.5)
        
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133) 
        
        #Scale by pixel size if necessary. 
        if per_pixel:
            px_size = self.xres*self.yres 
            scale_size = px_size
            cbar_title = r'Counts/deg$^2$'
        else:
            scale_size = 1
            cbar_title = 'Counts/pixel'
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax1, self.data['CTSMAP']/scale_size, self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=True, add_cbar=True)
        ax1.set_title('CTSMAP\n', fontsize=10) 
        
        if self.calc_cxfov.lower() == 'normal':
        
            make_image_axes.make_image_axes(ax2, self.data['BKGMAP']/scale_size, self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True)
            ax2.set_title('BKGMAP\n', fontsize=10) 
            
            make_image_axes.make_image_axes(ax3, self.data['CXFOV']/scale_size, self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title=cbar_title, ylabel=False, add_cbar=True)
            ax3.set_title('CXFOV\n', fontsize=10) 

            filename = f'SMILE_SXI_L3_SCIM{self.xres*60}x{self.yres*60}-SCI-CXF_{self.date_obs_str}-{self.date_end_str}_V01_key_ext.png'
            
        elif self.calc_cxfov.lower() == 'constrained':
        
            make_image_axes.make_image_axes(ax2, self.data['BKGCON']/scale_size, self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True)
            ax2.set_title('BKGCON\n', fontsize=10) 
            
            make_image_axes.make_image_axes(ax3, self.data['CXFOV']/scale_size, self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title=cbar_title, ylabel=False, add_cbar=True)
            ax3.set_title('CXFOV_CON\n', fontsize=10) 
            
            filename = f'SMILE_SXI_L3_SCIM{self.xres*60}x{self.yres*60}-SCI-CXFCON_{self.date_obs_str}-{self.date_end_str}_V01_key_ext.png'
            
        else:
            raise ValueError("calc_cxfov must be 'normal' or 'constrained'.") 
            
        #Add times to the plot. 
        fig.text(0.5, 0.95, f'{self.date_obs} - {self.date_end}', ha='center', fontsize=10)
        fig.text(0.5, 0.90, f'SMILE = ({self.pos[0]:.2f},{self.pos[1]:.2f},{self.pos[2]:.2f}), Aim = ({self.aim[0]:.2f},{self.aim[1]:.2f},{self.aim[2]:.2f}), Exposure = {self.expos}s', ha='center', fontsize=10)
        

        if save: 
            
            print ('Saving: ', self.fitspath+filename)
            fig.savefig(self.fitspath+filename)

    def plot_all_extensions(self, cmap='lundi', vmin=0, vmax=20, save=False):
        '''This will plot the final most important extensions, CTSMAP, BKGMAP and CXFOV.
        This does not support the constrained versions of background and foreground. 
        
        '''
        
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
        make_image_axes.make_image_axes(ax1, self.data['CTSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=True, add_cbar=True, xlabel=False)
        ax1.set_title('CTSMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax2, self.data['BKGMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True, xlabel=False)
        ax2.set_title('BKGMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax3, self.data['XBMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True, xlabel=False)
        ax3.set_title('XBMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax4, self.data['PSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True, xlabel=False)
        ax4.set_title('PSMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax5, self.data['PBMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='Counts/pixel', ylabel=False, add_cbar=True, xlabel=False)
        ax5.set_title('PBMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax6, self.data['CLMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=True, add_cbar=True, xlabel=True)
        ax6.set_title('CLMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax7, self.data['SPMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True)
        ax7.set_title('SPMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax8, self.data['VIGMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=1, cbar_title='', ylabel=False, add_cbar=True)
        ax8.set_title('VIGMAP\n', fontsize=10) 
        
        make_image_axes.make_image_axes(ax9, self.data['ERRFOV'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='', ylabel=False, add_cbar=True)
        ax9.set_title('ERRFOV\n', fontsize=10)       
        
        make_image_axes.make_image_axes(ax10, self.data['CXFOV'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=vmax, cbar_title='Counts/pixel', ylabel=False, add_cbar=True)
        ax10.set_title('CXFOV\n', fontsize=10)    
        
        if save: 
            filename = f'SMILE_SXI_L3_SCIM{self.xres*60}x{self.yres*60}-SCI-CXF_{self.date_obs_str}-{self.date_end_str}_V01_all_ext.png'
            print ('Saving: ', self.fitspath+filename)
            fig.savefig(self.fitspath+filename)
            
