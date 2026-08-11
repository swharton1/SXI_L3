#This will read in the total counts file. 

import os 
import numpy as np 
from astropy.io import fits 
import matplotlib.pyplot as plt 
import datetime as dt

#from . import paths 
from . import read_config 
from .SXI_Core import read_cmap 
from .SXI_Core import make_image_axes 

class read_tot_file():
    '''This will read in the background file and all its extensions.''' 
    
    def __init__(self, folder='L3_20260317T0240-20260317T0241/', filename='SMILE_SXI_L3_SCIM15-SCI-TOT_20260317T0240-20260317T0241_V01.fits'):
    
        #Get path to the data. 
        self.filename = filename 
        #self.datapath = paths.get_data_path()+folder 
        self.datapath = read_config.read_config(path_type='data_path')+folder
        #self.fitspath = paths.get_fits_path()
        self.fitspath = read_config.read_config(path_type='fits_path') 
        self.fullname = os.path.join(self.datapath, filename) 
        
        #Check the file exists. 
        assert os.path.isfile(self.fullname), f'{self.fullname} does not exist' 
        
        #Now you know it exists, open it. 
        print (f'Read {self.fullname}...')
        with fits.open(self.fullname) as hdul: 
            self.hdul = hdul 
            
            #Read out the entire file. Headers then data. 
            self.primary_header = self.hdul['PRIMARY'].header 
            #self.ctsmap_header = self.hdul['CTSMAP'].header
            #self.xbmap_header = self.hdul['XBMAP'].header
            #self.psmap_header = self.hdul['PSMAP'].header
            #self.pbmap_header = self.hdul['PBMAP'].header
            #self.clmap_header = self.hdul['CLMAP'].header
            #self.spmap_header = self.hdul['SPMAP'].header
            #self.vigmap_header = self.hdul['VIGMAP'].header
            
            #Now data. 
            self.data = {} 
            self.data['CTSMAP'] = self.hdul['PRIMARY'].data 
            #self.data['CTSMAP'] = self.hdul['CTSMAP'].data 
            #self.data['XBMAP'] = self.hdul['XBMAP'].data 
            #self.data['PSMAP'] = self.hdul['PSMAP'].data 
            #self.data['PBMAP'] = self.hdul['PBMAP'].data 
            #self.data['CLMAP'] = self.hdul['CLMAP'].data 
            #self.data['SPMAP'] = self.hdul['SPMAP'].data 
            #self.data['VIGMAP'] = self.hdul['VIGMAP'].data 

        #Further specific extractions of data. 
        self.get_orbit_info() 
        self.get_camera_info() 

    #PLOTTING FUNCTIONS. 
    ###################
        
    def plot_cts_extension(self, ext = 'CTSMAP', cmap='lundi', vmin=0, vmax=10, save=False):
        '''This will plot one of the extensions for you.''' 
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   
        
        #Use a scale from 0-1 for the vignetting map. Overwrite. 
        #if ext.upper() == 'VIGMAP': 
        #    vmin = 0
        #    vmax = 1 
               
        #Create the figure. 
        fig = plt.figure(figsize=(6,6))
        fig.subplots_adjust(top=0.8, left=0.20, right=0.85)
        ax = fig.add_subplot(111)
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax, self.data['CTSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title='Counts/pixel')
        ax.set_title('CTSMAP'+'\n\n')
        
        #Add figure title with key meta information. 
        #time_title = f'{self.date_obs}'
        #exp_title = f'Expos: {self.expos}s'
        #pos_title = 'SMILE: ({:.2f},{:.2f},{:.2f})'.format(*self.pos) 
        #aim_title = 'AIM: ({:.2f},{:.2f},{:.2f})'.format(*self.aim) 
        #metatitle = time_title+'\n'+exp_title+'\n'+pos_title+'\n'+aim_title 
        
        #fig.text(0.15, 0.95, metatitle, ha='left', va='top') 
        
        #Set filename as the title. 
        fig.text(0.5, 0.95, self.filename, ha='center', fontsize=10) 
               
        if save: 
            print ('Saving: ', self.fitspath+self.filename+'.png')
            fig.savefig(self.fitspath+self.filename+'.png') 
            
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
        self.dtime = dt.datetime.strptime(self.date_obs, '%Y-%m-%dT%H:%M:%S.%f') 
        
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
        self.xdeg_sep = self.primary_header['CDELT1']
        self.ydeg_sep = self.primary_header['CDELT2']
        
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

