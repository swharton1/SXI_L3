#This file will read in the mapL3 file. 

import os 
import numpy as np 
from astropy.io import fits 
import matplotlib.pyplot as plt 
import datetime as dt

from . import paths 

from SXI_Core import make_image_axes
from SXI_Core import read_cmap 
from dipole import dipole

class read_mapL3():
    '''This class reads in the mapL3 files that Steve has created.'''
    
    def __init__(self, folder='sim_0600/', filename='mapL3_0600.fits', calc_dipole=True):
        '''This reads in the file and constructs the filename. 
        
        Parameters
        ----------
        folder - name of folder. Currently in format 'sim_HHMM' 
        filename - Currently in format 'mapL3_HHMM' 
        calc_dipole - Boolean to calculate the dipole tilt angle in radians for CMEM. Def=True. 
        '''
        
        #Get path to the data. 
        self.datapath = paths.get_data_path()+folder 
        self.fullname = os.path.join(self.datapath, filename) 
        
        #Check the file exists. 
        assert os.path.isfile(self.fullname), f'{self.fullname} does not exist' 
        
        #Now you know it exists, open it. 
        print (f'Read {self.fullname}...')
        with fits.open(self.fullname) as hdul: 
            self.hdul = hdul 
            
            #Read out the entire file. Headers then data. 
            self.primary_header = self.hdul['PRIMARY'].header 
            self.ctsmap_header = self.hdul['CTSMAP'].header
            self.xbmap_header = self.hdul['XBMAP'].header
            self.psmap_header = self.hdul['PSMAP'].header
            self.pbmap_header = self.hdul['PBMAP'].header
            self.clmap_header = self.hdul['CLMAP'].header
            self.spmap_header = self.hdul['SPMAP'].header
            self.vigmap_header = self.hdul['VIGMAP'].header
            
            #Now data. 
            self.data = {} 
            self.data['CTSMAP'] = self.hdul['CTSMAP'].data 
            self.data['XBMAP'] = self.hdul['XBMAP'].data 
            self.data['PSMAP'] = self.hdul['PSMAP'].data 
            self.data['PBMAP'] = self.hdul['PBMAP'].data 
            self.data['CLMAP'] = self.hdul['CLMAP'].data 
            self.data['SPMAP'] = self.hdul['SPMAP'].data 
            self.data['VIGMAP'] = self.hdul['VIGMAP'].data 
            

        #Further specific extractions of data. 
        self.get_orbit_info() 
        self.get_camera_info() 
        
        if calc_dipole:
            self.get_dipole() 
        

        
            
    #FUNCTIONS TO EXTRACT KEY HEADER INFO FROM THE FILE, INCLUDING DIPOLE ANGLE.
    ########################################################################
           
    def get_orbit_info(self):
        '''This will extract the spacecraft position, aim point and time.''' 
        
        #Smile location 
        self.smile_loc = np.array([self.primary_header['POS_X'], self.primary_header['POS_Y'], self.primary_header['POS_Z']]) 
        
        #SXI Aim point 
        self.target_loc = np.array([self.primary_header['AIM_X'], self.primary_header['AIM_Y'], self.primary_header['AIM_Z']])    
        
        #Time  
        self.date_obs = self.ctsmap_header['DATE-OBS'] 
        self.date_end = self.ctsmap_header['DATE-END']
        
        #Get datetime object for the start. 
        self.dtime = dt.datetime.strptime(self.date_obs, '%Y-%m-%dT%H:%M:%S.%f') 
        
        #Get Energy Bands. 
        self.emin = self.ctsmap_header['EMIN']
        self.emax = self.ctsmap_header['EMAX']
        
        #Get Pointing in Sky coords. 
        self.ra = self.ctsmap_header['RA']
        self.dec = self.ctsmap_header['DEC'] 

    def get_camera_info(self):
        '''This will get information about the camera and resolution, including plotting arrays.''' 
            
        #Get information about the camera out. 
        #Number of pixels. 
        self.m_pixels = self.ctsmap_header['NAXIS1']
        self.n_pixels = self.ctsmap_header['NAXIS2'] 
        
        #Pixel widths. 
        self.xdeg_sep = self.ctsmap_header['CDELT1']
        self.ydeg_sep = self.ctsmap_header['CDELT2']
        
        #Lower bounds. 
        self.xdeg_min = self.ctsmap_header['CRVAL1']
        self.ydeg_min = self.ctsmap_header['CRVAL2']
        
        #Units. 
        self.x_unit = self.ctsmap_header['CTYPE1']
        self.y_unit = self.ctsmap_header['CTYPE2'] 
        
        #Calculate FOV. 
        self.phi_fov = -2*self.xdeg_min 
        self.theta_fov = -2*self.ydeg_min 
        
        #Get 1D pixel arrays for plotting. The edges of the pixels.  
        self.xarray = np.linspace(self.xdeg_min, -self.xdeg_min, self.m_pixels+1)
        self.yarray = np.linspace(self.ydeg_min, -self.ydeg_min, self.n_pixels+1)
        
        #Make 2D arrays for x and y. 
        self.X, self.Y = np.meshgrid(self.xarray, self.yarray)
        
        #Exposure 
        self.expos = self.ctsmap_header['EXPOS']   

    def get_dipole(self):
        'This gets the dipole angle in degrees.''' 
        
        #If calculating it from the time. 
        self.dipole = np.deg2rad(dipole.Dipole(self.dtime.year).tilt(self.dtime)) 


    
    #PLOTTING FUNCTIONS. 
    ###################
        
    def plot_raw_extension(self, ext = 'CTSMAP', cmap='lundi', vmin=0, vmax=10):
        '''This will plot one of the extensions for you.''' 
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   
            
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
        pos_title = 'SMILE: ({:.2f},{:.2f},{:.2f})'.format(*self.smile_loc) 
        aim_title = 'AIM: ({:.2f},{:.2f},{:.2f})'.format(*self.target_loc) 
        metatitle = time_title+'\n'+exp_title+'\n'+pos_title+'\n'+aim_title 
        
        fig.text(0.15, 0.95, metatitle, ha='left', va='top') 
