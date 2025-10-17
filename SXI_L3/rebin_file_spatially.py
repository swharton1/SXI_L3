#This will take the three types of file and produce a single file with the data rebinned spatially, then CXFOV calculated. 

import os 
import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits 

from . import paths 
from . import read_bkg_file
from . import read_tot_file
from . import read_vcy_file  
from . import rebin_funcs 
from .SXI_Core import read_cmap 
from .SXI_Core import make_image_axes

class rebin_file():
    '''This will take all three files and produce a single rebinned file.''' 
    
    def __init__(self, folder='L3_20260317T0240-20260317T0241/', bkg_file='SMILE_SXI_L3_SCIM15-SCI-BKG_20260317T0240-20260317T0241_V01.fits', tot_file='SMILE_SXI_L3_SCIM15-SCI-TOT_20260317T0240-20260317T0241_V01.fits', vcy_file='SMILE_SXI_L3_SCIM15-SCI-VCY_20260317T0240-20260317T0241_V01.fits', xres=1, yres=1):
    
        #Get path to the data. 
        self.folder = folder
        self.datapath = paths.get_data_path()+folder 
        self.fitspath = paths.get_fits_path()
        self.bkg_fullname = os.path.join(self.datapath, bkg_file)
        self.tot_fullname = os.path.join(self.datapath, tot_file)
        self.vcy_fullname = os.path.join(self.datapath, vcy_file)
        self.bkg_file = bkg_file
        self.tot_file = tot_file
        self.vcy_file = vcy_file 
        
        #Check the file exists. 
        assert os.path.isfile(self.bkg_fullname), f'{self.bkg_fullname} does not exist' 
        assert os.path.isfile(self.tot_fullname), f'{self.tot_fullname} does not exist' 
        assert os.path.isfile(self.vcy_fullname), f'{self.vcy_fullname} does not exist' 
        
        self.xres = xres
        self.yres = yres 
        
        #Read in each filetype. 
        self.bkg = read_bkg_file.read_bkg_file(folder=folder, filename=bkg_file)   
        self.tot = read_tot_file.read_tot_file(folder=folder, filename=tot_file)    
        self.vcy = read_vcy_file.read_vcy_file(folder=folder, filename=vcy_file) 
        
        #Rebin all the arrays spatially. 
        self.rebin_arrays() 
        
        #Now calculate CXFOV. 
        #print ('Calculating CXFOV...') 
        #self.rebin_data['CXFOV'] = (self.rebin_data['CTSMAP'] - self.rebin_data['BKGMAP'])/self.rebin_data['VIGMAP'] 
        
        #Calculate rebinning factor. i.e. how many times larger the box is. 
        self.rebin_factor = (self.xres/self.bkg.xdeg_sep)*(self.yres/self.bkg.ydeg_sep)
        
    def rebin_arrays(self):
        '''This will rebin all the arrays spatially and add them to a new dictionary.
        
        Parameters
        ----------
        bkg - file object returned from read_bkg_file. 
        tot - file object returned from read_tot_file. 
        vcy - file object returned from read_vcy_file. 
        
        Returns
        -------
        rebin_data - dictionary containing all the spatially rebinned extensions. 
        
        ''' 
        
        print ('Rebinning data spatially...') 
        
        #Check xres and yres are appropriate choices. 
        assert self.xres >= self.bkg.xdeg_sep, f'xres must be greater than {self.bkg.xdeg_sep}'
        assert self.yres >= self.bkg.ydeg_sep, f'yres must be greater than {self.bkg.ydeg_sep}'
        
        #Get new shape for the arrays. 
        self.m_pixels = self.bkg.phi_fov/self.xres 
        self.n_pixels = self.bkg.theta_fov/self.yres 
        
        #Check this gives an integer shape. 
        assert self.m_pixels.is_integer(), f'xres = {self.xres} did not give an integer number of pixels for m_pixels. Try something else.' 
        assert self.n_pixels.is_integer(), f'yres = {self.yres} did not give an integer number of pixels for n_pixels. Try something else.' 
        
        self.m_pixels = int(self.m_pixels)
        self.n_pixels = int(self.n_pixels) 
        
        self.xdeg_min = self.bkg.xdeg_min 
        self.ydeg_min = self.bkg.ydeg_min 
        self.x_unit = self.bkg.x_unit
        self.y_unit = self.bkg.y_unit
        
        
        self.emin = self.bkg.emin
        self.emax = self.bkg.emax
        self.ra = self.bkg.ra
        self.dec = self.bkg.dec
        self.expos = self.bkg.expos
        
        self.date_obs = self.bkg.date_obs
        self.date_end = self.bkg.date_end
        
        #Now get the new shape. 
        self.new_shape = (self.n_pixels, self.m_pixels)  
        
        #Now start rebinning all the individual data products.  
        rebin_data = {}  
        rebin_data['CTSMAP'] = rebin_funcs.rebin_sum(self.tot.data['CTSMAP'], self.new_shape) 
        rebin_data['BKGMAP'] = rebin_funcs.rebin_sum(self.bkg.data['BKGMAP'], self.new_shape)
        rebin_data['XBMAP'] = rebin_funcs.rebin_sum(self.bkg.data['XBMAP'], self.new_shape) 
        rebin_data['PSMAP'] = rebin_funcs.rebin_sum(self.bkg.data['PSMAP'], self.new_shape) 
        rebin_data['PBMAP'] = rebin_funcs.rebin_sum(self.bkg.data['PBMAP'], self.new_shape) 
        rebin_data['CLMAP'] = rebin_funcs.rebin_sum(self.bkg.data['CLMAP'], self.new_shape) 
        rebin_data['SPMAP'] = rebin_funcs.rebin_sum(self.bkg.data['SPMAP'], self.new_shape) 
        rebin_data['VIGMAP'] = rebin_funcs.rebin_mean(self.vcy.data['VIGMAP'], self.new_shape) 
        
        #Attach the dictionary to the object. 
        self.rebin_data = rebin_data 
        
    def plot_rebinned_data(self, cmap='lundi', vmin=0, vmax=10, save=False):
        '''This will plot the rebinned CTSMAP, BKGMAP and VIGMAP alongside the originals, as well as the CXFOV map. '''
    
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()   
            
        
        #Create the figure. 
        fig = plt.figure(figsize=(8,6))
        fig.subplots_adjust(top=0.80, wspace=0.5, hspace=0.4)
        ax1 = fig.add_subplot(231)
        ax2 = fig.add_subplot(232)
        ax3 = fig.add_subplot(233)
        ax4 = fig.add_subplot(234)
        ax5 = fig.add_subplot(235)
        ax6 = fig.add_subplot(236)
        #ax7 = fig.add_subplot(338)
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax1, self.tot.data['CTSMAP'], self.tot.xdeg_min, self.tot.ydeg_min, self.tot.n_pixels, self.tot.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax/self.rebin_factor, cbar_title='CTSMAP [cts/pix]')
        ax1.set_title('CTSMAP'+'\n')
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax2, self.bkg.data['BKGMAP'], self.bkg.xdeg_min, self.bkg.ydeg_min, self.bkg.n_pixels, self.bkg.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax/self.rebin_factor, cbar_title='BKGMAP [cts/pix]', ylabel=False)
        ax2.set_title('BKGMAP'+'\n')
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax3, self.vcy.data['VIGMAP'], self.vcy.xdeg_min, self.vcy.ydeg_min, self.vcy.n_pixels, self.vcy.m_pixels, cmap=cmap, vmin=0, vmax=1, cbar_title='VIGMAP', ylabel=False)
        ax3.set_title('VIGMAP'+'\n')  
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax4, self.rebin_data['CTSMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title='CTSMAP [cts/pix]')
        #ax4.set_title('CTSMAP'+'\n')  
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax5, self.rebin_data['BKGMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title='BKGMAP [cts/pix]', ylabel=False)
        #ax5.set_title('BKGMAP'+'\n')  
        
        #Make the axis. 
        make_image_axes.make_image_axes(ax6, self.rebin_data['VIGMAP'], self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=0, vmax=1, cbar_title='VIGMAP', ylabel=False)
        #ax6.set_title('VIGMAP'+'\n')  
        

        
        #Add titles from the files. 
        fig.text(0.5, 0.97, f'{self.tot_file}\n{self.bkg_file}\n{self.vcy_file}', ha='center', fontsize=10, va='top')
            
        if save: 
            filename = f'SMILE_SXI_L3_SCIM15-SCI-ALL_{self.folder[:-1]}_{self.xres}x{self.yres}_rebin.png'
            print ('Saving: ', self.fitspath+filename)
            fig.savefig(self.fitspath+filename)  
