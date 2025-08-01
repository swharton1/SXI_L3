#This function will fill in the image axes, adding histograms, scaling the images axes and adding colourbars in the standard way. 

import numpy as np 
import matplotlib.pyplot as plt

def make_image_axes(ax, data, xdeg_min, ydeg_min, n_pixels, m_pixels, cmap='hot', vmin=0, vmax=50, add_cbar=True, cbar_title='Counts/pixel'):
    '''This will take the data and axis and fill it in for you.
    
    Parameters
    ----------
    ax - axis object. 
    data - 2D array of data. Can be intensity, counts, etc. 
    xdeg_min - Minimum x value (phi pixel value) 
    ydeg_min - Minimum y value (theta pixel value) 
    n_pixels - Number of pixels in y direction 
    m_pixels - Number of pixels in x direction 
    cmap - colourmap object. def = 'hot'
    vmin - minimum value on colour scale. def = 0
    vmax - maximum value on colour scale. def = 50 
    add_cbar - boolean to add on a colour bar. 
    cbar_title - colourbar title. def = 'Counts/pixel'
    '''
    
    #Get 1D pixel arrays for plotting. The edges of the pixels.  
    xarray = np.linspace(xdeg_min, -xdeg_min, m_pixels+1)
    yarray = np.linspace(ydeg_min, -ydeg_min, n_pixels+1)
        
    #Make 2D arrays for x and y. 
    X, Y = np.meshgrid(xarray, yarray)
    
    #Make the image 
    mesh = ax.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax) 
    
    if add_cbar:
        #Add colourbars. 
        cbar = plt.colorbar(mesh, ax=ax, shrink=0.5, pad=0.15)
        cbar.set_label(cbar_title)
        
    #Histograms. 
    hist1 = data.sum(axis=0)
    hist2 = data.sum(axis=1) 
    
    phi_fov = -2*xdeg_min 
    theta_fov = -2*ydeg_min  
     
    xdeg_sep = phi_fov/m_pixels
    ydeg_sep = theta_fov/n_pixels
                  
    #Scale bar plot data. 
    n1 = 10*hist1.max()/theta_fov 
    n2 = 10*hist2.max()/phi_fov
        
    ax.bar(xarray[0:-1], hist1/n1, width=xdeg_sep, align='edge', bottom=theta_fov/2, clip_on=False, edgecolor='lightgrey', facecolor='lightgrey')
        
    ax.barh(yarray[0:-1], hist2/n2, height=ydeg_sep, align='edge', left=phi_fov/2, clip_on=False, edgecolor='lightgrey', facecolor='lightgrey')
    
    #Add labels. 
    ax.set_xlabel('deg')
    ax.set_ylabel('deg') 
    ax.set_xlim(-phi_fov/2, phi_fov/2) 
    ax.set_ylim(-theta_fov/2, theta_fov/2)
    ax.set_aspect('equal')
