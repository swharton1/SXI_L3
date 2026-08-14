#This will calculate the quality flag. 

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.ticker import MultipleLocator
import datetime as dt 
import matplotlib.dates as dates

def calc_quality_flag(aim, pos, expos, cxfov):
    '''This uses a points system based on spacecraft position and count rate.

    Parameters
    ----------
    aim - 3 element array for the aim point (ax, ay, az)
    pos - 3 element array for the spacecraft position (px, py, pz)
    expos - exposure time in seconds 
    cxfov - array representing the foreground emission. i.e. CXFOV.

    ''' 

    #Points from spacecraft position. 
    #This calculates the perpendicular angle lambda, named after Richard Hampson's favourite Greek letter. 

    tan_lda = (aim[0] - pos[0])/np.sqrt(pos[1]**2 + pos[2]**2) 
    lda = np.rad2deg(np.abs(np.arctan(tan_lda))) 

    if lda <= 25: ps = 0 
    elif (lda > 25) & (lda <= 50): ps = 1 
    else: ps = 2 

    #Now get a score based on image count rate. 
    total_counts = cxfov.sum() 
    cr = total_counts/expos

    if cr >= 150: pc = 0 
    elif (cr < 150) & (cr >= 75): pc = 1
    else: pc = 2 

    #Now get total quality flag. 
    qf = ps + pc 

    return qf 
    
def calc_quality_flag_2(aim, pos, cxfov, bkg):
    '''This quality flag is based on position and the signal to noise ratio. It also uses a binary system to encode the flag.
    
    Parameters
    ----------
    aim - 3 element array for the aim point (ax, ay, az)
    pos - 3 element array for the spacecraft position (px, py, pz)
    cxfov - array representing the foreground emission. i.e. CXFOV.
    bkg - array representing the foreground emission. i.e. BKGMAP.
    
    '''
    
    #Points from spacecraft position. 
    #This calculates the perpendicular angle lambda, named after Richard's favourite Greek letter. 

    tan_lda = (aim[0] - pos[0])/np.sqrt(pos[1]**2 + pos[2]**2) 
    lda = np.rad2deg(np.abs(np.arctan(tan_lda))) 
    
    #Assign position bits. 
    bit1 = 1 if (lda > 25) & (lda <= 50) else 0 
    bit2 = 2 if (lda > 50) else 0 
    
    #Calculate SNR from total foreground and background counts. 
    S = cxfov.sum() 
    B = bkg.sum() 
    
    SNR = S/(S**2 + B**2)**0.5
    
    #Assign SNR bits. 
    bit3 = 4 if (SNR > 0.33) & (SNR <= 0.67) else 0 
    bit4 = 8 if (SNR < 0.33) else 0 
    
    #Quality flag is the sum of all the bits. 
    qf = bit1 + bit2 + bit3 + bit4
    
    print (lda, SNR) 
    
    return qf 
    
def decompose_qf2(qf):
    '''This decomposes a quality flag bag into a binary number. Returns as a string from lowest to highest digit.
    
    Parameters
    ----------
    qf - quality flag total number. 
    
    Returns
    -------
    binary - binary array number. e.g. [1,0,0,1]
    
    ''' 
    
    #Check qf is a valid value for MY SYSTEM. 
    valid_qfs = [0,1,2,4,5,6,8,9,10]
    assert qf in valid_qfs, f'{qf} not a valid value in my system currently. Valid values are {valid_qfs}.' 
    
    binary = ''  # binary result

    while qf > 0:
        binary = str(qf & 1) + binary
        qf >>= 1
    
    #Pad to 4 digits if needed. Needed for smaller numbers.  
    binary = binary.zfill(4) 
    
    #Reorder from low to high. 
    binary = binary[::-1] 
    
    #Convert to a numpy array of integers. 
    binary = np.array(list(binary)).astype('int') 
    
    return binary

def interpret_qf(qf, plot=True):
    '''This interprets the binary number. Prints out what each digit means. 
    Option to make a simple barplot of it 
    
    Parameters
    ----------
    qf - quality flag total number. 
    plot - boolean to plot out the bits and interpret them. 
    
    '''  
    
    #Get binary number. 
    binary = decompose_qf2(qf) 
    
    #Get length of binary number.
    length = len(binary) 
     
    #Go through each digit to work out which flags are on. 
    for i in range(length):
        print (i)
        if (i == 0) & (binary[i] == 1): print ('lda between 25 and 50') 
        if (i == 1) & (binary[i] == 1): print ('lda more than 50') 
        if (i == 2) & (binary[i] == 1): print ('SNR between 0.33 and 0.67') 
        if (i == 3) & (binary[i] == 1): print ('SNR less than 0.33') 

    #Recalcu
    if plot: 
        #plt.close("all")        
        fig = plt.figure()
        ax1 = fig.add_subplot(111) 
        x = np.arange(length)
        ax1.bar(x, binary) 
        ax1.xaxis.set_major_locator(MultipleLocator(1))
        ax1.yaxis.set_major_locator(MultipleLocator(1)) 
        ax1.set_xticks([0,1,2,3])
         
        #Make labels. 
        labels = [r'25$^{\circ}$ < $\lambda$ $\leq$ 50$^{\circ}$', r'$\lambda$ > 50$^{\circ}$', r'0.33 < SNR $\leq$ 0.67', 'SNR < 0.33'] 
         
        ax1.set_xticklabels(labels) 
         
        ax1.set_title(f'QF Decomposition: {binary} - {qf}')
    

          
def interpret_qf_list(times, qf_list, ax1=None, cmap='Greys'):
    '''This will take a list of quality flags and the times they were taken and decompose each one into its bits. It then produces a plot to show the variation in the bits of the quality flags. Returns the axes. 
    
    Parameters
    ----------
    times - list/array of datetime objects (or just numbers).
    qf_list - list/array of total quality flag numbers.
    ax1 - If None, it will make one here. 
    
    Returns
    -------
    ax - Returns the axis with the QF decomposition displayed so it can be used in other plots. Further customisation of the axis can be done in another programme. 
    
    '''
    
    #Get length of qf_list. 
    qflen = len(qf_list) 
    
    #Make array to store all bits in. 
    bits = np.zeros((qflen,4)) 
    
    for q, qf in enumerate(qf_list):
        bits[q] = decompose_qf2(qf) 
    
    
    #Make array of indices for bits. 
    x = [-0.5,0.5,1.5,2.5,3.5]
    
    #Add to the end of the times array. 
    times_extra = np.concatenate((times, [times[-1]+dt.timedelta(seconds=300)]))  
    
    #Make 2D arrays for plotting. 
    X, TIMES = np.meshgrid(x, times_extra) 
    
    #Now create a plot.
    if ax1 is None: 
        fig = plt.figure(figsize=(8,4))
        fig.subplots_adjust(left=0.2)
        ax1 = fig.add_subplot(111)
        
    ax1.pcolormesh(TIMES, X, bits, cmap=cmap, vmin=-0.2, vmax=1)
    ax1.yaxis.set_major_locator(MultipleLocator(1)) 
    ax1.set_yticks([0,1,2,3]) 
    
    #Make labels. 
    labels = [r'25$^{\circ}$ < $\lambda$ $\leq$ 50$^{\circ}$', r'$\lambda$ > 50$^{\circ}$', r'0.33 < SNR $\leq$ 0.67', 'SNR < 0.33']      
    ax1.set_yticklabels(labels)   
    
    t_form = dates.DateFormatter('%H:%M')
    ax1.xaxis.set_major_formatter(t_form)
    
    return ax1 
    
    
def generate_times(expos=300, n=10):
    '''This is a simple function to generate some datetime objects for testing interpret_qf_list().'''
    
    start = dt.datetime.now() 
    diff = dt.timedelta(seconds=expos)
    
    dtime_list = np.array([start+diff*i for i in range(n)])   
    
    return dtime_list                     
