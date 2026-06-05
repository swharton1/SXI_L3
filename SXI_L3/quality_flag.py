#This will calculate the quality flag. 

import numpy as np 

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
    #This calculates the perpendicular angle lambda, named after Richard's favourite Greek letter. 

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
    bit3 = 4 if (SNR > 1) & (SNR <= 2) else 0 
    bit4 = 8 if (SNR < 1) else 0 
    
    #Quality flag is the sum of all the bits. 
    qf = bit1 + bit2 + bit3 + bit4
    
    print (lda, SNR) 
    
    return qf 
    
    
