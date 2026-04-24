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
