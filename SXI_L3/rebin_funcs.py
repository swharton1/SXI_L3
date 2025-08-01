#These functions are needed for spatial rebinning. 

def rebin_mean(arr, new_shape):
    '''This is a function I've got from
    https://scipython.com/blog/binning-a-2d-array-in-numpy/
    Use for all count maps. 
    
    Parameters
    ----------
    arr - The 2D array you wish to rebin. 
    new_shape - The new shape you wish to rebin it too. 

    Returns
    -------
    reshaped array. 
    '''

    shape = (new_shape[0], arr.shape[0] // new_shape[0],
     new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def rebin_sum(arr, new_shape):
    '''This is a function I've got from
    https://scipython.com/blog/binning-a-2d-array-in-numpy/
    This is adapted to find the sum instead of the mean. 
    Use for vignetting. 

    Parameters
    ----------
    arr - The 2D array you wish to rebin. 
    new_shape - The new shape you wish to rebin it too. 

    Returns
    -------
    reshaped array. 
    '''

    shape = (new_shape[0], arr.shape[0] // new_shape[0],
     new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).sum(-1).sum(1)   

