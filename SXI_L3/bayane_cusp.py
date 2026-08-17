#This is Bayane's cusp ID code. Taken from her dynamic time integration method: 
#https://github.com/BayaneMdW/SMILE-SXI-DTI-Method/blob/main/SXI_Dynamic_Time_Integration_method.py

import numpy as np 

#These functions are copied directly from her dynamic time integration method. 

def average_along_phi_axis(image: np.ndarray) -> np.ndarray:
    """Averages the image along the phi (azimuthal) axis (axis 0)."""
    return np.mean(image, axis=0)


def average_along_theta_axis(image: np.ndarray) -> np.ndarray:
    """Averages the image along the theta (elevation) axis (axis 1)."""
    return np.mean(image, axis=1)
    
def test2_is_cusp(image: np.ndarray) -> bool:
    """
    Test 2: Determines if the visible structure is a magnetospheric cusp.
    """
    itheta = average_along_theta_axis(image)
    theta_pic = itheta.max() - itheta.min()
    
    iphi = average_along_phi_axis(image)
    phi_pic = iphi.max() - iphi.min()

    return theta_pic >= phi_pic    
