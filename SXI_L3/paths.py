#This contains important paths. 
import os

#def get_plot_path(path='/home/s/sw682/Code/plots/SXI_L3_plots/'):
#    assert os.path.exists(path), f'{path} does not exist!' 
#    return path 
    
def get_data_path(path='/data/smile/shared/sims/GAC_V1/'):
    assert os.path.exists(path), f'{path} does not exist!'
    return path  
    
def get_fits_path(path=os.path.dirname(__file__)+'/binned_examples/'):
    assert os.path.exists(path), f'{path} does not exist!'
    return path  
