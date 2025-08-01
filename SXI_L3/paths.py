#This contains important paths. 
import os

def get_plot_path(path='/home/s/sw682/Code/plots/SXI_L3_plots3/'):
    assert os.path.exists(path), f'{path} does not exist!' 
    return path 
    
def get_data_path(path='/data/smile/shared/sims/GAC_OUTPUT/'):
    assert os.path.exists(path), f'{path} does not exist!'
    return path  
