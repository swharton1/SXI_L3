#This file will read the config file. 

from configparser import ConfigParser
import os 
from pathlib import Path 
import sys 

def read_config(path_type:str ='plot_path', section:str ='PATHINFO'):
    '''This will read the config file to extract a particular path.
    
    Parameters
    ----------
    path_type - which path to get. def = 'plot_path' 
    section - which section of the config file. def = 'PATHINFO'
    
    Returns
    -------
    desired_path - path as a string object to return. 
    
    ''' 
    
    #Create a config object. 
    config_object = ConfigParser() 

    #Get the path to the config file. 
    dir_to_file = Path(__file__)
    parent = dir_to_file.parents[1]
    
    #Get the config file name. 
    name = parent/'config.ini'
    
    #Check you can find the file. 
    assert os.path.exists(name), f'{name} does not exist.' 

    #Read out the paths into environment variables. 
    config_object.read(name) 
    
    try:
        pathinfo = config_object[section.upper()] 
    except KeyError as k: 
        print (f'{k} not a valid section in config.ini')
        sys.exit(1)
        
    try: 
        desired_path = pathinfo[path_type.lower()] 
    except KeyError as k:
        print (f'{k} not a valid key in {section} in config.ini') 
        sys.exit(1)

    
    return desired_path 

