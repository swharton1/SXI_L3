#This will make a proper python configuration file to store paths in. 

#See: https://configu.com/blog/working-with-python-configuration-files-tutorial-best-practices/#:~:text=Configuration%20files%20in%20Python%20are,can%20be%20read%20and%20modified.

from configparser import ConfigParser

#Make the config object. 
config_object = ConfigParser()

#Add information to the object. 
config_object["PATHINFO"] = {
    "data_path": "/data/smile/shared/sims/GAC_V2/",
    "fits_path": "/home/s/sw682/Code/SXI_L3/SXI_L3/binned_examples/",
    "plot_path": "/home/s/sw682/Code/plots/SXI_L3_plots/" 
}

#Write the configuration file. 
with open('config.ini', 'w') as conf:
    config_object.write(conf) 
    
