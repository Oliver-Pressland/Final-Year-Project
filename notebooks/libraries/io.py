import os
import numpy as np
import scipy.io
import pandas as pd
from pathlib import Path


class FileWizard():
       
    def __init__(self):
        pass
    
    def load_data(path, db):
        directory = Path(path)
        
        # Access file path in filesystem of given filepath
        
        all_files = directory.rglob('*.mat')
        
        # Fetch the folder with the ECG .mat files
        
        data_list = [i for i in all_files]
        
        # Data_list fetches all file names in database
        
        for item in data_list:
            directory_names = directory.parts
            
            # Splits up the directory into a text array
            # so that names of files, condition type can be appended to data table
            
            mat = scipy.io.loadmat(item)
            data = np.squeeze(np.asarray(mat['val']))
            
            # The raw 'mat' file ECG data is under a column called ['val']
            # so that needs to be extracted as an array
            
            new_row = {'condition':directory_names[-1], 
                       'name':os.path.splitext(item.name)[0], 'ecg':data}
            db = db.append(new_row, ignore_index=True)
            
            # Add ECG to table along with filename, condition type.
            
        return db
    
    def start(self, filepath, database):
        
        dirlist = os.listdir(filepath)
        
        # Fetch all folders with the different conditions
        
        for condition in dirlist:  
            database = FileWizard.load_data(filepath + condition, database)
            
            # Upload the ECG daa for every condition
            
        return database
    
    # Upload all ECG data from directory to database
    
