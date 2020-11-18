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
        all_files = directory.rglob('*.mat')
        data_list = [i for i in all_files]
        
        # Fetch all file names in database
        
        for item in data_list:
            directory_names = directory.parts
            mat = scipy.io.loadmat(item)
            data = np.squeeze(np.asarray(mat['val']))
            new_row = {'condition':directory_names[8], 'name':os.path.splitext(item.name)[0], 'ecg':data}
            db = db.append(new_row, ignore_index=True)
            
        return db
    
        # Upload all ECG data from directory to database
        
    def start(self, filepath, database):
        
        dirlist = os.listdir(filepath)
        
        for condition in dirlist:  
            database = FileWizard.load_data(filepath + condition, database)
            
        return database