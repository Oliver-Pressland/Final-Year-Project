# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 13:10:27 2020

@author: Oliver
"""
import os
import numpy as np
import scipy.io
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

database = pd.DataFrame(columns=['condition', 'name', 'ecg'])

from libraries.io import FileWizard

path1 = 'C:/Users/Oliver/Documents/FYP/code/database/MLII/'

fw = FileWizard()
database = fw.start(path1, database)
    
from libraries.noise_removal import BaselineNoiseRemover

# DC Notch filter to remove baseline noise from all signals

bnr = BaselineNoiseRemover(c = -0.99)

ecg_waves = database['ecg'].tolist()
ecg_filt = []

for wave in ecg_waves:
    filt = bnr.fit(wave)
    ecg_filt.append(filt)

database['ecg'] = pd.Series(ecg_filt)

from libraries.feature_extraction import LargeFrequencyExtractor

lfe = LargeFrequencyExtractor()
database = lfe.fit(database)
# Multilevel discrete decomposition to extract large frequencies from time series

from libraries.feature_extraction import PeakExtractor

LISTS2 = ['3 AFL', '4 AFIB', '5 SVTA', '6 WPW', 
         '7 PVC', '8 Bigeminy', '9 Trigeminy', '10 VT', '11 IVR', 
         '12 VFL', '13 Fusion', '14 LBBBB', '15 RBBBB', '16 SDHB', '17 PR']

for item in LISTS2:
    database = database[database['condition'] != item]

thresh = 20
pe = PeakExtractor(c=thresh)
database = pe.fit(database)

examples = database[database['condition'] == '1 NSR']
example1 = examples.iloc[1]
peaks1 = example1['peaks']
position1 = example1['peak position']
ecg1 = example1['coefficient 4']

from libraries.feature_extraction import MidPointExtractor

mpe = MidPointExtractor()
database = mpe.fit(database)

ecg = database.iloc[0]
print(ecg['midpoints'])

qrs_db = pd.DataFrame(columns=['condition', 'name', 'ecg'])

for i in range(0, len(database)):
    subject = database.iloc[i]
    midpoint = subject['midpoints']
    ecg = subject['coefficient 4']
    condition = subject['condition']
    name = subject['name']
    x = midpoint[0,:]
    for j in range(0, len(x)):
        if j == 0:
            wavelet = ecg[0:int(x[j])]
        else:
            wavelet = ecg[int(x[j - 1]):int(x[j])]
        print(j, wavelet)
        new_row = {'wavelet':wavelet, 'condition':condition, 'partof':name}
        qrs_db = qrs_db.append(new_row, ignore_index=True)

examples = qrs_db[qrs_db['condition'] == '1 NSR']
for i in range(0, 50):
    e1 = examples.loc[i]
    plt.plot(e1['wavelet'])