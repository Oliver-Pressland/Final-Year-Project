# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:33:50 2020

@author: Oliver
"""
from pywt import wavedec
from scipy.signal import savgol_filter, find_peaks
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class LargeFrequencyExtractor():
    
    def __init__(self):
        pass
    
    def fit(self, database):
        ecg_waves = database['ecg'].tolist()
        
        coeff1 = []
        coeff2 = []
        coeff3 = []
        coeff4 = []
        
        for wave in ecg_waves:
            coeffs = wavedec(wave, 'db1', level=4)
            coeff1.append(coeffs[0])
            coeff2.append(coeffs[1])
            coeff3.append(coeffs[2])
            coeff4.append(coeffs[3])
        
        database['coefficient 1'] = coeff1
        database['coefficient 2'] = coeff2
        database['coefficient 3'] = coeff3
        database['coefficient 4'] = coeff4
        # Multilevel discrete decomposition of ECG waves for compression and noise reduction.
        
        return database

class PeakExtractor():
    
    def __init__(self, c):
        self.c = c
    
    def fit(self, database):
        ecg_waves = database['coefficient 4'].tolist()
        peaks_list = []
        peak_position = []
        
        for wave in ecg_waves:
            peaks, position = find_peaks(wave, height=self.c)
            peaks_list.append(peaks)
            peak_position.append(position)
            
        database['peaks'] = peaks_list
        database['peak position'] = peak_position
        
        return database
    
class MidPointExtractor():
    
    def __init__(self):
        pass
    
    def midpoint(x1, y1, x2, y2):
        x = (x1 + x2) / 2
        y = (y1 + y2) / 2
        return x, y
    
    def find_midpoints(peaks, position):
        
        midpoints_x = np.ones(len(peaks) - 1, dtype=float)
        midpoints_y = np.ones(len(peaks) - 1, dtype=float)
        
        pos = position['peak_heights']
        for i in range(len(midpoints_x)):
            x1 = peaks[i]
            x2 = peaks[i + 1]
            y1 = pos[i]
            y2 = pos[i + 1]
            x, y = MidPointExtractor.midpoint(x1, y1, x2, y2)
            midpoints_x[i] = x
            midpoints_y[i] = y       
        return np.vstack([midpoints_x, midpoints_y])
        
    def fit(self, database):
        midpoints = []
        ecg_peaks = database['peaks'].tolist()
        ecg_positions = database['peak position'].tolist()
        i = 0
        for peak, pos in zip(ecg_peaks, ecg_positions):
            dim = len(peak) - 1
            if dim > 1:
                midpoint = MidPointExtractor.find_midpoints(peak, pos)
                midpoints.append(midpoint)
                #print(database.iloc[i])
            else:
                midpoints.append(np.nan)
            i = i + 1
        
        database['midpoints'] = midpoints
        database = database.dropna()
        
        return database
        
class WaveletSeparator():
    
    def __init__(self):
        pass
    
    def fit(self, database, new_database):
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
                new_row = {'wavelet':wavelet, 'condition':condition, 'partof':name}
                new_database = new_database.append(new_row, ignore_index=True)
            
        return new_database
    
class QRSHeightExtractor():
    
    def __init__(self, c):
        self.c = c
    
    def fit(self, database):
        heights = []
        
        for wave in database['wavelet']:
            high, pos1 = find_peaks(wave, height=self.c)
            low, pos2 = find_peaks(-wave, height=self.c)
            
            height1 = (high - low)
            
            pos1 = pos1['peak_heights']
            pos1 = pos1[0]
            pos2 = pos2['peak_heights']
            pos2 = pos2[0]
            
            pos2 = -pos2
            pos2 = abs(pos2)
            
            form = (pos2 + pos1)
            height = form
            heights.append(height)
            
        return heights

class IntervalLengthExtractor():
    
    def __init__(self):
        pass
        
    def fit(self, database):
        lengths = []
        
        for wave in database['wavelet']:
            length = len(wave)
            lengths.append(length)
            
        return lengths
    
class SQLengthExtractor():
    
    def __init__(self, c):
        self.c = c
    
    def fit(self, database):
        lengths = []
        wl = database['wavelet']
        
        for i, wave in enumerate(wl):     
            if i == 0:
                length = np.nan
                lengths.append(length)
            
            else:
                pastwave = wl[i - 1]
                highpast, pos1 = find_peaks(pastwave, height=self.c)
                low, pos2 = find_peaks(-wave, height=self.c)
                
                pos1 = pos1['peak_heights']
                pos1 = pos1[0]
                pos2 = pos2['peak_heights']
                pos2 = pos2[0]
                
                end = len(pastwave)
                post_qrs = pastwave[int(highpast[0]):int(end)]
                post_qrs = len(post_qrs)
                pre_qrs = len(wave[int(0):int(low[0])])
                length = pre_qrs + post_qrs
                lengths.append(length)
                
        return lengths