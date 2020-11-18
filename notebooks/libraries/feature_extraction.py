# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:33:50 2020

@author: Oliver
"""
from pywt import wavedec
from scipy.signal import savgol_filter, find_peaks
import numpy as np

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
        ecg_waves = database['ecg'].tolist()
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
    
    def midpoint(self, x1, y1, x2, y2):
        x = (x1 + x2) / 2
        y = (y1 + y2) / 2
        return x, y
    
    def find_midpoints(self, peaks, position):
        midpoints_x = np.ones(len(peaks) - 1, dtype=float)
        midpoints_y = np.ones(len(peaks) - 1, dtype=float)
        
        pos = position['peak_heights']
        for i, peak in enumerate(peaks):
            if i > 10:
                break
                
            x1 = peaks[i]
            x2 = peaks[i + 1]
            y1 = pos[i]
            y2 = pos[i + 1]
            x, y = MidPointExtractor.midpoint(x1, y1, x2, y2)
            midpoints_x[i] = x
            midpoints_y[i] = y
            
            
        return np.hstack((midpoints_x, midpoints_y))
        
    def fit(self, database):
        midpoints = []
        ecg_peaks = database['peaks'].tolist()
        ecg_positions = database['peak position'].tolist()
        
        for peak, pos in zip(ecg_peaks, ecg_positions):
            midpoint = MidPointExtractor.find_midpoints(peak, pos)
            midpoints.append(midpoint)
        
        database['midpoints'] = midpoints
        
        return database
        
