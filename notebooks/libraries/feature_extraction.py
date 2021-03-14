# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:33:50 2020

@author: Oliver
"""
from pywt import wavedec, waverec
from scipy.signal import savgol_filter, find_peaks
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



class FrequencyExtractor():
    
    # Extract and separate all large and small frequencies from the ECGs
    
    def __init__(self):
        pass
    
    def fit(self, database):
        
        def slope(x, w, c):
            
            # Method for extracting P and T waves
            
            # x is the ECG input
            # w is the window size 
            # c is a threshold for frequency removal
            # (how much of the ECG should be selected when calcing the slope)
            
            freqs = np.zeros_like(x)
            
            # Create zero array to fill with frequencies that are not too large
            
            for i in range(0, x.size):
                
                # This process iterates through the whole ECG
                
                try:
                    freqs[i:i+w] = x[i:i+w]
                    
                    # Select section of ECG according to the size of the window
                    
                    y1 = x[i]
                    y2 = x[i + w]
                    x1 = 0
                    x2 = 0 + w
                    
                    # Take coordinates of points at the start/end of window
    
                    slope = (y2 - y1)/(x2 - x1)
                    
                    # Calculate steepness between the two points 
                    
                    if slope > c or slope < -c:
                        freqs[i:i+w] = np.nan
                        
                    # Remove ECG points in areas where steepness is
                    # bigger than c or less than -c
                        
                    i = i + w
                    
                    # Advance to next w points in ECG
                    
                except IndexError:
                    break
            
            return freqs[~np.isnan(freqs)]
    
        ecg_waves = database['ecg'].tolist()
        
        small = []
        large = []
        
        # Small = P and T waves
        # Large = QRS complexes
        
        for wave in ecg_waves:
            coeffs = wavedec(wave, 'db4', level=6)
            
            # Perform the wave decomposition on the ECG
            
            coeffs[-0] = np.zeros_like(coeffs[-0])
            coeffs[-1] = np.zeros_like(coeffs[-1])
            coeffs[-2] = np.zeros_like(coeffs[-2])
            coeffs[-6] = np.zeros_like(coeffs[-6])
            
            # Coefficients are accessed in reverse here because the coefficients
            # from the python implementation are returned
            # in descending order (D3, D2, D1 instead of D1, D2, D3)
            
            large.append(waverec(coeffs, 'db4'))
            
            # Reconstruct the wave based on coeffs that represent QRS complex
            
            smallfreqs = slope(wave, 15, 3)
            
            # Perform slope operation on wave to remove QRS complex 
            # and keep P + T wave
            
            smallfreqs = savgol_filter(smallfreqs, 43, 7)
            small.append(smallfreqs)
            
        database['small frequencies'] = small
        database['large frequencies'] = large
        
        return database

class PeakExtractor():
    
    def __init__(self, c, p):
        self.c = c
        self.p = p
    
    def fit(self, database, column_name_in, column_name_out1, column_name_out2):
        ecg_waves = database[column_name_in].tolist()
        
        # Find all the peaks for QRS complex
        
        peaks_list = []
        
        # X position of peaks
        
        peak_position = []
        
        # Y position of peaks
        
        for i, wave in enumerate(ecg_waves):
            
            peaks, _ = find_peaks(wave, height=self.c, prominence=self.p, distance=100)
            
            
            # Find all peaks and their positions on one ECG
            #peaks = peaks[peaks > 100]
            #peaks = peaks[peaks < 3500]
            
            # Remove all peaks that do not lie on a complete interval
            position = np.zeros(len(peaks))
            
            for i, peak in enumerate(peaks):
                position[i] = wave[peak]
                
            
            peaks_list.append(peaks)
            peak_position.append(position)
            
            
        database[column_name_out1] = peaks_list
        database[column_name_out2] = peak_position
        
        return database
    
class MidPointExtractor():
    
    def __init__(self):
        pass
    
    def midpoint(x1, y1, x2, y2):
        
        # Calculation of the midpoint
        
        x = (x1 + x2) / 2
        y = (y1 + y2) / 2
        return x, y
    
    def find_midpoints(peaks, position):
        
        # Create empty arrays size of input array minus one
        
        midpoints_x = np.ones(len(peaks) - 1, dtype=float)
        midpoints_y = np.ones(len(peaks) - 1, dtype=float)
        
        # Empty array containing positions of all midpoints on the ECG
        
        pos = position
        # Y position of where midpoint lies
        
        for i in range(len(midpoints_x)):
            
            x1 = peaks[i]
            # Find first peak x position along ECG
            
            x2 = peaks[i + 1]
            # Find next peak x position along ECG
            
            y1 = pos[i]
            # Find first peak y position along ECG
            y2 = pos[i + 1]
            # Find next peak y position along ECG
            x, y = MidPointExtractor.midpoint(x1, y1, x2, y2)
            # Calculate midpoint of the peaks
            midpoints_x[i] = x
            midpoints_y[i] = y       
            # Append midpoint value to array
        return np.vstack([midpoints_x, midpoints_y])
        
    def fit(self, database, in_column_1, in_column_2, out_column):
        midpoints = []
        
        # Empty list to fill midpoints for every ECG up
        
        ecg_peaks = database[in_column_1].tolist()
        ecg_positions = database[in_column_2].tolist()
        
        # Take existing peaks and their positions found from ECG database 
        # to calculate midpoints
        
        i = 0
        
        for peak, pos in zip(ecg_peaks, ecg_positions):
            dim = len(peak) - 1
            if dim > 1:
                midpoint = MidPointExtractor.find_midpoints(peak, pos)
                midpoints.append(midpoint)
            else:
                midpoints.append(np.nan)
            i = i + 1
        
        database[out_column] = midpoints
        database = database.dropna()
        
        return database
    
class WaveletSeparator():
    
    def __init__(self):
        pass
    
    def fit(self, database, new_database):
        for i in range(0, len(database)):
            subject = database.iloc[i]
            
            # Locate ECG record from database
            
            midpoint = subject['midpoints']
            ecg = subject['large frequencies']
            condition = subject['condition']
            name = subject['name']
            
            # Locate all these values from the record
            
            x = midpoint[0,:]
            
            # Only X is necessary because we are splitting along X axis 
            for j in range(0, len(x)):
                if j == 0:
                    wavelet = ecg[0:int(x[j])]
                else:
                    wavelet = ecg[int(x[j - 1]):int(x[j])]
                
                # Append new row for new database containing this information
                
                new_row = {'wavelet':wavelet, 'condition':condition, 'partof':name}
                new_database = new_database.append(new_row, ignore_index=True)
            
        return new_database
    
    