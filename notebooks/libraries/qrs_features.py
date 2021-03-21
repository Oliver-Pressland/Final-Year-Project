# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 12:58:10 2021

@author: Oliver
"""
from pywt import wavedec, waverec
from scipy.signal import savgol_filter, find_peaks
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class QRSfeatures():
    
    def __init__(self, c):
        self.c = c
    
    def fit(self, database):
        
        column_name_in = database['wavelet']
        
        ## DIFFERENT FEATURE EXTRACTION METHODS ##
        
        # 1. Time length of the sample
        
        def wave_timelength(wave, sample_rate):
            try:
                length = len(wave) 
                time_length = length * sample_rate
                
            except (IndexError, ValueError) as e:
                time_length = np.nan
                
            return time_length
        
        # 2. QR height, RS height, QS time length
            
        def QRS_extraction(wave, sample_rate, gain):
            try:
                high, pos1 = find_peaks(wave, height=self.c)
                low, pos2 = find_peaks(-wave, height=20)
                
                pos1 = wave[high]
                pos2 = wave[low]

                lowmin1_y = min(pos2)
                lowmin1_x = low[pos2.argmin(axis=0)]
                pos2_new = np.delete(pos2, pos2.argmin())

                lowmin2_y = min(pos2_new)
                lowmin2_x = low[np.where(pos2==lowmin2_y)]

                highmax_y = max(pos1)
                highmax_x = high[pos1.argmax(axis=0)]

                q = min(lowmin2_x, lowmin1_x)
                r = highmax_x
                s = max(lowmin2_x, lowmin1_x)
                
                # Find QRS by finding the peaks and troughs of sample
                
                qr_height = (wave[r] - wave[q]) / gain
                rs_height = (wave[r] - wave[s]) / gain
                qs_timelength = (s - q) * sample_rate
                
            except (IndexError, ValueError) as e:
                qr_height = np.nan
                rs_height = np.nan
                qs_timelength = np.nan
                
                # This catches all anomalous values in case of an exception
                
            return qr_height, rs_height, qs_timelength
        
        # 3. Time gap between two samples
        
        def timegap_length(this_wave, next_wave, sample_rate):
            try:
                
                total = np.hstack((this_wave, next_wave))
                
                peaks, pos = find_peaks(total, height=self.c, prominence=30)

                r1 = peaks[0]
                r2 = peaks[-1]
                
                timegap = (r2 - r1) * sample_rate
            
            except (IndexError, ValueError) as e:
                timegap = np.nan
                
                # This catches all anomalous values in case of an exception
                
            return timegap
        
        def arraycheck(value):
            
            # Convert array value to scalar if this is the case
            
            if type(value) is np.ndarray:
                new_value = value.item()
            else:
                new_value = value
            return new_value
        
        ## START OF THE METHOD ##
        
        time_lengths = []
        qr_heights = []
        rs_heights = []
        qs_timelengths = []
        timegaps = []
        
        f = 0.002777777777777777777777
        g = 200
        
        for i, this_item in enumerate(column_name_in):
            try:
            
                # Extract time length 
                time_length = wave_timelength(this_item, f)
                time_lengths.append(time_length)
                
                # Extract QRS info
                qr_height, rs_height, qs_timelength = QRS_extraction(this_item, f, g)
                
                qr_height = arraycheck(qr_height)
                rs_height = arraycheck(rs_height)
                qs_timelength = arraycheck(qs_timelength)
                
                qr_heights.append(qr_height)
                rs_heights.append(rs_height)
                qs_timelengths.append(qs_timelength)
            
                # Timegaps extractions
    
                next_item = column_name_in[i + 1]
        
                timegap = timegap_length(this_item, next_item, f)
                timegap = arraycheck(timegap)

                timegaps.append(timegap)
                
            except KeyError as e:
                timegaps.append(np.nan)

         
        database['time length'] = time_lengths
        database['qr height'] = qr_heights
        database['rs height'] = rs_heights
        database['qs time length'] = qs_timelengths
        database['time gap'] = timegaps
        
        return database
        
                