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
            
        def QRS_extraction(wave, sample_rate):
            try:
                high, pos1 = find_peaks(wave, height=self.c)
                low, pos2 = find_peaks(-wave, height=20)
                
                #high = high[high > 100]
                #low = low[low > 100]
            
                q = low[0]
                r = high[0]
                s = low[-1]
                
                # Find QRS by finding the peaks and troughs of sample
                
                qr_height = wave[r] - wave[q]
                rs_height = wave[r] - wave[s]
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
                
                peaks, pos = find_peaks(total, height=self.c)
            
                r1 = peaks[0]
                r2 = peaks[1]
                
                timegap = (r2 - r1) * sample_rate
            
            except (IndexError, ValueError) as e:
                timegap = np.nan
                
                # This catches all anomalous values in case of an exception
                
            return timegap
        
        ## START OF THE METHOD ##
        
        time_lengths = []
        qr_heights = []
        rs_heights = []
        qs_timelengths = []
        timegaps = []
        
        f = 0.002777777777
        
        for i, this_item in enumerate(column_name_in):
            try:
            
                # Extract time length 
                time_length = wave_timelength(this_item, f)
                time_lengths.append(time_length)
                
                # Extract QRS info
                qr_height, rs_height, qs_timelength = QRS_extraction(this_item, f)
                qr_heights.append(qr_height)
                rs_heights.append(rs_height)
                qs_timelengths.append(qs_timelength)
            
                # Timegaps extractions
    
                next_item = column_name_in[i + 1]
                timegap = timegap_length(this_item, next_item, f)
                timegaps.append(timegap)
                
            except KeyError as e:
                timegaps.append(np.nan)

         
        database['time length'] = time_lengths
        database['qr height'] = qr_heights
        database['rs height'] = rs_heights
        database['qs time length'] = qs_timelengths
        database['time gap'] = timegaps
        
        return database
        
                