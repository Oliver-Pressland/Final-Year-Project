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
    
    def __init__(self):
        pass
    
    def fit(self, database):
        ecg_waves = database['ecg'].tolist()
        
        a6 = []
        d6 = []
        d5 = []
        d4 = []
        d3 = []
        d2 = []
        d1 = []
        small = []
        large = []
        
        for wave in ecg_waves:
            coeffs = wavedec(wave, 'sym6', level=6)
            a6.append(coeffs[-0])
            d6.append(coeffs[-1])
            d5.append(coeffs[-2])
            d4.append(coeffs[-3])
            d3.append(coeffs[-4])
            d2.append(coeffs[-5])
            d1.append(coeffs[-6])
            
            small.append(waverec(coeffs[:-6] + [None] * 6, 'sym6'))
            
            coeffs[-0] = np.zeros_like(coeffs[-0])
            coeffs[-1] = np.zeros_like(coeffs[-1])
            coeffs[-2] = np.zeros_like(coeffs[-2])
            coeffs[-3] = np.zeros_like(coeffs[-3])
            coeffs[-6] = np.zeros_like(coeffs[-6])
            
            large.append(waverec(coeffs, 'sym6'))
        
        database['coefficient a6'] = a6
        database['coefficient d6'] = d6
        database['coefficient d5'] = d5
        database['coefficient d4'] = d4
        database['coefficient d3'] = d3
        database['coefficient d2'] = d2
        database['coefficient d1'] = d1
        database['small frequencies'] = small
        database['large frequencies'] = large
        # Multilevel discrete decomposition of ECG waves for compression and noise reduction.
        
        
        
        return database

class PeakExtractor():
    
    def __init__(self, c):
        self.c = c
    
    def fit(self, database):
        ecg_waves = database['large frequencies'].tolist()
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
            ecg = subject['large frequencies']
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
    
class QRSFeatureExtractor():
    
    def __init__(self, c):
        self.c = c
    
    def fit(self, database):
        qr_heights = []
        rs_heights = []
        qs_lengths = []
        
        for wave in database['wavelet']:
            high, pos1 = find_peaks(wave, height=self.c)
            low, pos2 = find_peaks(-wave, height=self.c)
            
            pos1 = pos1['peak_heights']
            pos2 = pos2['peak_heights']

            q = pos2[0]
            r = pos1[0]
            s = pos2[-1]
            
            q2 = low[0]
            s2 = low[-1]
            
            q = -q
            q = abs(q)
            
            qr_height = (q + r)
            rs_height = (r + s)
            qs_length = (q + s)
            
            qr_heights.append(qr_height)
            rs_heights.append(rs_height)
            qs_lengths.append(qs_length)
            
        return qr_heights, rs_heights, qs_lengths

class IntervalLengthExtractor():
    
    def __init__(self):
        pass
        
    def fit(self, database):
        lengths = []
        
        for wave in database['wavelet']:
            length = len(wave)
            lengths.append(length)
            
        return lengths
    
class TimeGapExtractor():
    
    def __init__(self, c):
        self.c = c
    
    def fit(self, database):
        timegaps = []
        wl = database['wavelet']
        
        for i, wave in enumerate(wl):     
            if i == 0:
                timegap = np.nan
                timegaps.append(timegap)
            
            else:
                previous_wave = wl[i - 1]
                previous_peaks, pos1 = find_peaks(-previous_wave, height=self.c)
                next_peaks, pos2 = find_peaks(-wave, height=self.c)
                
                pos1 = pos1['peak_heights']
                pos2 = pos2['peak_heights']
                
                q1 = previous_peaks[0]
                s1 = previous_peaks[-1]
                q2 = next_peaks[0]
                s2 = next_peaks[-1]
                
                end = len(previous_wave)
                post_qrs = previous_wave[int(s1):int(end)]
                post_qrs = len(post_qrs)
                pre_qrs = len(wave[int(0):int(q2)])
                timegap = pre_qrs + post_qrs
                timegaps.append(timegap)
                
        return timegaps

class SmallFrequencySeparator():
    
    def __init__(self):
        pass
        
    def fit(self, database):  
        p_waves = []
        t_waves = []
        freqs = database['small frequencies']
        
        for wave in freqs:
            pt, pos_pt = find_peaks(-wave, 5)
            p_waves.append(wave[0: pt[0]])
            plt.title('0')
            plt.plot(wave[0: pt[0]])
            plt.show()
            
            for i, pts in enumerate(pt, 1):
                try:
                    print(i)
                    separation = wave[pt[i-1]:pt[i]]
                    plt.title(i)
                    plt.plot(separation)
                    plt.show()
                    
                    startpoint = separation[0]
                    endpoint = separation[-1]
                    pt2, pos_pt2 = find_peaks(-separation, -5)
                    print("pt2 before: " + str(pt2))
                    if pt2.size > 1:
                        pt2 = pt2[-1]
                    
                    elif pt2.size == 0:
                        break
                    
                    midpoint = int(pt2)
                    print("pt2 now: " + str(midpoint))
                    #endpoint = int(endpoint)
                    #print(midpoint)
                    #t_waves.append(separation[0:midpoint])
                    #p_waves.append(separation[midpoint:endpoint])
                except IndexError:
                    print("An anomaly occurred with data subject: " + str(i))
                    print(IndexError)
                finally:
                    continue
            #return p_waves, t_waves
                
        
        
            
        
    