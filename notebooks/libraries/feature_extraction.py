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
        def slope(x, w, c):
            slopes = []
            freqs = np.zeros_like(x)
            for i in range(0, x.size):
                try:
                    freqs[i:i+w] = x[i:i+w]
                    y1 = x[i]
                    y2 = x[i + w]
                    x1 = 0
                    x2 = 0 + w
        
                    #print("y1: " + str(y1) + "y2: " + str(y2) + "x1: " + str(x1) + "x2: " + str(x2))
        
                    slope = (y2 - y1)/(x2 - x1)
                    slopes.append(slope)
                    if slope > c or slope < -c:
                        freqs[i:i+w] = np.nan
                        
                    i = i + w
                    
                except IndexError:
                    break
            
            return freqs[~np.isnan(freqs)], slopes
    
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
            coeffs = wavedec(wave, 'db4', level=6)
            a6.append(coeffs[-0])
            d6.append(coeffs[-1])
            d5.append(coeffs[-2])
            d4.append(coeffs[-3])
            d3.append(coeffs[-4])
            d2.append(coeffs[-5])
            d1.append(coeffs[-6])
            
            coeffs[-0] = np.zeros_like(coeffs[-0])
            coeffs[-1] = np.zeros_like(coeffs[-1])
            coeffs[-2] = np.zeros_like(coeffs[-2])
            #coeffs[-3] = np.zeros_like(coeffs[-3])
            #coeffs[-5] = np.zeros_like(coeffs[-5])
            coeffs[-6] = np.zeros_like(coeffs[-6])
            
            large.append(waverec(coeffs, 'db4'))
            
            #coeffs[-1] = np.zeros_like(coeffs[-1])
            #coeffs[-2] = np.zeros_like(coeffs[-2])
            #coeffs[-3] = np.zeros_like(coeffs[-3])
            #coeffs[-4] = np.zeros_like(coeffs[-4])
            #coeffs[-5] = np.zeros_like(coeffs[-5])
            #coeffs[-6] = np.zeros_like(coeffs[-6])
            
            #small.append(waverec(coeffs, 'db4'))
            
            #coeffs = wavedec(wave, 'db4', level=6)
            
            smallfreqs, s = slope(wave, 15, 3)
            smallfreqs = savgol_filter(smallfreqs, 43, 7)
            small.append(smallfreqs)
            
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
            peaks, position = find_peaks(wave, height=self.c, distance=50)
            peaks_list.append(peaks)
            peak_position.append(position)
            
        database['peaks'] = peaks_list
        database['peak position'] = peak_position
        
        return database
    
class SmallPeakExtractor():
    
    def __init__(self):
        pass
    
    def fit(self, database):
        small_waves = database['small frequencies'].tolist()
        small_peaks = []
        small_peaks_position = []
        
        for wave in small_waves:
            peaks, position = find_peaks(wave, -5, distance=30)
            small_peaks.append(peaks)
            small_peaks_position.append(position)
            
        database['small peaks'] = small_peaks
        database['small peak position'] = small_peaks_position
        
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

class SmallMidpointExtractor():
    def __init__(self):
        pass
    
    def midpoint(x1, y1, x2, y2):
        x = (x1 + x2) / 2
        y = (y1 + y2) / 2
        x = round(x)
        y = round(y)
        return x, y

    def find_midpoints(peaks, pos):
            
        midpoints_x = np.ones(len(peaks) - 1, dtype=float)
        midpoints_y = np.ones(len(peaks) - 1, dtype=float)
    
        for i in range(len(midpoints_x)):
            x1 = peaks[i]
            x2 = peaks[i + 1]
            y1 = pos[i]
            y2 = pos[i + 1]
            x, y = SmallMidpointExtractor.midpoint(x1, y1, x2, y2)
            midpoints_x[i] = x
            midpoints_y[i] = y       
        return np.vstack([midpoints_x.astype(int), midpoints_y.astype(int)])
    
    def fit(self, database):
        small_midpoints = []
        small_peaks = database['small peaks']
        small_positions = database['small peak position']
        
        for peak, pos in zip(small_peaks, small_positions):
            try:
                #print(pos)
                midpoints = SmallMidpointExtractor.find_midpoints(peak, pos['peak_heights'])
                small_midpoints.append(midpoints)
            except ValueError:
                small_midpoints.append(np.nan)
        
        database['small midpoints'] = small_midpoints
        database = database.dropna()
        return database

class SmallMidpointExtractor2():
    
    def __init__(self):
        pass
    def fit(self, database):
        small_waves = database['small frequencies'].tolist()
        small_peaks = []
        small_peaks_position = []
            
        for wave in small_waves:
            peaks, position = find_peaks(-wave, -20, distance=30)
            small_peaks.append(peaks)
            small_peaks_position.append(position)
            
        database['small midpoints2'] = small_peaks
        database['small midpoints2 pos'] = small_peaks_position
    
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
    
class SmallWaveletSeparator():
    
    def __init__(self):
        pass
    
    def fit(self, database, new_database):
        for i in range(0, len(database)):
            subject = database.iloc[i]
            midpoint = subject['small midpoints']
            ecg = subject['small frequencies']
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
    
class SmallWaveletSeparator2():
    
    def __init__(self):
        pass
    
    def fit(self, database, new_database):
        for i in range(0, len(database)):
            subject = database.iloc[i]
            midpoint = database['small midpoints2']
            condition = subject['condition']
            name = subject['name']
            x = midpoint
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
        
        
        for i, wave in enumerate(database['wavelet']):
            try:
                high, pos1 = find_peaks(wave, height=self.c)
                low, pos2 = find_peaks(-wave, height=10)
                
                pos1 = pos1['peak_heights']
                pos2 = pos2['peak_heights']
    
                q = pos2[0]
                r = pos1[0]
                s = pos2[-1]
                
                q2 = low[0]
                r2 = high[0]
                s2 = low[-1]
                
                q = -q
                q = abs(q)
                
                qr_height = (q + r)
                rs_height = (r + s)
                qs_length = (s2 - q2)
                
                qr_heights.append(qr_height)
                rs_heights.append(rs_height)
                qs_lengths.append(qs_length)
                
            except IndexError:
                qr_heights.append(np.nan)
                rs_heights.append(np.nan)
                qs_lengths.append(np.nan)
                
                continue
                
        print(len(qr_heights))
        print(len(rs_heights))
        print(len(qs_lengths))
        
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
        
        for i, this_wave in enumerate(wl):    
            try:
                next_wave = wl[i + 1]
           
                troughs_this, pos_this = find_peaks(-this_wave, self.c)
                troughs_next, pos_next = find_peaks(-next_wave, self.c)
                
                pos_this = pos_this['peak_heights']
                pos_next = pos_next['peak_heights']
                
                s1 = troughs_this[1]
                q2 = troughs_next[0]
                
                q2 = s1 + q2
                
                timegap = q2 - s1
                timegaps.append(timegap)
            except KeyError:
                timegaps.append(np.nan)
                continue
            except IndexError:
                timegaps.append(np.nan)
                continue
             
        return timegaps

class SmallFrequencySeparator():
    
    def __init__(self):
        pass
        
    def fit(self, database):  
        db_pt = pd.DataFrame(columns=['condition', 'name', 'ecg', 'no'])
        small_waves = []
        p_waves = []
        t_waves = []
        names = []
        cond = []
        no = []
        
        for i in range(0, database.size):
            #print(i)
            try:
                r = database.iloc[i]
            except IndexError:
                continue
            wave = r['small frequencies']
            name = r['name']
            condition = r['condition']
            pt, pos_pt = find_peaks(-wave, 12, distance=20)
            small_waves.append(wave[0:pt[0]])
            
            for i in range(2, pt.size - 2):
                co = 0
                a = pt[i - 1]
                b = pt[i]
                i = i + 1
                co = co + 1
                
                c = wave[a:b]
                size = c.size
                if size < 50 or np.all((c > 0)):
                    continue
                small_waves.append(c)
                names.append(name)
                cond.append(condition)
                no.append(co)
        
        try:
            db_pt['ecg'] = small_waves
            db_pt['condition'] = cond
            db_pt['name'] = name
            db_pt['no'] = no
        except ValueError:
            return db_pt
            
class PTLabeller():
    
    def __init__(self):
        pass
        
    def fit(self, database):
        types = []
        smallfreqs = database['wavelet']
        
        for s in smallfreqs:
            start = s[0]
            end = s[-1]
            if start > end:
                types.append('P')
            else:
                types.append('T')
        database['type'] = types
        return database
            
class PTFeatureExtractor():
    
    def __init__(self):
        pass

    def fit(self, database):
        size = []
        inclineheight = []
        declineheight = []
        inclinelength = []
        declinelength = []
        smallfreqs = database['wavelet'].tolist()
        
        for s in smallfreqs:
            
            length = s.size
            size.append(length)
            
            start = s[0]
            end = s[-1]
            s_peak, s_pos = find_peaks(s, -20)
            s_pos = s_pos['peak_heights']
            
            pointer = np.where(s_peak == np.amax(s_peak))
            s_peak = np.amax(s_peak)
            s_pos = s_pos[pointer]
            
            incline = s[0:s_peak]
            decline = s[s_peak:-1]
                
            inclineheight.append(s_pos - start)
            declineheight.append(end - s_pos)
            inclinelength.append(incline.size)
            declinelength.append(decline.size)
         
            
        print(len(size))
        print(len(inclineheight))
        print(len(declineheight))
        print(len(inclinelength))
        print(len(declinelength))
        database['size'] = size
        database['incline height'] = inclineheight
        database['decline height'] = declineheight
        database['incline length'] = inclinelength
        database['decline length'] = declinelength
        
        return database
        
        