# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:34:18 2020

@author: Oliver
"""
from scipy.signal import lfilter, filtfilt, savgol_filter
class BaselineNoiseRemover():
    
    def __init__(self, c):
        self.c = c

# DC Notch filter to remove baseline noise

    def fit(self, sig):
        b = [1, -1];
        a = [1, self.c];
        filt = filtfilt(b, a, sig);
        filt = savgol_filter(filt, 25, 2)
        return filt
