# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class DynamicTimeWarper():
    
    def __init__(self):
        pass
    
    def fit(self, x, y):
        #Calculate a dynamic distance matrix
        distmat = np.zeros((len(x), len(y)))
        for i in range(0, len(x)):
            for j in range(0, len(y)):
                #Algorithm fills the matrix by calculating the dynamic distance between i and j
                if i == 0:
                    min_val = distmat[0, j - 1]
                elif j == 0:
                    min_val = distmat[i - 1, 0]
                elif (i and j) == 0:
                    min_val = distmat[0, 0]
                else:
                    min_val = min(distmat[i - 1, j - 1],
                                  distmat[i - 1, j],
                                  distmat[i, j - 1])
                    
                distmat[i, j] = abs(x[i] - y[j]) + min_val
        
        path = []
        position = []
        
        i = len(x) - 1
        j = len(y) - 1
        path.append(distmat[-1, -1])
        position.append([i, j])
        
        while (i or j) >= 0:
            conditions = np.zeros(3)
            conditions = ([distmat[i - 1, j - 1],
                           distmat[i - 1, j],
                           distmat[i, j - 1]])
        
            nextvalue = np.amin(conditions)
            path.append(nextvalue)
            condition = np.where(conditions == np.amin(conditions))
            if np.any(condition[0] == 0) == True:
                i = i - 1
                j = j - 1
            elif np.any(condition[0] == 1) == True:
                i = i - 1
            elif np.any(condition[0] == 2) == True:
                j = j - 1
            position.append([i, j])
            
        return distmat, (np.asarray(path), np.asarray(position))
    
    def example(self):
        x = np.arange(10)
        y = np.arange(10) + 2
        z = np.arange(10) + 4
        
        plt.plot(x)
        plt.plot(y)
        plt.show()
        signals = []
        
        signals.append(x)
        signals.append(y)
        signals.append(z)
        
        mat, path = DynamicTimeWarper.fit(x,y)
        positions = path[1]
        ax = plt.axes()
        ax = sns.heatmap(mat, annot=True, cmap='coolwarm_r', ax=ax)
        ax.set_title("Dynamic Time Warping example")
        ax.set(xlabel='ECG 1', ylabel='ECG 2')
        ax.invert_yaxis()
        
        a = positions[:,1]
        b = positions[:,0]
        ax.plot(a + 0.5, b + 0.5, color='blue', linewidth=5)
