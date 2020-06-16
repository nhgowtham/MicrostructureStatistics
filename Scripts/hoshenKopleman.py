import numpy as np
import matplotlib.pyplot as plt
import pymks
import numpy as np
from PIL import Image
from Scripts import SpatialCorrelations as corr
import matplotlib.pyplot as plt
import pandas as pd
import os
from Scripts import velocityCalculations as vel

def unionArray(arr,eq1,eq2):
    arr_ = arr
    val1 = eq1
    val2 = eq2
    while arr_[val2]!=val2:
         val2 = arr_[val2]
    arr_[eq1] = val2
    return arr_

def hoshenKoplemanLabels(img_):
    properLabels = []
    for i in range(200000):
        properLabels.append(i)
        
        A = []
    B = []
    largestLabel = 0
    labels = np.zeros(img_.shape).astype('int32')
    for j in range(img_.shape[0]):
        for i in range(img_.shape[1]):
            if img_[j][i] == 1:
                if i>0 and j>0:
                    left = img_[j][i-1]
                    above = img_[j-1][i]
        
                    if left==0 and above== 0:
                        largestLabel=largestLabel+1
                        labels[j][i] = largestLabel
                    
                    if left==1 and above== 0:
                        labels[j][i] = properLabels[labels[j][i-1]]
                    if left==0 and above== 1:
                        labels[j][i] = properLabels[labels[j-1][i]]
                    if left==1 and above== 1:
                    
                        if labels[j][i-1]!= labels[j-1][i]:
                            if labels[j][i-1]> labels[j-1][i]:
                                A.append(labels[j][i-1])
                                B.append(labels[j-1][i])
                                properLabels=unionArray(properLabels,int(labels[j][i-1]),int(labels[j-1][i]))
                                labels[j][i] = properLabels[int(labels[j][i-1])]
                            if labels[j][i-1]< labels[j-1][i]:
                                A.append(int(labels[j-1][i]))
                                B.append(int(labels[j][i-1]))
                                properLabels=unionArray(properLabels,int(labels[j-1][i]),int(labels[j][i-1]))
                                labels[j][i] = properLabels[int(labels[j][i-1])]
                        else:
                            labels[j][i] = properLabels[labels[j][i-1]]
                        
                if i == 0 and j == 0:
                    if img_[j][i] == 1:
                        largestLabel=largestLabel+1
                        labels[j][i] = largestLabel
                if i == 0 and j>0:
                    above = img_[j-1][i]
                    if img_[j][i] == 1:
                        if above==1:
                            labels[j][i] = properLabels[labels[j-1][i]]
                        if above==0:
                            largestLabel=largestLabel+1
                            labels[j][i] = largestLabel
                if i>0 and j==0:
                    left = img_[j][i-1]
                    if img_[j][i] == 1:
                        if left==1:
                            labels[j][i] = properLabels[labels[j][i-1]]
                    
                        if left==0 :
                            largestLabel=largestLabel+1
                            labels[j][i] = largestLabel
    for i in range(labels.shape[0]):
        for j in range(labels.shape[1]):
            labels[i][j] = properLabels[labels[i][j]]
    
    L = {}
    a = np.unique(labels).tolist()
    for i in a:
        L[i] = a.index(i)
    
    for i in range(labels.shape[0]):
        for j in range(labels.shape[1]):
            labels[i][j] = L[labels[i][j]]
    
    
    return labels
                

    
    
def precipitateCentres(labelImage, labelNumber):
    img = (labelImage==labelNumber)*1
    row =0
    column =0
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if img[i][j]==1:
                row = row+i
                column = column+j
    
    return int(row/np.sum(img)) , int(column/np.sum(img))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
    