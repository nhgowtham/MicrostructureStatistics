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
                
                # part 1
                if 1==1:
         
                    if i==0:
                        iminus1 = img_.shape[1]-1
                    else:
                        iminus1 = i-1
                    if j==0:  
                        jminus1 = img_.shape[0]-1
                    else:
                        jminus1 = j-1
                    
                    
                    
                    left = img_[j][iminus1]
                    above = img_[jminus1][i]
        
                    if left==0 and above== 0:
                        largestLabel=largestLabel+1
                        labels[j][i] = largestLabel
                    
                    if left==1 and above== 0:
                        labels[j][i] = properLabels[labels[j][iminus1]]
                    if left==0 and above== 1:
                        labels[j][i] = properLabels[labels[jminus1][i]]
                    if left==1 and above== 1:
                    
                        if labels[j][iminus1]!= labels[jminus1][i]:
                            if labels[j][iminus1]> labels[jminus1][i]:
                                A.append(labels[j][iminus1])
                                B.append(labels[jminus1][i])
                                properLabels=unionArray(properLabels,int(labels[j][iminus1]),int(labels[jminus1][i]))
                                labels[j][i] = properLabels[int(labels[j][iminus1])]
                            if labels[j][iminus1]< labels[jminus1][i]:
                                A.append(int(labels[jminus1][i]))
                                B.append(int(labels[j][iminus1]))
                                properLabels=unionArray(properLabels,int(labels[jminus1][i]),int(labels[j][iminus1]))
                                labels[j][i] = properLabels[int(labels[j][iminus1])]
                        else:
                            labels[j][i] = properLabels[labels[j][iminus1]]
                            
    # Restarting the algorithm from scratch ro account for edge values                     
    for j in range(img_.shape[0]):
        for i in range(img_.shape[1]):
            if img_[j][i] == 1:
                
                # part 1
                if 1==1:
         
                    if i==0:
                        iminus1 = img_.shape[1]-1
                    else:
                        iminus1 = i-1
                    if j==0:  
                        jminus1 = img_.shape[0]-1
                    else:
                        jminus1 = j-1
                    
                    
                    
                    left = img_[j][iminus1]
                    above = img_[jminus1][i]
                    
                    if left==1 and above== 0:
                        labels[j][i] = properLabels[labels[j][iminus1]]
                    if left==0 and above== 1:
                        labels[j][i] = properLabels[labels[jminus1][i]]
                    if left==1 and above== 1:
                    
                        if labels[j][iminus1]!= labels[jminus1][i]:
                            if labels[j][iminus1]> labels[jminus1][i]:
                                A.append(labels[j][iminus1])
                                B.append(labels[jminus1][i])
                                properLabels=unionArray(properLabels,int(labels[j][iminus1]),int(labels[jminus1][i]))
                                labels[j][i] = properLabels[int(labels[j][iminus1])]
                            if labels[j][iminus1]< labels[jminus1][i]:
                                A.append(int(labels[jminus1][i]))
                                B.append(int(labels[j][iminus1]))
                                properLabels=unionArray(properLabels,int(labels[jminus1][i]),int(labels[j][iminus1]))
                                labels[j][i] = properLabels[int(labels[j][iminus1])]
                        else:
                            labels[j][i] = properLabels[labels[j][iminus1]]
               
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
    