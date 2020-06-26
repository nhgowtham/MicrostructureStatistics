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
    for i in range(20000):
        properLabels.append(i)
        
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
                                properLabels=unionArray(properLabels,int(labels[j][i-1]),int(labels[j-1][i]))
                                labels[j][i] = properLabels[int(labels[j][i-1])]
                            if labels[j][i-1]< labels[j-1][i]:
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
                            
    
            
            
    # Restarting the algorithm from scratch ro account for edge values   
        
    
    for j in range(img_.shape[0]):
        for i in range(img_.shape[1]):
            if img_[j][i] == 1:
                
                # part 1
                if 1:
         
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
                        properLabels=unionArray(properLabels,int(labels[j][i]),int(labels[j][iminus1]))
                    if left==0 and above== 1:
                        properLabels=unionArray(properLabels,int(labels[j][i]),int(labels[jminus1][i]))
                    if left==1 and above== 1:
                        if labels[j][iminus1]!= labels[jminus1][i]:
                            if 1:
                                properLabels=unionArray(properLabels,int(labels[j][iminus1]),int(labels[jminus1][i]))
                                properLabels=unionArray(properLabels,int(labels[j][i]),int(labels[jminus1][i]))
                        else:
                            properLabels=unionArray(properLabels,int(labels[j][i]),int(labels[jminus1][i]))
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
    def precipitateCentres(labelImage, labelNumber):
    img = (labelImage==labelNumber)*1
    m_y =[]
    m_x =[]
    x = []
    y = []

    for i in range(img.shape[0]):
        m_y.append(np.sum(img[i]))
        y.append(i/(img.shape[0]-1)*np.pi*2)
    for j in range(img.shape[1]):
        m_x.append(np.sum(img[:,j]))
        x.append(j/(img.shape[1]-1)*np.pi*2)
    m_x = np.array(m_x)
    m_y = np.array(m_y)
    x = np.array(x)
    y = np.array(y)
    

    alpha_x = np.cos(x)
    beta_x = np.sin(x)
    alpha_mean = np.sum(m_x*np.cos(x))/(np.sum(m_x))
    beta_mean = np.sum(m_x*np.sin(x))/(np.sum(m_x))
    theta_avg_x = np.arctan2(-1*beta_mean, -1*alpha_mean)+np.pi
    cog_x = img.shape[1]*theta_avg_x/(2*np.pi)
    
    alpha_y = np.cos(y)
    beta_y = np.sin(y)
    alpha_mean = np.sum(m_y*np.cos(y))/(np.sum(m_y))
    beta_mean = np.sum(m_y*np.sin(y))/(np.sum(m_y))
    theta_avg_y = np.arctan2(-1*beta_mean, -1*alpha_mean)+np.pi
    cog_y = img.shape[0]*theta_avg_y/(2*np.pi)
    
    
        
    return int(cog_y),int(cog_x)

def areaDistribution(labels):
    A =[]
    noOfLabels = np.max(labels)
    for i in range(1,noOfLabels+1):
        A.append(np.sum((labels==i)*1))
    return A