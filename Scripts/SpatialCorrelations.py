import pymks
import numpy as np
from PIL import Image

from pymks.stats import autocorrelate
from pymks import PrimitiveBasis
from pymks.tools import draw_microstructures
from pymks.tools import draw_autocorrelations
from pymks.stats import crosscorrelate
from pymks.tools import draw_crosscorrelations

import matplotlib.pyplot as plt


def png_to_numpy(image_path, threshold, x_min=72, x_max=328, y_min=139, y_max=395):
    
    '''
    Input : A .png image path for B&W image
            x and y range values of image to be subsetted 
    
    Output: Numpy array of image
    
    
    '''
    
    image_in_numpy = np.array(Image.open(image_path).convert('L'))
    im_bool = image_in_numpy > threshold
    im_bin = (image_in_numpy > threshold) * 255
    
    # Subsetting to remove whitespace
    image_subset =im_bin[x_min:x_max,y_min:y_max]
    
    imagefinal = np.uint8(image_subset)
    
    
    return imagefinal

def binarize_image(image_as_numpy):
    
    '''
    Input : Image in numpy array format of 0 or 255 cell values
    
    Output : Binary Image as numpy
    
    '''
    return (image_as_numpy==255)*1
    
def auto_corr_from_pymks(img_binary):
    '''
    Input : Binary image with shape (X,Y)
    
    Output : Image of auto Correlation for both states
    '''
    immg = np.zeros((1,img_binary.shape[0],img_binary.shape[1]))
    immg[0] = img_binary
    p_basis = PrimitiveBasis(n_states=2)
    X_auto = autocorrelate(immg, p_basis, periodic_axes=(0, 1))
    correlations = [('black', 'black'), ('white', 'white')]
    draw_autocorrelations(X_auto[0], autocorrelations=correlations)
    
    return X_auto

def cross_corr_from_pymks(img_binary):
    '''
    Input : Binary image with shape (X,Y)
    
    Output : Image of cross correlation
    '''
    p_basis = PrimitiveBasis(n_states=2)
    immg = np.zeros((1,img_binary.shape[0],img_binary.shape[1]))
    immg[0] = img_binary
    X_cross = crosscorrelate(immg, p_basis, periodic_axes=(0, 1))
    correlations = [('black', 'white')]
    draw_crosscorrelations(X_cross[0], crosscorrelations=correlations)
    
    return X_cross

def probability_matrix(img_binary):
    '''
    Input : Binary Image
    
    Output : Probability of states on array
    
    
    '''
    
    m_white = np.zeros(img_binary.shape)
    m_black = np.zeros(img_binary.shape)


    for i in range(img_binary.shape[0]):
        for j in range(img_binary.shape[1]):
            ei = i+1
            ej = j+1
            wi = i-1
            wj = j-1
        
            if(ei>img_binary.shape[0]-1):
                ei = ei - img_binary.shape[0]
            if(wi<0):
                wi = wi + img_binary.shape[0]
            if(ej>img_binary.shape[1]-1):
                ej = ej - img_binary.shape[1]
            if(wj<0):
                wj = wj + img_binary.shape[1]
            
            p_white = img_binary[ei][j]+img_binary[wi][j]+img_binary[i][ej]+img_binary[i][wj]
            p_white = p_white/(4.0)
            p_black = 1.0 - p_white
        
            m_white[i][j] = p_white
            m_black[i][j] = p_black
        
    return m_white, m_black

def get_2_point_statistics(m1,m2):
    
    Fourier1 = np.fft.fft2(m1)
    Fourier2 = np.fft.fft2(m2)
    Fourier2 = Fourier2/(m1.shape[0]*m1.shape[1])
    Fourier2_conj =np.conjugate(Fourier2)
    NetFourier = np.multiply(Fourier1,Fourier2_conj)
    
    back_to_time =np.fft.ifft2(NetFourier)
    img_back = np.real(back_to_time)
    
    image =np.reshape(img_back,(m1.shape[0],m2.shape[1]))
    plt.imshow(np.fft.fftshift(image), cmap='jet')
    plt.colorbar()
    plt.grid(b=None)
    plt.show()
    return np.fft.fftshift(image) 
    
def auto_corr_from_code(img_binary):
    
    '''
    Input : Binary image with shape (X,Y)
    
    Output : Images of auto correlation
    '''
    m_white,m_black =probability_matrix(img_binary)
    
    
    print('white region : ')
    auto1 = get_2_point_statistics(m_white,m_white)
    
    print('black region : ')
    auto2 = get_2_point_statistics(m_black,m_black)
    
    
    return auto1, auto2

def cross_corr_from_code(img_binary):
    '''
    Input : Binary image with shape (X,Y)
    
    Output : Image of cross correlation
    '''
    m_white,m_black =probability_matrix(img_binary)
    
    auto1 = get_2_point_statistics(m_white,m_black)
    
    
    
    return auto1
    
    

    
    
    
        
    
    
    