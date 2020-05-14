import pymks
import numpy as np
from PIL import Image

from pymks.stats import autocorrelate
from pymks import PrimitiveBasis
from pymks.tools import draw_microstructures
from pymks.tools import draw_autocorrelations
from pymks.stats import crosscorrelate
from pymks.tools import draw_crosscorrelations
import pandas as pd
import matplotlib.pyplot as plt



def dat_to_numpy(image_path):
    a = pd.read_csv(image_path,sep=" ", header=None)
    row_img =np.max(a[0])+1
    column_img = np.max(a[1])+1
    img = np.zeros((row_img,column_img))
    a_numpy = a.to_numpy()
    for i in range(a_numpy.shape[0]):
        img[int(a_numpy[i][0])][int(a_numpy[i][1])] = a_numpy[i][2]
    
    return img
def binarize_image_dat(image_as_numpy):
    
    '''
    Input : Image in numpy array format of 0 or 0 cell values
    
    Output : Binary Image as numpy
    
    '''
    return (image_as_numpy>=0.5)*1

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
            
            p_white = img_binary[i][j]#+img_binary[ei][j]+img_binary[wi][j]+img_binary[i][ej]+img_binary[i][wj]
            p_white = p_white#/(5.0)
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
    img_back = np.abs(back_to_time)
    
    image =np.reshape(img_back,(m1.shape[0],m2.shape[1]))
    return np.fft.fftshift(image) 
    
def auto_corr_from_code(img_binary):
    
    '''
    Input : Binary image with shape (X,Y)
    
    Output : Images of auto correlation
    '''
    m_white,m_black =probability_matrix(img_binary)
    auto1 = get_2_point_statistics(m_white,m_white)
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
def radialDistribution(cross):
    '''
    Input : Correlation Vector
    
    Output : probability radially distributed
    '''   
    
    radiusVector = np.zeros(cross.shape)
    for i in range(cross.shape[0]):
        for j in range(cross.shape[1]):
            radiusVector[i][j] = np.sqrt((i-cross.shape[0]/2)**2 +(j-cross.shape[1]/2)**2)
    r = np.zeros(int(cross.shape[0]/2))
    for i in range(r.shape[0]):
        r1 = (radiusVector>=i)*1
        r2 = (radiusVector<=i+1)*1
        r_net = r1*r2
        sum_val =np.sum(r_net)
        convolved = np.multiply(r_net,cross)
        r[i] =np.sum(convolved)/sum_val
    r[0] = cross[int(cross.shape[0]/2)][int(cross.shape[1]/2)]
    return r

def giveAvailablePoints(cross, theta1, theta2, angularShift = 0):
    '''
    Input : Correlation Vector, 2 angles between which you need prob distribution
    
    Output : probability radially distributed between angles theta1 and theta2
    '''
    start = (theta1+angularShift)*np.pi/180
    end = (theta2+angularShift)*np.pi/180
    
    
    angleVector = np.zeros(cross.shape)
    for i in range(cross.shape[0]):
        for j in range(cross.shape[1]):
            if i ==cross.shape[0]/2:
                angleVector[i][j] = 0
            if j == cross.shape[1]/2:
                angleVector[i][j] = np.pi/2
            else:
                angleVector[i][j] = np.arctan((i-cross.shape[0]/2)/(cross.shape[1]/2-j))
    
    for i in range(angleVector.shape[0]):
        for j in range(angleVector.shape[1]):
            if i < cross.shape[0]/2 and j < cross.shape[1]/2 :
                angleVector[i][j] = angleVector[i][j] + np.pi
            if i > cross.shape[0]/2 and j < cross.shape[1]/2:
                angleVector[i][j] = angleVector[i][j] + np.pi
            if i > cross.shape[0]/2 and j > cross.shape[1]/2:
                angleVector[i][j] = angleVector[i][j] + np.pi*2
    
    
    r1 = (angleVector>=start)*1
    r2 = (angleVector<=end)*1
    r_net_angle = r1*r2
    
    
    radiusVector = np.zeros(cross.shape)
    for i in range(cross.shape[0]):
        for j in range(cross.shape[1]):
            radiusVector[i][j] = np.sqrt((i-cross.shape[0]/2)**2 +(j-cross.shape[1]/2)**2)
    r = np.zeros(int(cross.shape[0]/2))
    for i in range(r.shape[0]):
        r1 = (radiusVector>=i)*1
        r2 = (radiusVector<=i+1)*1
        r_net = r1*r2
        r_net = r_net * r_net_angle
        sum_val =np.sum(r_net)
        convolved = np.multiply(r_net,cross)
        r[i] =np.sum(convolved)/sum_val
    r[0] = cross[int(cross.shape[0]/2)][int(cross.shape[1]/2)]
    
        
    return r_net_angle, r

def levelSetWithoutSmoothening(image_as_numpy):
    '''
    INPUT : raw image as numpy
    
    Output : image as level set
    
    
    '''
    
    
    return image_as_numpy-0.5

def regionsOfBoundary(levelSetImage, c):
    
    
    '''
    Input : Level Set image
    
    Output : Only boundaries
    
    '''
    filter1 = (levelSetImage>-1*c)*1
    filter2 = (levelSetImage<c)*1
    filter3 =filter1*filter2
   
    return filter3

def volumeSmoothenedMicrostructure(levelSetImage):
    
    '''
    Input : Image already level set
    
    Output: Smoothened with 3x3 matrix
    
    '''
    smoothened_image = np.zeros(levelSetImage.shape)
    for i in range(levelSetImage.shape[0]):
        for j in range(levelSetImage.shape[1]):
            ei = i+1
            ej = j+1
            wi = i-1
            wj = j-1
        
            if(ei>levelSetImage.shape[0]-1):
                ei = ei - levelSetImage.shape[0]
            if(wi<0):
                wi = wi + levelSetImage.shape[0]
            if(ej>levelSetImage.shape[1]-1):
                ej = ej - levelSetImage.shape[1]
            if(wj<0):
                wj = wj + levelSetImage.shape[1]
            
            smoothened_image[i][j] = levelSetImage[i][j]+levelSetImage[i][ej]+levelSetImage[i][wj]+levelSetImage[wi][j]+levelSetImage[wi][ej]+levelSetImage[wi][wj]+levelSetImage[ei][j]+levelSetImage[ei][ej]+levelSetImage[ei][wj]
            smoothened_image[i][j]  = smoothened_image[i][j]/9
            
        
    return smoothened_image


def gradientMagnitude(img_,c):
    
    '''
    Input : raw image as numpy, size of boundary pixels to consider
    
    Output : Gradient Magnitudes Calculated on image
    
    '''
    img_ =levelSetWithoutSmoothening(img_)
    img_ = volumeSmoothenedMicrostructure(img_)
    boundary = regionsOfBoundary(img_,c)
    gradients = np.gradient(img_)
    x_grad = gradients[0]
    y_grad = gradients[1]
    gradient_abs = np.sqrt(x_grad*x_grad+y_grad*y_grad)
    gradient_abs = gradient_abs*boundary
    return gradient_abs

def dphidt(listOffiles, fileNumber, timestepForDt):
    '''
    Input : String list of files in order, file you want derivative for, time step for dt
    
    Output : dpho by dt of that image
    
    '''
    
    path = listOffiles[fileNumber]
    path_plus = listOffiles[fileNumber+timestepForDt]
    path_minus = listOffiles[fileNumber-timestepForDt]
    path_plus2 = listOffiles[fileNumber+2*timestepForDt]
    path_minus2 = listOffiles[fileNumber-2*timestepForDt]
    
    img_ = dat_to_numpy(path)
    img_ =levelSetWithoutSmoothening(img_)
    img_ = volumeSmoothenedMicrostructure(img_)
    
    img_plus = dat_to_numpy(path_plus)
    img_plus =levelSetWithoutSmoothening(img_plus)
    img_plus = volumeSmoothenedMicrostructure(img_plus)
    
    
    img_minus = dat_to_numpy(path_minus)
    img_minus =levelSetWithoutSmoothening(img_minus)
    img_minus = volumeSmoothenedMicrostructure(img_minus)
    
    img_plus2 = dat_to_numpy(path_plus2)
    img_plus2 =levelSetWithoutSmoothening(img_plus2)
    img_plus2 = volumeSmoothenedMicrostructure(img_plus2)
    
    
    img_minus2 = dat_to_numpy(path_minus2)
    img_minus2 =levelSetWithoutSmoothening(img_minus2)
    img_minus2 = volumeSmoothenedMicrostructure(img_minus2)
    
    
    return (-img_plus2+8*img_plus-8*img_minus+img_minus2)/(12*500*timestepForDt)
    
    
    
def velocity(listOffiles, fileNumber, timestepForDt):
    '''
    Input : String list of files in order, file you want derivative for, time step for dt

    Output : Velocity in boundary regions
    '''
    
    path = listOffiles[fileNumber]
    img_ = dat_to_numpy(path)
    gradient_magnitude =gradientMagnitude(img_,0.3)
    dphidt_ =dphidt(listOffiles, fileNumber, timestepForDt)
    
    img_ =levelSetWithoutSmoothening(img_)
    img_ = volumeSmoothenedMicrostructure(img_)
    img_boundary = regionsOfBoundary(img_,0.3)
    velocity = np.abs(np.divide(dphidt_*img_boundary,gradient_magnitude))
    velocity[np.isnan(velocity)] = 0
    
    return velocity

def plotHistogramOfVelocity(velocity):
    
    '''
    Input : velocity of image boundaries
    
    Output : Histogram
    
    '''
    listy =(velocity.flatten().tolist())
    A =[]
    for i in listy:
        if i!=0.0:
            A.append(i)
        
    A.sort(reverse = True)
    B = plt.boxplot(A)
    #B =plt.hist(A, bins=100, alpha=0.8,histtype=u'step', density=True, color = 'r')
    plt.show()
    
    return A





        
