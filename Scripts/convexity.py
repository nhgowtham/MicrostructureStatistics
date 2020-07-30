import numpy as np

def monteCarloConvexity(images):
    counter = 0
    counterPositive12 = 0
    A12=[]

    pointX=[]
    pointY=[]
    for i in range(images.shape[0]):
        for j in range(images.shape[1]):
            if images[i][j]==1:
                pointX.append(i)
                pointY.append(j)
           
    for i in range(50000):
        p1 = np.random.randint(0,len(pointX))
        p2 = np.random.randint(0,len(pointX))
        x1 = pointX[p1]
        y1 = pointY[p1]
        x2 = pointX[p2]
        y2 = pointY[p2] 
        x12 = int((x1+x2)/2)
        y12 = int((y1+y2)/2)

        counter=counter+1
        if images[x12][y12]==1:
            counterPositive12 = counterPositive12+1
            

    return (counterPositive12/counter)
            
    