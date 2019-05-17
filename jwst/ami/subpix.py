#! /usr/bin/env python

import numpy as np

def weightpixels(array, weightarray):
    """
    Short Summary
    -------------
    Weight.... wait for Deep's input (I asked)

    Parameters
    ----------
    array: 3D float array
        oversampled model       
 
    weightarray: 2D float array
        pixel weight array used to weight

    if np.shape(weightarray)[0] !=np.shape(weightarray)[1]:
        raise ValueError("Pixel Weight Array Is Not Square")

    oversample = np.shape(weightarray)[0]
    shapex = np.shape(array)[0]
    shapey = np.shape(array)[1]
    b = array.copy()
    b = b.reshape(shapex//oversample, oversample, shapey//oversample, oversample)
    d = b.copy()

    for i in range(np.shape(weightarray)[0]):
        for j in range(np.shape(weightarray)[1]):
            d[:,i,:,j] = d[:,i,:,j]* weightarray[i,j]
   """

    # e.g 1 in the center, 0.8 on the edges:
    d[:,1, :, 1] = d[:,1, :, 1]

    d[:,0, :, :] = 0.8 * d[:,0, :, :] 
    d[:,2, :, :] = 0.8 * d[:,2, :, :] 
    d[:,:, :, 0] = 0.8 * d[:, :, :, 0] 
    d[:,:, :, 2] = 0.8 * d[:, :, :, 2] 

    for i in range(shapex):
        for j in range(shapey):
            d[i,:,j,:] = b[i,:,j,:] * weightarray

    return d.sum(-1).sum(1)
