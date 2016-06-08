"""
MODULE RECEIVES A NUMARRAY OBJECT AS INPUT, REPRESENTING A 2-D IMAGE,
AND COMPUTES THE ABSOLUTE DERIVATE OF THAT IMAGE.  A NUMARRY OBJECT
IS RETURNED BY :py:func:`quickderiv`.

:Authors: CHRISTOPHER HANLEY

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
#
# VERSION:
#   Version 0.1.0: created -- CJH
#
from __future__ import absolute_import, division # confidence high

# IMPORT EXTERNAL MODULES
import numpy as np

def qderiv(array): # TAKE THE ABSOLUTE DERIVATIVE OF A NUMARRY OBJECT
    """Take the absolute derivate of an image in memory."""

    #Create 2 empty arrays in memory of the same dimensions as 'array'
    tmpArray = np.zeros(array.shape,dtype=np.float64)
    outArray = np.zeros(array.shape, dtype=np.float64)

    # Get the length of an array side
    (naxis1,naxis2) = array.shape
    #print "The input image size is (",naxis1,",",naxis2,")."

#Main derivate loop:
    #Shift images +/- 1 in Y.
    for y in range(-1,2,2):
        if y == -1:
            #shift input image 1 pixel right
            tmpArray[0:(naxis1-1),1:(naxis2-1)] = array[0:(naxis1-1),0:(naxis2-2)]
            #print "Y shift = 1"
        else:
            #shift input image 1 pixel left
            tmpArray[0:(naxis1-1),0:(naxis2-2)] = array[0:(naxis1-1),1:(naxis2-1)]
            #print "Y shift = -1"
        #print "call _absoluteSubtract()"
        (tmpArray,outArray) = _absoluteSubtract(array,tmpArray,outArray)

    #Shift images +/- 1 in X.
    for x in range(-1,2,2):
        if x == -1:
            #shift input image 1 pixel right
            tmpArray[1:(naxis1-1),0:(naxis2-1)] = array[0:(naxis1-2),0:(naxis2-1)]
            #print "X shift = 1"
        else:
            #shift input image 1 pixel left
            tmpArray[0:(naxis1-2),0:(naxis2-1)] = array[1:(naxis1-1),0:(naxis2-1)]
            #print "X shift = -1"
        #print "call _absoluteSubtract()"
        (tmpArray,outArray) = _absoluteSubtract(array,tmpArray,outArray)



    return outArray.astype(np.float32)


def _absoluteSubtract(array,tmpArray,outArray):
    #subtract shifted image from imput image
    tmpArray = array - tmpArray
    #take the absolute value of tmpArray
    tmpArray = np.fabs(tmpArray)
    #save maximum value of outArray or tmpArray and save in outArray
    outArray = np.maximum(tmpArray,outArray)
    #zero out tmpArray before reuse
    tmpArray = tmpArray * 0.

    return (tmpArray,outArray)

# END MODULE
