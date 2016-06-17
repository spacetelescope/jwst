#
#  Module for  applying straylight correction.
# The routine correct_MRS applues a straylight correction to MRS science
# slope images. The straylight mask contains 0's for science regions and
# 1's for gaps between slices.
#
# Remaining questions:
# Need to read in the DQ array and find out where there are bad pixel or
# slope not defined. We do not want to smooth data that is bad or poorly
# determined.


import numpy as np
import logging
import math
from jwst import datamodels
from astropy.convolution import convolve, Box2DKernel
#from matplotlib import pyplot as plt


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def correct_MRS(input_model, straylight_model):
    """
    Short Summary
    -------------
    Corrects the MIRI MRS data for straylight

    Parameter
    ----------
    input_model: data model object
        science data to be corrected

    straylight_model: holds the straylight mask for the correction

    sci_ngroups:int
        number of groups in input data


    Returns
    -------
    output: data model object
        straylight-subtracted science data

    """

    # Save some data parameterss for easy use later
    #nrows   = input_model.data.shape[0]
    #ncols = input_model.data.shape[1]

    nrows, ncols = input_model.data.shape

    mask = straylight_model.data
 #   plt.imshow(mask)
 #   plt.show()

    x = float('nan')
    # The straylight mask has values of 1 and 0. The straylight task uses data
    # in the between the slices of the MRS data to correct the science data
    # in the slices. In the mask the pixels found in the gaps (between the
    #slices) have a value of 1 and the science pixels have a value of 0.


    # Create output as a copy of the input science data model
    # sci_mask is the input science image * mask
    # science regions  = 0 (reference pixel are  also = 0)

    output = input_model.copy() # this is used in algorithm to
    # find the straylight correction.
    output2_data = input_model.data.copy()
    # final result = output2_data - straylight correction
    # Output2 is the orginal (no NANs have been removed)

    #_________________________________________________________________
    # is there are nans remove them because they mess up the correction
    index_inf = np.isinf(output.data).nonzero()
    index_nan = np.isnan(mask).nonzero()

    output.data[index_inf] = 0.0
    mask[index_inf] = 0

    output.data[index_nan] = 0.0
    mask[index_nan] = 0
    #_________________________________________________________________

    sci_mask = output.data * mask
    straylight_image = output.data * 0.0

    #sci_mask now contains 0's in science regions of detector.

    debug_row = 1
    #print 'input data ',input_model.data[debug_row,0:50]
    #print 'data',input_model.data[debug_row,0:50] #row 11 columns 1 to 49
    #print 'mask ',mask[debug_row,0:50] #row 11 columns 1 to 49
    #print 'mask data',sci_mask[debug_row,0:50]

    #We Want Sci mask smoothed for GAP region with 3 X 3 box car filter
    #
    # Handle edge cases for boxcar smoothing, by determining the
    # boxcar smoothing of the mask.

    sci_ave = convolve(sci_mask, Box2DKernel(3))
    mask_ave = convolve(mask, Box2DKernel(3))

    #print 'sci ave',sci_ave[debug_row,0:50]
    #print 'mask ave',mask_ave[debug_row,0:50]

    # catch /0 cases
    index = np.where(mask_ave == 0) # zero catches cases that would be #/0
                                   # near edges values are 0.3333 0.6667 1
    sci_ave[index] = 0
    mask_ave[index] = 1
    sci_smooth = sci_ave / mask_ave

    #print 'sci_smooth',sci_smooth[debug_row,0:50]
    x = np.arange(ncols)

    # Loop over each row (0 to 1023)
    for j in range(nrows):

        row_mask = mask[j, :]

        if(np.sum(row_mask) > 0):
            # find the locations of slice gaps
            #determine the data in the slice gaps
            yuse = sci_smooth[j, np.where(row_mask == 1)]

            #find the x locations of the slice gaps
            xuse = x[np.where(row_mask == 1)]
            #if(j ==1):
            #    print 'row ',j+1,yuse.shape,yuse
            #    print 'xuse',xuse,xuse.shape,type(xuse)

            #find data in same gap area
            nn = len(xuse)
            #dx difference in adjacent slice gaps pixels --> used
            # to find x limits of each gap
            dx = xuse[1:nn] - xuse[0:nn - 1]

            # Find the number of slice gaps in  row
            idx = np.asarray(np.where(dx > 1))
            ngroups = len(idx[0]) + 1

            # xlimits are the x values that mark limits of a slice gaps
            xlimits = np.zeros(ngroups + 1)
            i = 1
            xlimits[0] = xuse[0]
            for index in idx[0]:
                xlimits[i] = xuse[index]
                i = i + 1

            xlimits[ngroups] = xuse[nn - 1]#+1
            #if(j ==1):
            #    print '# xlimits',ngroups
            #    print 'xlimits',xlimits


            xg = np.zeros(ngroups)
            yg = np.zeros(ngroups)
            # loop over all slice gaps in a row
            # find the mean y straylight value for each slice gaps
            # also find the mean x value for this slice gap
            for i in range(ngroups):
                lower_limit = xlimits[i]
                upper_limit = xlimits[i + 1]
                igap = np.asarray(np.where(np.logical_and(xuse > lower_limit, xuse <= upper_limit)))
                ymean = np.mean(yuse[0, igap[0]])
                xmean = np.mean(xuse[igap[0]])
                yg[i] = ymean
                xg[i] = xmean
                #if(j ==1):
                #    print 'xg yg ',i,xg[i],yg[i]
        # else  entire row is zero
        else:
            xg = np.array([0, 1032])
            yg = np.array([0, 0])
       #using mean y value in slice gaps and x location of ymean
       #interpolate what straylight contribution based on yg and xg
       # for all points in row

        for k in x:
            if(x[k] >= xg[0] and x[k] <= xg[-1]):
                ynew = np.interp(x[k], xg, yg);
            else:
                if x[k] < xg[0]:
                    ynew = yg[0] + (x[k] - xg[0]) * (yg[1] - yg[0]) / (xg[1] - xg[0])
                elif x[k] > xg[-1]:
                    ynew = yg[-1] + (x[k] - xg[-1]) * (yg[-1] - yg[-2]) / (xg[-1] - xg[-2])
            straylight_image[j, k] = ynew

        # end loop over rows


    #straylight_image[np.where(straylight_image<0)] = 0
    straylight_image[straylight_image < 0] = 0

    #print 'straylight image',straylight_image[debug_row,0:40]

    # pull out the science region (1024 pixel/row) to do boxcar smoothing on

    simage = convolve(straylight_image, Box2DKernel(25))

    #print 'simage',simage[debug_row,0:50]
    # remove the straylight correction for the reference pixels
    simage[:, 1028:1032] = 0.0
    simage[:, 0:4] = 0.0

    output.data = output2_data - simage
    #print 'out',output.data[debug_row,0:50]

    return output
