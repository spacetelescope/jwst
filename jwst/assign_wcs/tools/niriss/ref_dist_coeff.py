#!/usrs/bin/env python

import itertools
import numpy as np
import matplotlib.pyplot as plt
import sys

from astropy.table import Table, Column
from astropy.io import ascii
from astropy.modeling import models
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from asdf import AsdfFile

#------------------------------------------------------------------------------

def apply_coeff(x, y, cx, cy):
	""" transform coordinates

	Parameters:
	-----------
	x : int
	   x position
	y : int
       y position
    cx : list
       x coefficients
    cy : list
       y coefficients
    """

	x_out = cx[0] + cx[1]*x + cx[2]*y + cx[3]*x**2 + cx[4]*x*y + cx[5]*y**2 + cx[6]*x**3 + cx[7]*x**2*y + cx[8]*x*y**2 + cx[9]*y**3 + cx[10]*x**4 + \
	cx[11]*x**3*y + cx[12]*x**2*y**2 + cx[13]*x*y**3 + cx[14]*y**4
	y_out = cy[0] + cy[1]*x + cy[2]*y + cy[3]*x**2 + cy[4]*x*y + cy[5]*y**2 + cy[6]*x**3 + cy[7]*x**2*y + cy[8]*x*y**2 + cy[9]*y**3 + cy[10]*x**4 + \
	cy[11]*x**3*y + cy[12]*x**2*y**2 + cy[13]*x*y**3 + cy[14]*y**4

	return x_out, y_out

#------------------------------------------------------------------------------

def dist_coeff():
	""" Create distortion correction reference file.

	Parameters:
	___________
	coefficients: float
		coefficients

	The order of the input text files with coefficients is from det2sky and then sky2det.
	"""

	"Put coefficients into variables by reading in text files"
	coeffxr, coeffyr = np.loadtxt(sys.argv[1], skiprows=6, usecols=(1,2), unpack=True)
	print coeffxr
	print coeffyr

	coeffxi, coeffyi = np.loadtxt(sys.argv[2], skiprows=6, usecols=(1,2), unpack=True)
	print coeffxi
	print coeffyi

	"Transform coefficients to get undistorted, real, sky x."
	x0_0r = coeffxr[0]
	print 'Coefficient for x0_0: ', x0_0r
	x1_0r = coeffxr[1]
	print 'Coefficient for x1_0: ', x1_0r
	x0_1r = coeffxr[2]
	print 'Coefficient for x0_1: ', x0_1r
	x2_0r = coeffxr[3]
	print 'Coefficient for x2_0: ', x2_0r
	x1_1r = coeffxr[4]
	print 'Coefficient for x1_1: ', x1_1r
	x0_2r = coeffxr[5]
	print 'Coefficient for x0_2: ', x0_2r
	x3_0r = coeffxr[6]
	print 'Coefficient for x3_0: ', x3_0r
	x2_1r = coeffxr[7]
	print 'Coefficient for x2_1: ', x2_1r
	x1_2r = coeffxr[8]
	print 'Coefficient for x1_2: ', x1_2r
	x0_3r = coeffxr[9]
	print 'Coefficient for x0_3: ', x0_3r
	x4_0r = coeffxr[10]
	print 'Coefficient for x4_0: ', x4_0r
	x3_1r = coeffxr[11]
	print 'Coefficient for x3_1: ', x3_1r
	x2_2r = coeffxr[12]
	print 'Coefficient for x2_2: ', x2_2r
	x1_3r = coeffxr[13]
	print 'Coefficient for x1_3: ', x1_3r
	x0_4r = coeffxr[14]
	print 'Coefficient for x0_4: ', x0_4r

	"Transform coefficients to get undistorted, real, sky y."
	y0_0r = coeffyr[0]
	print 'Coefficient for y0_0: ', y0_0r
	y1_0r = coeffyr[1]
	print 'Coefficient for y1_0: ', y1_0r
	y0_1r = coeffyr[2]
	print 'Coefficient for y0_1: ', y0_1r
	y2_0r = coeffyr[3]
	print 'Coefficient for y2_0: ', y2_0r
	y1_1r = coeffyr[4]
	print 'Coefficient for y1_1: ', y1_1r
	y0_2r = coeffyr[5]
	print 'Coefficient for y0_2: ', y0_2r
	y3_0r = coeffyr[6]
	print 'Coefficient for y3_0: ', y3_0r
	y2_1r = coeffyr[7]
	print 'Coefficient for y2_1: ', y2_1r
	y1_2r = coeffyr[8]
	print 'Coefficient for y1_2: ', y1_2r
	y0_3r = coeffyr[9]
	print 'Coefficient for y0_3: ', y0_3r
	y4_0r = coeffyr[10]
	print 'Coefficient for y4_0: ', y4_0r
	y3_1r = coeffyr[11]
	print 'Coefficient for y3_1: ', y3_1r
	y2_2r = coeffyr[12]
	print 'Coefficient for y2_2: ', y2_2r
	y1_3r = coeffyr[13]
	print 'Coefficient for y1_3: ', y1_3r
	y0_4r = coeffyr[14]
	print 'Coefficient for y0_4: ', y0_4r

	"Transform coefficients to get distorted, ideal, detector x."
	x0_0i = coeffxi[0]
	print 'Coefficient for x0_0: ', x0_0i
	x1_0i = coeffxi[1]
	print 'Coefficient for x1_0: ', x1_0i
	x0_1i = coeffxi[2]
	print 'Coefficient for x0_1: ', x0_1i
	x2_0i = coeffxi[3]
	print 'Coefficient for x2_0: ', x2_0i
	x1_1i = coeffxi[4]
	print 'Coefficient for x1_1: ', x1_1i
	x0_2i = coeffxi[5]
	print 'Coefficient for x0_2: ', x0_2i
	x3_0i = coeffxi[6]
	print 'Coefficient for x3_0: ', x3_0i
	x2_1i = coeffxi[7]
	print 'Coefficient for x2_1: ', x2_1i
	x1_2i = coeffxi[8]
	print 'Coefficient for x1_2: ', x1_2i
	x0_3i = coeffxi[9]
	print 'Coefficient for x0_3: ', x0_3i
	x4_0i = coeffxi[10]
	print 'Coefficient for x4_0: ', x4_0i
	x3_1i = coeffxi[11]
	print 'Coefficient for x3_1: ', x3_1i
	x2_2i = coeffxi[12]
	print 'Coefficient for x2_2: ', x2_2i
	x1_3i = coeffxi[13]
	print 'Coefficient for x1_3: ', x1_3i
	x0_4i = coeffxi[14]
	print 'Coefficient for x0_4: ', x0_4i

	"Transform coefficients to get distorted, ideal, detector y."
	y0_0i = coeffyi[0]
	print 'Coefficient for y0_0: ', y0_0i
	y1_0i = coeffyi[1]
	print 'Coefficient for y1_0: ', y1_0i
	y0_1i = coeffyi[2]
	print 'Coefficient for y0_1: ', y0_1i
	y2_0i = coeffyi[3]
	print 'Coefficient for y2_0: ', y2_0i
	y1_1i = coeffyi[4]
	print 'Coefficient for y1_1: ', y1_1i
	y0_2i = coeffyi[5]
	print 'Coefficient for y0_2: ', y0_2i
	y3_0i = coeffyi[6]
	print 'Coefficient for y3_0: ', y3_0i
	y2_1i = coeffyi[7]
	print 'Coefficient for y2_1: ', y2_1i
	y1_2i = coeffyi[8]
	print 'Coefficient for y1_2: ', y1_2i
	y0_3i = coeffyi[9]
	print 'Coefficient for y0_3: ', y0_3i
	y4_0i = coeffyi[10]
	print 'Coefficient for y4_0: ', y4_0i
	y3_1i = coeffyi[11]
	print 'Coefficient for y3_1: ', y3_1i
	y2_2i = coeffyi[12]
	print 'Coefficient for y2_2: ', y2_2i
	y1_3i = coeffyi[13]
	print 'Coefficient for y1_3: ', y1_3i
	y0_4i = coeffyi[14]
	print 'Coefficient for y0_4: ', y0_4i

	"Generate ideal or detector coordinates."
	x_pix = np.arange(0, 2048, 1)
	y_pix = np.arange(0, 2048, 1)

	x_det_orig = []
	y_det_orig = []

	for x, y in itertools.product(x_pix, y_pix):
		x_det_orig.append(x)
		y_det_orig.append(y)

	print 'Done with first for loop!!!!'

	print len(x_det_orig)
	print len(y_det_orig)

	"To go from ideal, distorted to real, undistorted."

	x_sky = []
	y_sky = []

	#-------
	for x, y in zip(x_det_orig, y_det_orig):
		x_real, y_real = apply_coeff(x, y, coeffxr, coeffyr)
		x_sky.append(x_real)
		y_sky.append(y_real)
	#-------

	print ' '
	print 'Done with second for loop!!!!'

	print len(x_sky)
	print len(y_sky)
	x_sky = np.array(x_sky)
	y_sky = np.array(y_sky)
	print 'Above is for transformation to sky pixels.'

	"To go from real, undistorted to ideal, distorted."

	x_det = []
	y_det = []

	for x, y in zip(x_sky, y_sky):
		x_ideal, y_ideal = apply_coeff(x, y, coeffxi, coeffyi)
		x_det.append(x_ideal)
		y_det.append(y_ideal)

	print ' '
	print 'Done with third for loop!!!!'

	print len(x_det)
	print len(y_det)
	print 'Above is for transformation to detector pixels.'

	"Derive residuals from the complete transformation."
	x_det = np.array(x_det)
	y_det = np.array(y_det)
	x_det_orig = np.array(x_det_orig)
	y_det_orig = np.array(y_det_orig)

	resid_x = x_det - x_det_orig
	resid_y = y_det - y_det_orig
	print ' '
	print 'Residual in x: ', resid_x
	print np.min(resid_x), np.max(resid_x)
	print ' '
	print 'Residual in y: ', resid_y
	print np.min(resid_y), np.max(resid_y)

	'Open asdf reference file and transfrom from det to sky.'
	f = AsdfFile.open('niriss_ref_distortion_image.asdf')
	det2sky_trans= f.tree['model']
	sky_x, sky_y = det2sky_trans(x_det_orig, y_det_orig)

	'Transform from sky to det.'
	sky2det_trans = det2sky_trans.inverse
	inv_x, inv_y = sky2det_trans(sky_x, sky_y)

	'Calculate residuals from reference file.'
	x_resid_ref = inv_x - x_det_orig
	y_resid_ref = inv_y - y_det_orig

	'Calculate residuals of residuals.'
	x_resid_resid = np.absolute(resid_x - x_resid_ref)
	y_resid_resid = np.absolute(resid_y - y_resid_ref)

	print ' '
	print x_resid_resid
	print 'The minimum in x_resid_resid is: ', np.min(x_resid_resid)
	print 'The maximum in x_resid_resid is: ', np.max(x_resid_resid)
	print 'The mean in x_resid_resid is: ', np.mean(x_resid_resid)
	print 'The standard deviation in x_resid_resid is: ', np.std(x_resid_resid)
	print 'The median in x_resid_resid is: ', np.median(x_resid_resid)

	print ' '
	print y_resid_resid
	print 'The minimum in y_resid_resid is: ', np.min(y_resid_resid)
	print 'The maximum in y_resid_resid is: ', np.max(y_resid_resid)
	print 'The mean in y_resid_resid is: ', np.mean(y_resid_resid)
	print 'The standard deviation in y_resid_resid is: ', np.std(y_resid_resid)
	print 'The median in y_resid_resid is: ', np.median(y_resid_resid)

	print ' '
	print 'Done, Done, and Done. Git the bleep along little doggies.'

if __name__ == '__main__':
	dist_coeff()