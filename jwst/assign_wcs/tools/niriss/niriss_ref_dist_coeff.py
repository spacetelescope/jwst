#!/usrs/bin/env python

import numpy as np
import sys

from astropy.modeling import models
from asdf import AsdfFile

from ref_dist_coeff import dist_coeff

def ref_dist_coeff():
	""" Create distortion correction reference file.

	Parameters:
	___________
	coefficients: float
		coefficients

    The order of the input text files with coefficients is from det2sky and then sky2det.
	"""

	"Put coefficients into variables by reading in text files"
	coeffxr, coeffyr = np.loadtxt(sys.argv[1], skiprows=6, usecols=(1,2), unpack=True)
	print(coeffxr)
	print(coeffyr)

	coeffxi, coeffyi = np.loadtxt(sys.argv[2], skiprows=6, usecols=(1,2), unpack=True)
	print(coeffxi)
	print(coeffyi)

	"Transform coefficients to get undistorted, real, sky x."
	x0_0r = coeffxr[0]
	print('Coefficient for x0_0: ', x0_0r)
	x1_0r = coeffxr[1]
	print('Coefficient for x1_0: ', x1_0r)
	x0_1r = coeffxr[2]
	print('Coefficient for x0_1: ', x0_1r)
	x2_0r = coeffxr[3]
	print('Coefficient for x2_0: ', x2_0r)
	x1_1r = coeffxr[4]
	print('Coefficient for x1_1: ', x1_1r)
	x0_2r = coeffxr[5]
	print('Coefficient for x0_2: ', x0_2r)
	x3_0r = coeffxr[6]
	print('Coefficient for x3_0: ', x3_0r)
	x2_1r = coeffxr[7]
	print('Coefficient for x2_1: ', x2_1r)
	x1_2r = coeffxr[8]
	print('Coefficient for x1_2: ', x1_2r)
	x0_3r = coeffxr[9]
	print('Coefficient for x0_3: ', x0_3r)
	x4_0r = coeffxr[10]
	print('Coefficient for x4_0: ', x4_0r)
	x3_1r = coeffxr[11]
	print('Coefficient for x3_1: ', x3_1r)
	x2_2r = coeffxr[12]
	print('Coefficient for x2_2: ', x2_2r)
	x1_3r = coeffxr[13]
	print('Coefficient for x1_3: ', x1_3r)
	x0_4r = coeffxr[14]
	print('Coefficient for x0_4: ', x0_4r)

	"Transform coefficients to get undistorted, real, sky y."
	y0_0r = coeffyr[0]
	print('Coefficient for y0_0: ', y0_0r)
	y1_0r = coeffyr[1]
	print('Coefficient for y1_0: ', y1_0r)
	y0_1r = coeffyr[2]
	print('Coefficient for y0_1: ', y0_1r)
	y2_0r = coeffyr[3]
	print('Coefficient for y2_0: ', y2_0r)
	y1_1r = coeffyr[4]
	print('Coefficient for y1_1: ', y1_1r)
	y0_2r = coeffyr[5]
	print('Coefficient for y0_2: ', y0_2r)
	y3_0r = coeffyr[6]
	print('Coefficient for y3_0: ', y3_0r)
	y2_1r = coeffyr[7]
	print('Coefficient for y2_1: ', y2_1r)
	y1_2r = coeffyr[8]
	print('Coefficient for y1_2: ', y1_2r)
	y0_3r = coeffyr[9]
	print('Coefficient for y0_3: ', y0_3r)
	y4_0r = coeffyr[10]
	print('Coefficient for y4_0: ', y4_0r)
	y3_1r = coeffyr[11]
	print('Coefficient for y3_1: ', y3_1r)
	y2_2r = coeffyr[12]
	print('Coefficient for y2_2: ', y2_2r)
	y1_3r = coeffyr[13]
	print('Coefficient for y1_3: ', y1_3r)
	y0_4r = coeffyr[14]
	print('Coefficient for y0_4: ', y0_4r)

	"Transform coefficients to get distorted, ideal, detector x."
	x0_0i = coeffxi[0]
	print('Coefficient for x0_0: ', x0_0i)
	x1_0i = coeffxi[1]
	print('Coefficient for x1_0: ', x1_0i)
	x0_1i = coeffxi[2]
	print('Coefficient for x0_1: ', x0_1i)
	x2_0i = coeffxi[3]
	print('Coefficient for x2_0: ', x2_0i)
	x1_1i = coeffxi[4]
	print('Coefficient for x1_1: ', x1_1i)
	x0_2i = coeffxi[5]
	print('Coefficient for x0_2: ', x0_2i)
	x3_0i = coeffxi[6]
	print('Coefficient for x3_0: ', x3_0i)
	x2_1i = coeffxi[7]
	print('Coefficient for x2_1: ', x2_1i)
	x1_2i = coeffxi[8]
	print('Coefficient for x1_2: ', x1_2i)
	x0_3i = coeffxi[9]
	print('Coefficient for x0_3: ', x0_3i)
	x4_0i = coeffxi[10]
	print('Coefficient for x4_0: ', x4_0i)
	x3_1i = coeffxi[11]
	print('Coefficient for x3_1: ', x3_1i)
	x2_2i = coeffxi[12]
	print('Coefficient for x2_2: ', x2_2i)
	x1_3i = coeffxi[13]
	print('Coefficient for x1_3: ', x1_3i)
	x0_4i = coeffxi[14]
	print('Coefficient for x0_4: ', x0_4i)

	"Transform coefficients to get distorted, ideal, detector y."
	y0_0i = coeffyi[0]
	print('Coefficient for y0_0: ', y0_0i)
	y1_0i = coeffyi[1]
	print('Coefficient for y1_0: ', y1_0i)
	y0_1i = coeffyi[2]
	print('Coefficient for y0_1: ', y0_1i)
	y2_0i = coeffyi[3]
	print('Coefficient for y2_0: ', y2_0i)
	y1_1i = coeffyi[4]
	print('Coefficient for y1_1: ', y1_1i)
	y0_2i = coeffyi[5]
	print('Coefficient for y0_2: ', y0_2i)
	y3_0i = coeffyi[6]
	print('Coefficient for y3_0: ', y3_0i)
	y2_1i = coeffyi[7]
	print('Coefficient for y2_1: ', y2_1i)
	y1_2i = coeffyi[8]
	print('Coefficient for y1_2: ', y1_2i)
	y0_3i = coeffyi[9]
	print('Coefficient for y0_3: ', y0_3i)
	y4_0i = coeffyi[10]
	print('Coefficient for y4_0: ', y4_0i)
	y3_1i = coeffyi[11]
	print('Coefficient for y3_1: ', y3_1i)
	y2_2i = coeffyi[12]
	print('Coefficient for y2_2: ', y2_2i)
	y1_3i = coeffyi[13]
	print('Coefficient for y1_3: ', y1_3i)
	y0_4i = coeffyi[14]
	print('Coefficient for y0_4: ', y0_4i)

	"Create reference file using asdf formats."
	polyxr = models.Polynomial2D(4, c0_0=x0_0r, c1_0=x1_0r, c0_1=x0_1r, c2_0=x2_0r, c1_1=x1_1r, c0_2=x0_2r, c3_0=x3_0r, c2_1=x2_1r, c1_2=x1_2r, c0_3=x0_3r, c4_0=x4_0r, c3_1=x3_1r, c2_2=x2_2r, c1_3=x1_3r, c0_4=x0_4r)
	polyyr = models.Polynomial2D(4, c0_0=y0_0r, c1_0=y1_0r, c0_1=y0_1r, c2_0=y2_0r, c1_1=y1_1r, c0_2=y0_2r, c3_0=y3_0r, c2_1=y2_1r, c1_2=y1_2r, c0_3=y0_3r, c4_0=y4_0r, c3_1=y3_1r, c2_2=y2_2r, c1_3=y1_3r, c0_4=y0_4r)
	mapping = models.Mapping([0, 1, 0, 1])
	det2sky = mapping | polyxr & polyyr

	polyxi = models.Polynomial2D(4, c0_0=x0_0i, c1_0=x1_0i, c0_1=x0_1i, c2_0=x2_0i, c1_1=x1_1i, c0_2=x0_2i, c3_0=x3_0i, c2_1=x2_1i, c1_2=x1_2i, c0_3=x0_3i, c4_0=x4_0i, c3_1=x3_1i, c2_2=x2_2i, c1_3=x1_3i, c0_4=x0_4i)
	polyyi = models.Polynomial2D(4, c0_0=y0_0i, c1_0=y1_0i, c0_1=y0_1i, c2_0=y2_0i, c1_1=y1_1i, c0_2=y0_2i, c3_0=y3_0i, c2_1=y2_1i, c1_2=y1_2i, c0_3=y0_3i, c4_0=y4_0i, c3_1=y3_1i, c2_2=y2_2i, c1_3=y1_3i, c0_4=y0_4i)
	mapping = models.Mapping([0, 1, 0, 1])
	sky2det = mapping | polyxi & polyyi

	det2sky.inverse = sky2det

	f = AsdfFile()
	f.tree['model'] = det2sky
	f.tree['TITLE'] = "NIRISS CDP5 distortion reference data."
	f.tree['FILENAME'] = "niriss_ref_dist_coeff_image.asdf"
	f.tree['EXP_TYPE'] = "NIS_IMAGE"
	f.tree['REFTYPE'] = "DISTORTION"
	f.tree['AUTHOR'] = "Michael A. Wolfe, Alex Fullerton"
	f.tree['PEDIGREE'] = "GROUND"
	f.tree['INSTRUMENT'] = "NIRISS"
	f.tree['DETECTOR'] = "NIS"
	f.tree['USEAFTER'] = "Jan 01 2015 00:00:00.00"
	f.tree['DESCRIP'] = "This is a distortion correction reference file."
	f.tree['HISTORY'] = "This file is being delivered because it is the distortion correction file. This file was created from the python script niriss_ref_dist_coeff.py and data provied by Alex Fullerton. There is no documentation yet and the python script is titled: niriss_ref_dist_coeff.py. The data used to create the file can be found in the text files: NIS_Ideal2Real.txt & NIS_Real2Ideal.txt."
	f.write_to('niriss_ref_distortion_image.asdf')

	f = AsdfFile()
	f.tree['model'] = det2sky
	f.tree['TITLE'] = "NIRISS CDP5 distortion reference data."
	f.tree['FILENAME'] = "jwst_niriss_distortion_0001.asdf"
	f.tree['EXP_TYPE'] = "NIS_IMAGE"
	f.tree['REFTYPE'] = "DISTORTION"
	f.tree['AUTHOR'] = "Michael A. Wolfe, Alex Fullerton"
	f.tree['PEDIGREE'] = "GROUND"
	f.tree['INSTRUMENT'] = "NIRISS"
	f.tree['DETECTOR'] = "NIS"
	f.tree['USEAFTER'] = "Jan 01 2015 00:00:00.00"
	f.tree['DESCRIP'] = "This is a distortion correction reference file."
	f.tree['HISTORY'] = "This file is being delivered because it is the distortion correction file. This file was created from the python script niriss_ref_dist_coeff.py and data provied by Alex Fullerton. There is no documentation yet and the python script is titled: niriss_ref_dist_coeff.py. The data used to create the file can be found in the text files: NIS_Ideal2Real.txt & NIS_Real2Ideal.txt."
	f.write_to('jwst_niriss_distortion_0001.asdf')

	print('Done, done, and done with creating the disrtortion correction reference file. Git the bleep along little doggies!!!!!')

	dist_coeff()

if __name__ == '__main__':
	ref_dist_coeff()