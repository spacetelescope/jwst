# Routines to read the D2C files

import sys
import numpy as np
import math
from .. import datamodels
from astropy.io import fits
from . import cube

#________________________________________________________________________________
def ReadD2CMapFile(self, channel, subchannel,):
    """
    Short Summary
    -------------
    Uses old style of distortion files. A 1033 X 1025 file containing the alpha,beta and lambda for
    each pixel corner. There is a  seperate fits file for each Detector and subchannel (so 6 files).

    Parameters
    ----------
    channel: Channel working with
    subchannel: subchannel working with

    Returns
    -------
    returns self.wcs holding a mapping information


    """
    if(channel == '1' or channel == '2'):
        if(subchannel == 'LONG'):
            self.wcs['Filename'] = 'MIRI_FM_SW_C_D2C_01.00.00.fits'
        elif (subchannel == 'MEDIUM'):
            self.wcs['Filename'] = 'MIRI_FM_SW_B_D2C_01.00.00.fits'
        elif (subchannel == 'SHORT'):
            self.wcs['Filename'] = 'MIRI_FM_SW_A_D2C_01.00.00.fits'
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elif(channel == '3' or channel == '4'):
        if(subchannel == 'LONG'):
            self.wcs['Filename'] = 'MIRI_FM_LW_C_D2C_01.00.00.fits'
        elif(subchannel == 'MEDIUM'):
            self.wcs['Filename'] = 'MIRI_FM_LW_B_D2C_01.00.00.fits'
        elif(subchannel == 'SHORT'):
            self.wcs['Filename'] = 'MIRI_FM_LW_A_D2C_01.00.00.fits'

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    f = fits.open(self.wcs['Filename'])

    print('Reading D2C file', self.wcs['Filename'])

    self.wcs['Wavelength'] = f[1].data
    self.wcs['Alpha'] = f[2].data
    self.wcs['SliceNo'] = f[3].data

    header = f[0].header
    if(channel == '1'):
        self.wcs['L_MIN'] = header['L_MIN1']
        self.wcs['L_MAX'] = header['L_MAX1']

        self.wcs['A_MIN'] = header['A_MIN1']
        self.wcs['A_MAX'] = header['A_MAX1']

        self.wcs['B_MIN'] = header['B_MIN1']
        self.wcs['B_DELTA'] = header['B_DEL1']
    elif(channel == '2'):
        self.wcs['L_MIN'] = header['L_MIN2']
        self.wcs['L_MAX'] = header['L_MAX2']

        self.wcs['A_MIN'] = header['A_MIN2']
        self.wcs['A_MAX'] = header['A_MAX2']

        self.wcs['B_MIN'] = header['B_MIN2']
        self.wcs['B_DELTA'] = header['B_DEL2']
    elif(channel == '3'):
        self.wcs['L_MIN'] = header['L_MIN3']
        self.wcs['L_MAX'] = header['L_MAX3']

        self.wcs['A_MIN'] = header['A_MIN3']
        self.wcs['A_MAX'] = header['A_MAX3']

        self.wcs['B_MIN'] = header['B_MIN3']
        self.wcs['B_DELTA'] = header['B_DEL3']
    elif(channel == '4'):
        self.wcs['L_MIN'] = header['L_MIN4']
        self.wcs['L_MAX'] = header['L_MAX4']

        self.wcs['A_MIN'] = header['A_MIN4']
        self.wcs['A_MAX'] = header['A_MAX4']

        self.wcs['B_MIN'] = header['B_MIN4']
        self.wcs['B_DELTA'] = header['B_DEL4']

#    print(' Header wcs',self.wcs['L_MIN'],self.wcs['L_MAX'],self.wcs['A_MIN'],self.wcs['A_MAX'],self.wcs['B_MIN'],self.wcs['B_DELTA'])


def ReadDistortionFile(self, channel, subchannel,):
    """
    Short Summary
    -------------
    Read the distortion files
    Extension 1: Slice Number
    Extension 4 Alpha coefficients 1 or 3
    Extension 5 Lambda coefficients 1 or 3

    Extension 6 Alpha coefficients 2 or 4
    Extension 7 Lambda coefficients 2 or 4

    Extension 12 v2,v3 coefficients

    Parameters
    ----------
    channel: Channel working with
    subchannel: subchannel working with

    Returns
    -------
    returns self.wcs holding a mapping information


    """

    if(channel == '1' or channel == '2'):
        if(subchannel == 'LONG'):
            self.wcs['Filename'] = 'MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_5B.02.00.fits'
        elif (subchannel == 'MEDIUM'):
            self.wcs['Filename'] = 'MIRI_FM_MIRIFUSHORT_12MEDIUM_DISTORTION_5B.02.00.fits'
        elif (subchannel == 'SHORT'):
            self.wcs['Filename'] = 'MIRI_FM_MIRIFUSHORT_12SHORT_DISTORTION_5B.02.00.fits'
    elif(channel == '3' or channel == '4'):
        if(subchannel == 'LONG'):
            self.wcs['Filename'] = 'MIRI_FM_MIRIFULONG_34LONG_DISTORTION_5B.02.00.fits'
        elif(subchannel == 'MEDIUM'):
            self.wcs['Filename'] = 'MIRI_FM_MIRIFULONG_34MEDIUM_DISTORTION_5B.02.00.fits'
        elif(subchannel == 'SHORT'):
            self.wcs['Filename'] = 'MIRI_FM_MIRIFULONG_34SHORT_DISTORTION_5B.02.00.fits'


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    f = fits.open(self.wcs['Filename'])

    print('Reading D2C file', self.wcs['Filename'])

    header = f[0].header
    channel_type = f[0].header['CHANNEL']
    band = f[0].header['BAND']
    detector = f[0].header['DETECTOR']
    ch1 = 'CH{0}'.format(channel_type[0])
    ch2 = 'CH{0}'.format(channel_type[1])
    slices = f[1].data
    alpha1 = f[('Alpha-' + ch1, 1)].data
    lam1 = f[('Lambda-' + ch1, 1)].data
    alpha2 = f[('Alpha-' + ch2, 1)].data
    lam2 = f[('Lambda-' + ch2, 1)].data
    ab_v23 = f[('al,be->V2/V3', 1)].data.copy()
    v23_ab = f[('V2/V3->al,be', 1)].data.copy()
    b0_ch1 = f[0].header['B_ZERO' + ch1[2]]
    bdel_ch1 = f[0].header['B_DEL' + ch1[2]]
    b0_ch2 = f[0].header['B_ZERO' + ch2[2]]
    bdel_ch2 = f[0].header['B_DEL' + ch2[2]]

    ich = 0


    #coeff_names  = build_coeff_names(alpha1.names)
    xa1, acoeff1 = pullout_poly_models(alpha1)
    xl1, lcoeff1 = pullout_poly_models(lam1)

    xa2, acoeff2 = pullout_poly_models(alpha2)
    xl2, lcoeff2 = pullout_poly_models(lam2)


# hard code these for now - this code is just to check assign_wcs not to be used in pipeline
    self.wcs['NUM_COEFF'] = 5
    self.wcs['NUM_V2V3_COEFF'] = 2

    v2coeff = np.zeros((2, 2))
    v3coeff = np.zeros((2, 2))
    nslices = 1
    if(channel == '1'):
        self.wcs['B_MIN'] = b0_ch1
        self.wcs['B_DELTA'] = bdel_ch1
        self.wcs['ALPHA_COEFFS'] = acoeff1
        self.wcs['XALPHA'] = xa1
        self.wcs['LAMBDA_COEFFS'] = lcoeff1
        self.wcs['XLAMBDA'] = xl1

        nslices = len(xa1)
        v2coeff[0][0] = ab_v23[0][1]
        v2coeff[0][1] = ab_v23[0][2]
        v2coeff[1][0] = ab_v23[0][3]
        v2coeff[1][1] = ab_v23[0][4]

        v3coeff[0][0] = ab_v23[1][1]
        v3coeff[0][1] = ab_v23[1][2]
        v3coeff[1][0] = ab_v23[1][3]
        v3coeff[1][1] = ab_v23[1][4]

        ich = 1
    elif(channel == '2'):
        self.wcs['B_MIN'] = b0_ch2
        self.wcs['B_DELTA'] = bdel_ch2
        self.wcs['XALPHA'] = xa2
        self.wcs['ALPHA_COEFFS'] = acoeff2
        self.wcs['XLAMBDA'] = xl2
        self.wcs['LAMBDA_COEFFS'] = lcoeff2

        nslices = len(xa2)
        v2coeff[0][0] = ab_v23[2][1]
        v2coeff[0][1] = ab_v23[2][2]
        v2coeff[1][0] = ab_v23[2][3]
        v2coeff[1][1] = ab_v23[2][4]

        v3coeff[0][0] = ab_v23[3][1]
        v3coeff[0][1] = ab_v23[3][2]
        v3coeff[1][0] = ab_v23[3][3]
        v3coeff[1][1] = ab_v23[3][4]
        ich = 2
    elif(channel == '3'):
        self.wcs['B_MIN'] = b0_ch1
        self.wcs['B_DELTA'] = bdel_ch1

        self.wcs['ALPHA_COEFFS'] = acoeff1
        self.wcs['XALPHA'] = xa1
        self.wcs['LAMBDA_COEFFS'] = lcoeff1
        self.wcs['XLAMBDA'] = xl1

        nslices = len(xa1)
        v2coeff[0][0] = ab_v23[0][1]
        v2coeff[0][1] = ab_v23[0][2]
        v2coeff[1][0] = ab_v23[0][3]
        v2coeff[1][1] = ab_v23[0][4]

        v3coeff[0][0] = ab_v23[1][1]
        v3coeff[0][1] = ab_v23[1][2]
        v3coeff[1][0] = ab_v23[1][3]
        v3coeff[1][1] = ab_v23[1][4]
        ich = 3
    elif(channel == '4'):
        self.wcs['B_MIN'] = b0_ch2
        self.wcs['B_DELTA'] = bdel_ch2

        self.wcs['XALPHA'] = xa2
        self.wcs['ALPHA_COEFFS'] = acoeff2
        self.wcs['XLAMBDA'] = xl2
        self.wcs['LAMBDA_COEFFS'] = lcoeff2

        nslices = len(xa2)
#        nslices = len(alpha1)

        v2coeff[0][0] = ab_v23[2][1]
        v2coeff[0][1] = ab_v23[2][2]
        v2coeff[1][0] = ab_v23[2][3]
        v2coeff[1][1] = ab_v23[2][4]

        v3coeff[0][0] = ab_v23[3][1]
        v3coeff[0][1] = ab_v23[3][2]
        v3coeff[1][0] = ab_v23[3][3]
        v3coeff[1][1] = ab_v23[3][4]
        ich = 4
#    print(' Header wcs',self.wcs['B_MIN'],self.wcs['B_DELTA'])

    self.wcs['V2_COEFFS'] = v2coeff
    self.wcs['V3_COEFFS'] = v3coeff

    alpha_min = []
    alpha_max = []
    beta_min = []
    beta_max = []
    lambda_min = []
    lambda_max = []

    v2_min = []
    v2_max = []
    v3_min = []
    v3_max = []


    print('number of slices determining min and max ', nslices)
    for i in range(nslices):
        #print('on slice',i)
        sliceno = i + 1 + ich * 100
        x, y = find_slices(slices, sliceno) # returns x,y starting at 1
        #print('x y',x[0:10],y[0:10])

        alpha, beta, lam = xy2abl(self, i, x, y) # x,y starts at 1,1
        #print('alpha,beta',alpha[0:10],beta[0:10],lam[0:10])

        alpha_min.append(min(alpha));
        alpha_max.append(max(alpha));

        beta_min.append(min(beta));
        beta_max.append(max(beta));

        lambda_min.append(min(lam));
        lambda_max.append(max(lam));

        if(self.coord_system == 'v2-v3'):

            xan, yan = ab2xyan(self, alpha, beta)
            v2, v3 = xyan2v23(self, xan, yan)

            v2 = v2 * 60.0
            v3 = v3 * 60.0

            v2_min.append(min(v2));
            v2_max.append(max(v2));

            v3_min.append(min(v3));
            v3_max.append(max(v3));


    lmin = min(lambda_min)
    lmax = max(lambda_max)
    if(self.coord_system == 'alpha-beta'):
        amin_final = min(alpha_min)
        amax_final = max(alpha_max)

        bmin_final = min(beta_min)
        bmax_final = max(beta_max)

    elif(self.coord_system == 'v2-v3'):
        amin_final = min(v2_min)
        amax_final = max(v2_max)

        bmin_final = min(v3_min)
        bmax_final = max(v3_max)


    print('min and max alpha, beta, lambda results ', amin_final, amax_final, bmin_final, bmax_final, lmin, lmax)

    return amin_final, amax_final, bmin_final, bmax_final, lmin, lmax


def build_coeff_names(names):
    names = names[1:]
    names = [name.replace('VAR2(', "c") for name in names]
    names = [name.replace(')', "") for name in names]
    names = [name.replace(',', "_") for name in names]
    return names


def pullout_poly_models(data):
    """
    Create a 2D polynomial model for the transformation
    detector --> local MIRI frame
    Works for alpha and lambda coordinates.
    """
    nslices = len(data)

    xa = []
    acoeffs = []
    for i in range(nslices):
        al = data[i]
        xa.append(al[0])
        acoeffs.append(al[1:])

    return xa, acoeffs


def find_slices(slices, sliceno):
    index = np.where(slices == sliceno)
    #y = index[0]+1
    #x = index[1]+1

    y = index[0]
    x = index[1]
    return x, y


def xy2abl(self, slice_no, x, y):


    xa_slice = self.wcs['XALPHA'][slice_no]
    acoeff_slice = self.wcs['ALPHA_COEFFS'][slice_no]

    xl_slice = self.wcs['XLAMBDA'][slice_no]
    lcoeff_slice = self.wcs['LAMBDA_COEFFS'][slice_no]

    beta_zero = self.wcs['B_MIN']
    beta_delta = self.wcs['B_DELTA']


#    print('xa ',xa_slice)
#    print('coeff ',acoeff_slice)

#    print('xl ',xl_slice)
#    print('l coeff ',lcoeff_slice)



    b = beta_zero + (slice_no) * beta_delta

    xa = x - xa_slice
    xl = x - xl_slice


#________________________________________________________________________________

#    if (isinstance(x,np.ndarray)):
#        alpha = []
#        beta = []
#        lam = []
#        for k in range (len(x)):
#            a = 0;
#            l = 0;
#            beta.append(b)
#            h = 0
#            for i in range (self.wcs['NUM_COEFF']):
#                for j in range (self.wcs['NUM_COEFF']):

#                    a = a  + acoeff_slice[h] * pow(xa[k],j) * pow(y[k],i)
#                    l = l  + lcoeff_slice[h] * pow(xl[k],j) * pow(y[k],i)
#                    h = h + 1

#            alpha.append(a)
#            lam.append(l)
#________________________________________________________________________________
#    else:
#        a = 0;
#        l = 0;
#        beta=b
#        h = 0
#        for i in range (self.wcs['NUM_COEFF']):
#            for j in range (self.wcs['NUM_COEFF']):

#                a = a  + acoeff_slice[h] * pow(xa,j) * pow(y,i)
#                l = l  + lcoeff_slice[h] * pow(xl,j) * pow(y,i)
#                h = h + 1

#        alpha=a
#        lam=l




    alpha = x * 0;

    lam = alpha
    beta = alpha + b
    h = 0
    for i in range(self.wcs['NUM_COEFF']):
        for j in range(self.wcs['NUM_COEFF']):
            alpha = alpha + acoeff_slice[h] * pow(xa, j) * pow(y, i)
            lam = lam + lcoeff_slice[h] * pow(xl, j) * pow(y, i)
            h = h + 1

    return alpha, beta, lam




def ab2xyan(self, alpha, beta):

    v2coeff = self.wcs['V2_COEFFS']
    v3coeff = self.wcs['V3_COEFFS']


    xan = alpha * 0
    yan = xan


    for i in range(self.wcs['NUM_V2V3_COEFF']):
        for j in range(self.wcs['NUM_V2V3_COEFF']):
            xan = xan + v2coeff[i][j] * pow(alpha, j) * pow(beta, i)
            yan = yan + v3coeff[i][j] * pow(alpha, j) * pow(beta, i)

    return xan, yan


def xyan2v23(self, xan, yan):
    v2 = xan
    v3 = -(yan + 7.8)

    return v2, v3
