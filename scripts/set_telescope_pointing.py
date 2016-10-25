#!/usr/bin/env python

# Copyright (C) 2010-2011 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

'''
This script adds absolute pointing information to the FITS files provided
to it on the command line (one or more).

Currently it only uses a constant value for the engineering keywords
since the Engineering Database does not yet contain them.

It assumes the following keywords are present in the file header:

V2_REF (arcseconds)
V3_REF (arcseconds)
VPARITY (+1 or -1)
V3I_YANG (decimal degrees)

The keywords added are:

RA_V1
DEC_V1
PA_V3
CRVAL1
CRVAL2
PC1_1
PC1_2
PC2_1
PC2_2

It does not currently place the new keywords in any particular location
in the header other than what is required by the standard.
'''

from __future__ import print_function, division

import numpy as np
from numpy import cos, sin
import sys
import astropy.io.fits as fits


def m_v_to_siaf(ya, v3, v2, vidlparity):  # This is a 321 rotation
    mat = np.array([[cos(v3)*cos(v2),
                    cos(v3)*sin(v2),
                    sin(v3)],
                   [-cos(ya)*sin(v2)+sin(ya)*sin(v3)*cos(v2),
                    cos(ya)*cos(v2)+sin(ya)*sin(v3)*sin(v2),
                    -sin(ya)*cos(v3)],
                   [-sin(ya)*sin(v2)-cos(ya)*sin(v3)*cos(v2),
                    sin(ya)*cos(v2)-cos(ya)*sin(v3)*sin(v2),
                    cos(ya)*cos(v3)]])
    pmat = np.array([[0., vidlparity, 0.],
                     [0., 0., 1.],
                     [1., 0., 0.]])
    return np.dot(pmat, mat)


def vector_to_ra_dec(v):
    """Returns tuple of spherical angles from unit direction Vector """
    ra = np.arctan2(v[1], v[0])
    dec = np.arcsin(v[2])
    if ra < 0.:
        ra += 2. * np.pi
    return(ra, dec)

R2D = 180./np.pi
D2R = np.pi/180.
A2R = D2R/3600.
R2A = 3600.*R2D

# Define the FGS1 to SI-FOV DCM,  Transpose of DCM in SE-20 section 5.8.4.2.

m1 = np.array(
    [[0.9999994955442, 0.0000000000000, 0.0010044457459],
     [0.0000011174826, 0.9999993811310, -0.0011125359826],
     [-0.0010044451243, 0.0011125365439, 0.9999988766756]])
m2 = np.array(
    [[0, 0, 1],
     [1, 0, 0],
     [0, 1, 0]])

m_fgs1_to_sifov = np.dot(m2, m1)
m_fgs1_to_sifovT = m_fgs1_to_sifov.transpose()

# Define the SI-FOV to V-frame DCM,  From Luis' IOC
m_sifov_to_v = np.array(
    [[0.99999742598, 0., 0.00226892608],
     [0., 1., 0.],
     [-0.00226892608, 0., 0.99999742598]])


def calc_wcs(v2ref, v3ref, v3idlyang, vidlparity,
             q, j2fgs_matrix, fsmcorr):
    '''
    v2ref (arcsec), v3ref (arcsec), v3idlyang (deg), vidlparity (+1 or -1),
    are the relevant siaf parameters. The assumed units are shown in
    parentheses.

    It is assumed that the siaf ref position is the corresponding WCS
    reference position.

    Parameter q is the SA_ZATTEST<n> engineering parameters where
    n ranges from 1 to 4.

    Parameter j2fgs_matrix is the transformation matrix specified by
    engineering parameters SA_ZRFGS2J<n><m> where both n and m range
    from 1 to 3. This is to be provided as a 1d list using this order:
    11, 21, 31, 12, 22, 32, 13, 23, 33

    Parameter fsmcorr are two values provided as a list consisting of:
    [SA_ZADUCMDX, SA_ZADUCMDY]

    This routine returns two tuples:
    The first is of (CRVAL1, CRVAL2, PA_y_axis)
    The second is of (V1ra, V1dec, V3pa)
    All angles are in decimal degrees.
    '''

    q1, q2, q3, q4 = q
    m_eci2j = np.array(
        [[1. - 2.*q2*q2 - 2.*q3*q3,
            2.*(q1*q2 + q3*q4),
            2.*(q3*q1 - q2*q4)],
         [2.*(q1*q2 - q3*q4),
            1. - 2.*q3*q3 - 2.*q1*q1,
            2.*(q2*q3 + q1*q4)],
         [2.*(q3*q1 + q2*q4),
            2.*(q2*q3 - q1*q4),
            1. - 2.*q1*q1 - 2.*q2*q2]])

    mj2fgs1 = np.array(j2fgs_matrix).reshape((3, 3)).transpose()

    m_sifov_fsm_delta = np.array(
        [[1., fsmcorr[0]/22.01, fsmcorr[1]/21.68],
         [-fsmcorr[0]/22.01, 1., 0.],
         [-fsmcorr[1]/21.68, 0., 1.]])

    mpartial = np.dot(m_sifov_to_v,
                      np.dot(m_sifov_fsm_delta,
                             np.dot(m_fgs1_to_sifov, mj2fgs1)))

    m_eci2v = np.dot(mpartial, m_eci2j)

    v1pt = m_eci2v[0]
    xeci_ra, xeci_dec = vector_to_ra_dec(v1pt)

    # V1 is given by the first row.
    v1_ra, v1_dec = vector_to_ra_dec(m_eci2v[0])
    # V3 is given by the third row
    v3_ra, v3_dec = vector_to_ra_dec(m_eci2v[2])

    # The V3PA @ V1 is given by
    y = cos(v3_dec) * sin(v3_ra-v1_ra)
    x = sin(v3_dec) * cos(v1_dec) - \
        cos(v3_dec) * sin(v1_dec) * cos((v3_ra - v1_ra))
    V3PA = np.arctan2(y, x)

    m_eci2siaf = np.dot(m_v_to_siaf(v3idlyang * D2R,
                                    v3ref * A2R,
                                    v2ref * A2R,
                                    vidlparity), m_eci2v)
    siaf_x = 0. * A2R
    siaf_y = 0. * A2R
    refpos = np.array(
               [siaf_x,
                siaf_y,
                np.sqrt(1.-siaf_x * siaf_x - siaf_y * siaf_y)])
    msky = np.dot(m_eci2siaf.transpose(), refpos)
    vaper_ra, vaper_dec = vector_to_ra_dec(msky)

    vysiaf = np.array([0., 1., 0.])

    myeci = np.dot(m_eci2siaf.transpose(), vysiaf)
    # The Y axis of the aperture is given by
    vy_ra, vy_dec = vector_to_ra_dec(myeci)
    # The VyPA @ xref,yref is given by
    y = cos(vy_dec) * sin(vy_ra-vaper_ra)
    x = sin(vy_dec) * cos(vaper_dec) - \
        cos(vy_dec) * sin(vaper_dec) * cos((vy_ra - vaper_ra))
    vypa = np.arctan2(y, x)
    wcsinfo = (vaper_ra*R2D, vaper_dec*R2D, vypa*R2D)
    vinfo = (v1_ra*R2D, v1_dec*R2D, V3PA*R2D)
    return wcsinfo, vinfo


def add_wcs(filename):
    '''
    Given the name of a valid partially populated level 1b JWST file,
    determine the simple WCS parameters from the SIAF keywords in that
    file and the engineering parameters that contain information about
    the telescope pointing.

    It presumes all the accessed keywords are present (see first block).
    '''
    hdulist = fits.open(filename, 'update')
    pheader = hdulist[0].header
    obsstart = float(pheader['EXPSTART'])
    obsend = float(pheader['EXPEND'])
    q, j2fgs_matrix, fsmcorr = get_pointing(obsstart, obsend)
    v2ref = float(pheader['V2_REF'])
    v3ref = float(pheader['V3_REF'])
    v3idlyang = float(pheader['V3I_YANG'])
    vparity = int(pheader['VPARITY'])

    # compute relevant WCS information
    wcsinfo, vinfo = calc_wcs(v2ref, v3ref, v3idlyang, vparity,
                              q, j2fgs_matrix, fsmcorr)
    crval1, crval2, pa_aper_deg = wcsinfo
    v1_ra_deg, v1_dec_deg, v3_pa_deg = vinfo
    # update header
    pheader['RA_V1'] = v1_ra_deg
    pheader['DEC_V1'] = v1_dec_deg
    pheader['PA_V3'] = v3_pa_deg
    pheader['CRVAL1'] = crval1
    pheader['CRVAL2'] = crval2
    pheader['PC1_1'] = np.cos(pa_aper_deg * D2R)
    pheader['PC1_2'] = -np.sin(pa_aper_deg * D2R)
    pheader['PC2_1'] = np.sin(pa_aper_deg * D2R)
    pheader['PC2_2'] = np.cos(pa_aper_deg * D2R)
    hdulist.flush()
    hdulist.close()


def get_pointing(obstart, obsend):
    '''
    For the time being this simply returns the same damn values regardless of
    the input time (awaiting the time that these parameters are actually
    in the engineering database)
    '''
    # The following values of q correspond to the engineering keyword values:
    # SA_ZATEST1, SA_ZATEST2, SA_ZATEST3, SA_ZATEST4
    q = np.array([-0.36915286,  0.33763282,  0.05758533,  0.86395264])
    # The following values of j2fgs_matrix correspond to the engineering
    #  keyword values of:
    # SA_ZRFGS2K11 SA_ZRFGS2K21 SA_ZRFGS2K31
    # SA_ZRFGS2K21 SA_ZRFGS2K22 SA_ZRFGS2K32
    # SA_ZRFGS2K31 SA_ZRFGS2K32 SA_ZRFGS2K33
    j2fgs_matrix = np.array(
        [-1.00444000e-03,  3.38145836e-03, 9.99993778e-01,
         9.99999496e-01, -3.90000000e-14, 1.00444575e-03,
         3.39649146e-06,  9.99994283e-01, -3.38145665e-03])
    # The following values of fsmcorr correspond to the engineering keywords:
    # SA_ZADUCMDX, SA_ZADUCMDY
    fsmcorr = np.array([0., 0.])
    return q, j2fgs_matrix, fsmcorr

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        raise ValueError('missing filename argument(s)')
    for filename in sys.argv[1:]:
        add_wcs(filename)
