
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import sys
import numpy as np
import math


#_______________________________________________________________________
def V2V32RADEC(ra_ref,dec_ref,roll_ref,v2_ref,v3_ref,v2, v3):

        # v2_ref,v3_ref given in arc mins
        # v2 and v3 are given in arc mins
        # ra_ref,dec_ref, roll_ref given in degrees
        # it is assumed that the V2,V3 coordinates have the effects of dithering included
        
    d2r = math.pi/180.0

    #v2_ref = -8.3942412 # in arcmins
    #v3_ref = -5.3123744 # in arcmins
    #roll_ref = 0.0     # degress
    #ra_ref = 45.0      # degrees
    #dec_ref = 0.0     # degrees

    v2 = v2/60.0      # convert to degrees
    v3 = v3/60.0      # convert to degrees

    v2_ref = v2_ref/60.0 # covert to degrees
    v3_ref = v3_ref/60.0 # convert to degrees
    v3_ref_rad = v3_ref*d2r
    roll_ref_rad = roll_ref * d2r
        
    delta_v2 = (v2 - v2_ref) * math.cos(v3_ref_rad)
    delta_v3 = (v3-v3_ref)
    delta_ra = delta_v2 * math.cos(roll_ref_rad) + delta_v3* math.sin(roll_ref_rad)
    delta_dec = -delta_v2* math.sin(roll_ref_rad) + delta_v3*math.cos(roll_ref_rad)

    ra = ra_ref + delta_ra/math.cos(dec_ref*d2r)
    dec = delta_dec + dec_ref

    #print('v2,v3',v2[0,15:20],v3[0,15:20])
    #print('ra dec',ra[0,15:20],dec[0,15:20])


    return ra,dec
#_______________________________________________________________________

def radec2std(RA0,DEC0,ra,dec):

# Compute the standard coordinates xi,eta from CRVAL1,CRVAL2

# tangent projection, ra dec to xi, eta
  deg2rad = math.pi/180.0
  rad2arcsec = (180.0*3600.0)/math.pi

  crval1 = RA0*deg2rad
  crval2 = DEC0*deg2rad
  radiff = ra*deg2rad - crval1;
  decr = dec*deg2rad;

  h = np.sin(decr) *math.sin(crval2) + np.cos(decr)*math.cos(crval2)*np.cos(radiff);

  xi = np.cos(decr)*np.sin(radiff)/h;
  eta = ( np.sin(decr)*math.cos(crval2) - np.cos(decr)*math.sin(crval2)*np.cos(radiff) )/h;

  #xi = xi/deg2rad
  #eta = eta/deg2rad;
  xi = xi * rad2arcsec
  eta = eta * rad2arcsec


  return xi,eta
