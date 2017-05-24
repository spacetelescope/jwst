# Cube Class
# Spaxel Class

import sys
import numpy as np
import math
import logging
from .. import datamodels
from . import coord

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)



class CubeInfo(object):
# Array of classes, CubeType[i]  defines type of Cube being created and the list
# of files are are used in making the cube.

    def __init__(self, instrument,detector, parameter1, parameter2, output_name,coord_system):

        self.channel = []
        self.subchannel = []

        self.file = []

        self.transform_v23toab = []
        self.transform_worldtov23 = []

        self.filter = []
        self.grating = []

        self.output_name = ''
        self.detector = detector
        self.instrument = instrument
        self.output_name = output_name
        self.coord_system = coord_system
        if(instrument == 'MIRI'):
            self.channel = parameter1
            self.subchannel = parameter2
            
        elif(instrument == 'NIRSPEC'):
            self.filter = parameter1
            self.grating = parameter2

#_______________________________________________________________________
    def SetScale(self, a_scale, b_scale, wscale):
        self.Cdelt1 = a_scale
        self.Cdelt2 = b_scale
        self.Cdelt3 = wscale
#_______________________________________________________________________
# cube in alpha-beta space (single exposure cube - small FOV assume rectangular coord system
    def SetGeometryAB(self, footprint):
        self.a_min, self.a_max, self.b_min, self.b_max, self.lambda_min, self.lambda_max = footprint

        #set up the a (x) coordinates of the cube
        range_a = self.a_max - self.a_min
        self.naxis1 = int(math.ceil(range_a / self.Cdelt1))

        # adjust min and max based on integer value of naxis1
        a_center = (self.a_max + self.a_min) / 2.0
        self.a_min = a_center - (self.naxis1 / 2.0) * self.Cdelt1
        self.a_max = a_center + (self.naxis1 / 2.0) * self.Cdelt1

        self.xcoord = np.zeros(self.naxis1)
        self.Crval1 = self.a_min
        self.Crpix1 = 0.5
        xstart = self.a_min + self.Cdelt1 / 2.0
        for i in range(self.naxis1):
            self.xcoord[i] = xstart
            xstart = xstart + self.Cdelt1

#_______________________________________________________________________
        #set up the lambda (z) coordinate of the cube

        range_lambda = self.lambda_max - self.lambda_min
        self.naxis3 = int(math.ceil(range_lambda / self.Cdelt3))

         # adjust max based on integer value of naxis3
        lambda_center = (self.lambda_max + self.lambda_min) / 2.0

        self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.Cdelt3
        self.lambda_max = lambda_center + (self.naxis3 / 2.0) * self.Cdelt3

        self.lambda_max = self.lambda_min + (self.naxis3) * self.Cdelt3

        self.zcoord = np.zeros(self.naxis3)
        self.Crval3 = self.lambda_min
        self.Crpix3 = 1.0
        zstart = self.lambda_min + self.Cdelt3 / 2.0

        for i in range(self.naxis3):
            self.zcoord[i] = zstart
            zstart = zstart + self.Cdelt3
#_______________________________________________________________________

        range_b = self.b_max - self.b_min

        self.naxis2 = int(math.ceil(range_b / self.Cdelt2))
        b_center = (self.b_max + self.b_min) / 2.0

        # adjust min and max based on integer value of naxis2
        self.b_max = b_center + (self.naxis2 / 2.0) * self.Cdelt2
        self.b_min = b_center - (self.naxis2 / 2.0) * self.Cdelt2


        self.ycoord = np.zeros(self.naxis2)
        self.Crval2 = self.b_min
        self.Crpix2 = 0.5
        ystart = self.b_min + self.Cdelt2 / 2.0
        for i in range(self.naxis2):
            self.ycoord[i] = ystart
            ystart = ystart + self.Cdelt2

#_______________________________________________________________________
# Footprint values are RA,DEC values on the sky
# Values are given in degrees 

    def SetGeometry(self, footprint):

        deg2rad = math.pi/180.0
        ra_min, ra_max, dec_min, dec_max,lambda_min, lambda_max = footprint # in degrees
        dec_ave = (dec_min + dec_max)/2.0

        # actually this is hard due to converenge of hour angle
        # improve determining ra_ave in the future - do not just average (BAD) 
        ra_ave = ((ra_min + ra_max)/2.0 )#* math.cos(dec_ave*deg2rad) 
#        range_ra = (ra_max - ra_min) * 3600.0 * math.cos(dec_ave*deg2rad)
#        range_dec = (dec_max - dec_min) * 3600.0

        self.Crval1 = ra_ave 
        self.Crval2 = dec_ave
        xi_center,eta_center = coord.radec2std(self.Crval1, self.Crval2,ra_ave,dec_ave)
        
        xi_min,eta_min = coord.radec2std(self.Crval1, self.Crval2,ra_min,dec_min)
        xi_max,eta_max = coord.radec2std(self.Crval1, self.Crval2,ra_max,dec_max)
#________________________________________________________________________________
        # find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
        # to find location of center abs of min values is how many pixels 

        n1a = int(math.ceil(math.fabs(xi_min) / self.Cdelt1)) 
        n2a = int(math.ceil(math.fabs(eta_min) / self.Cdelt2)) 

        n1b = int(math.ceil(math.fabs(xi_max) / self.Cdelt1)) 
        n2b = int(math.ceil(math.fabs(eta_max) / self.Cdelt2)) 

        xi_min = 0.0 - (n1a * self.Cdelt1) - self.Cdelt1/2.0
        xi_max = (n1b * self.Cdelt1) + self.Cdelt1/2.0

        eta_min = 0.0 - (n2a * self.Cdelt2) - self.Cdelt2/2.0
        eta_max = (n2b * self.Cdelt2) + self.Cdelt2/2.0
        
        self.Crpix1 = n1a
        self.Crpix2 = n2a

        self.naxis1 = n1a + n1b
        self.naxis2 = n2a + n2b

        self.a_min  = xi_min
        self.a_max = xi_max
        self.b_min = eta_min
        self.b_max = eta_max

# center of spaxels 
        self.xcoord = np.zeros(self.naxis1)
        xstart = xi_min + self.Cdelt1 / 2.0
        for i in range(self.naxis1):
            self.xcoord[i] = xstart
            xstart = xstart + self.Cdelt1
#            print('xcoord',self.xcoord[i],i)

        self.ycoord = np.zeros(self.naxis2)
        ystart = eta_min + self.Cdelt2 / 2.0

        for i in range(self.naxis2):
            self.ycoord[i] = ystart
            ystart = ystart + self.Cdelt2
#            print('ycoord',self.ycoord[i],i)

#        print('ycoord xcoord shape',self.ycoord.shape,self.xcoord.shape)
        
#_______________________________________________________________________
        
#        ystart = self.ycoord[0]
#        yend = self.ycoord[0] + self.Cdelt2*(self.naxis2)
        
#        xstart = self.xcoord[0]
#        xend = self.xcoord[0] + self.Cdelt1*(self.naxis1)

#        yy,xx = np.mgrid[ystart:yend:self.Cdelt2,
#                         xstart:xend:self.Cdelt1]

        ygrid = np.zeros(self.naxis2*self.naxis1)
        xgrid = np.zeros(self.naxis2*self.naxis1)

        k = 0 
        ystart = self.ycoord[0]
        for i in range(self.naxis2):
            xstart = self.xcoord[0]
            for j in range(self.naxis1): 
                xgrid[k] = xstart
                ygrid[k] = ystart
                xstart = xstart + self.Cdelt1
                k = k + 1
            ystart = ystart + self.Cdelt2
        

#        print('y start end',ystart,yend)
#        print('x start end',xstart,xend)

#        print('yy shape',yy.shape,self.ycoord.shape)
#        print('xx shape',xx.shape,self.xcoord.shape)

#        self.Ycenters = np.ravel(yy)
#        self.Xcenters = np.ravel(xx)

        self.Xcenters = xgrid
        self.Ycenters = ygrid
#_______________________________________________________________________
        #set up the lambda (z) coordinate of the cube

        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        range_lambda = self.lambda_max - self.lambda_min
        self.naxis3 = int(math.ceil(range_lambda / self.Cdelt3))

         # adjust max based on integer value of naxis3
        lambda_center = (self.lambda_max + self.lambda_min) / 2.0
        self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.Cdelt3
        self.lambda_max = self.lambda_min + (self.naxis3) * self.Cdelt3

        self.zcoord = np.zeros(self.naxis3)
        self.Crval3 = self.lambda_min
        self.Crpix3 = 1.0
        zstart = self.lambda_min + self.Cdelt3 / 2.0

        for i in range(self.naxis3):
            self.zcoord[i] = zstart
            zstart = zstart + self.Cdelt3

#_______________________________________________________________________
    def PrintCubeGeometry(self, instrument):
        log.info('Cube Geometry:')
        blank = '  '
        if (self.coord_system == 'alpha-beta'):
            log.info('axis# Naxis  CRPIX    CRVAL      CDELT(arc sec)  MIN & Max (alpha,beta arc sec)')
        else:
            log.info('axis# Naxis  CRPIX    CRVAL      CDELT(arc sec)  MIN & Max (xi,eta arc sec)')
        log.info('Axis 1 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f', 
                 self.naxis1, self.Crpix1, self.Crval1, self.Cdelt1, self.a_min, self.a_max)
        log.info('Axis 2 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f', 
                 self.naxis2, self.Crpix2, self.Crval2, self.Cdelt2, self.b_min, self.b_max)
        log.info('Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f', 
                 self.naxis3, self.Crpix3, self.Crval3, self.Cdelt3, self.lambda_min, self.lambda_max)

        if(instrument == 'MIRI'):
            # length of channel and subchannel are the same 
            number_bands = len(self.channel)

            for i in range(number_bands):
                this_channel = self.channel[i]
                this_subchannel = self.subchannel[i]
                log.info('Cube covers channel, subchannel: %s %s ', this_channel,this_subchannel)
        elif(instrument == 'NIRSPEC'):
            # number of filters and gratings are the same
            number_bands = len(self.filter)

            for i in range(number_bands):
                this_fwa = self.filter[i]
                this_gwa = self.grating[i]
                log.info('Cube covers grating, filter: %s %s ', this_gwa,this_fwa)

##################################################################################
class Spaxel(object):


    __slots__ = ['flux', 'error','flux_weight','iflux']

    def __init__(self):
        self.flux = 0.0
        self.flux_weight = 0.0
        self.iflux = 0
        self.error = 0

class SpaxelAB(object):

    __slots__ = ['flux', 'error', 'flux_weight','iflux']

    def __init__(self):

        self.flux = 0
        self.error = 0
        self.flux_weight = 0.0
        self.iflux = 0.0
