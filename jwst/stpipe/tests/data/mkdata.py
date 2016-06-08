#!/usr/bin/env python
import numpy
from astropy.io import fits as pyfits

name = 'science.fits'
hdu1 = pyfits.PrimaryHDU()
data = numpy.random.random((16, 16))
hdu2 = pyfits.ImageHDU(data, name='SCI')
data = numpy.random.random((16, 16))
hdu3 = pyfits.ImageHDU(data, name='SCI')
hdulist = pyfits.HDUList([hdu1, hdu2, hdu3])
hdulist.writeto(name, clobber=True)

name = 'flat.fits'
data = numpy.random.random((16, 16))
hdu1 = pyfits.PrimaryHDU(data)
hdulist = pyfits.HDUList([hdu1])
hdulist.writeto(name, clobber=True)
