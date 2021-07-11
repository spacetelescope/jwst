import numpy as numpy
from astropy import units
import math
from . import Tools

__author__ = "Do Kester"
__year__ = 2020
__license__ = "GPL3"
__version__ = "2.5.3"
__url__ = "https://www.bayesicfitting.nl"
__status__ = "Perpetual Beta"

#  *
#  * This file is part of the BayesicFitting package.
#  *
#  * BayesicFitting is free software: you can redistribute it and/or modify
#  * it under the terms of the GNU Lesser General Public License as
#  * published by the Free Software Foundation, either version 3 of
#  * the License, or ( at your option ) any later version.
#  *
#  * BayesicFitting is distributed in the hope that it will be useful,
#  * but WITHOUT ANY WARRANTY; without even the implied warranty of
#  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  * GNU Lesser General Public License for more details.
#  *
#  * The GPL3 license can be found at <http://www.gnu.org/licenses/>.
#  *
#  * A JAVA version of this code was part of the Herschel Common
#  * Science System (HCSS), also under GPL3.
#  *
#  *    2006 - 2014 Do Kester, SRON (Java code)
#  *    2017 - 2020 Do Kester


class ImageAssistant( object ):
    """
    ImageAssistant contains 2 methods to assist with more dimensional
    fitting.

    1. getIndices Generates indices for data arrays of any dimension.
       To be used as input in the Fitter classes.
    2. resizeData Resizes the data arrays into a 1-dimensional array.
       To be used as data in the Fitter.




    Example
    -------
    >>> ymap = numpy.arange( 6, dtype=float ).reshape( 2, 3 )
    >>> ias = ImageAssistant()
    >>> ky = ias.getIndices( ymap )
    >>> print( ky.shape )
        (6,2)
    >>> print( ky[4,0], ky[4,1], ymap[ ky[4,0], ky[4,1] ] )
        1 0 4
    >>> ias = ImageAssistant( order='F')
    >>> ky = ias.getIndices( ymap )
    >>> print( ky.shape )
        (6,2)
    >>> print( ky[4,0], ky[4,1], ymap[ ky[4,1], ky[4,0] ] )
        0 1 4

    ## Suppose y is a 2-dimensional map of something
    >>> aass = ImageAssistant( )
    >>> input = aass.getIndices( y )
    >>> fitter = Fitter( input, some2dModel )
    >>> pars = fitter.fit( aass.resizeData( y ) )
    >>> yfit = some2dModel.result( input )                  # Double1d
    >>> yfit2d = aass.resizeData( yfit, shape=y.shape )     # Double2d


    Author       Do Kester

    """


    def __init__( self, order='C' ):
        """
        Helper class to construct from an image, the input arrays
        needed for the Fitters.

        Parameters
        ----------
        order : 'C' or 'F'
            set index view according to character
            'C' orders from slow to fast
            'F' orders from fast to slow
        """
        self.order = order
        self.shape = None

    def getIndices( self, ya, order='C' ):
        """
        Generates indices for data arrays of any dimension.

        To be used as input in the Fitter classes.

        Parameters
        ----------
        ya : map
            array of y ( data ) values for which an indexed array
        order : 'C' or 'F'
            set index view according to character

        Returns
        -------
        numpy.array of ints : the indices of the pixels

        """
        self.shape = ya.shape

        size = ya.size
        rank = ya.ndim

        if rank == 1:
            return numpy.arange( size, dtype=int )
        else:
            kdata = numpy.zeros( ( size, rank ), dtype=int )

            dd = 1
            kr = rank - 1
            for k in range( rank ) :
                for i in range( size ) :
                    kk = kr if self.order == 'C' else k
                    kdata[i,kk] = ( i // dd ) % self.shape[kr]
                dd *= self.shape[kr]
                kr -= 1

            return kdata

    def getPositions( self, ymap, order='C', center=True, deproject=None ) :
        """
        Return the (x,y) positions of the pixels in the map.

        Parameters
        ----------
        ya : map
            array of y ( data ) values for which an indexed array
        order : 'C' or 'F'
            set index view according to character
        center : bool
            if True, return the positions of the center of the pixels.
            otherwise the (left,lower) corner
        deproject : callable
            Deprojection method: from projected map to sky position,
            returning (x,y,...) position given the map indices (ix,iy,...)
            Default: returning the indices as floats (+0.5 if center)

        Returns
        -------
        numpy.array of floats : the positions of the pixels
        """
        xdata = numpy.asarray( self.getIndices( ymap, order=order ), dtype=float )
        if center :
            xdata += 0.5
        if deproject is not None :
            xdata = deproject( xdata )
        return xdata


    def getydata( self, ya ):
        """
        Return a copy of ya as a 1 dim array.

        Parameters
        ----------
        ya : array_like
            map to be reshaped
        """
        return numpy.reshape( ya.copy(), ya.size, order=self.order )

    def resizeData( self, res, shape=None ):
        """
        Reshape the data (res) into the same shape as the map (ya)

        Parameters
        ----------
        res : array_like
            result of the fit as a 1-dim array
        shape : tuple of int
            dimensional lengths of the reconstructable map
            default remembered from a call to getIndices
        """
        if shape is None :
            shape = self.shape

        return numpy.reshape( res, shape, order=self.order )



