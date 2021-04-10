import numpy as numpy

from .Kernel import Kernel

__author__ = "Do Kester"
__year__ = 2017
__license__ = "GPL3"
__version__ = "0.9"
__maintainer__ = "Do"
__status__ = "Development"

#  *
#  *    This file is part of the BayesicFitting package.
#  *
#  *    BayesicFitting is free software: you can redistribute it and/or modify
#  *    it under the terms of the GNU Lesser General Public License as
#  *    published by the Free Software Foundation, either version 3 of
#  *    the License, or ( at your option ) any later version.
#  *
#  *    BayesicFitting is distributed in the hope that it will be useful,
#  *    but WITHOUT ANY WARRANTY; without even the implied warranty of
#  *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  *    GNU Lesser General Public License for more details.
#  *
#  *    The GPL3 license can be found at <http://www.gnu.org/licenses/>.
#  *
#  *    2017 - 2018 Do Kester

class Tophat( Kernel ):
    """
    Tophat is a Kernel function which is 1.0 between [-0.5,0.5]; it is 0 elsewhere.

    Attributes
    ----------
    nconv : int
        successive autoconvolutions of the tophat. max=6.

    Thanks to Romke Bontekoe and Mathematica for providing the analytic expressions.

    """
    FWHM = [0.5, 0.5, 0.6339745962155612, 0.7223517244643762, 0.7971951696335494,
            0.8660920722545018, 0.9301994777857857]
    RANGE = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
    INTEGRAL = [1.0, 1.0, 4/3, 3/2, 192/115, 20/11, 11520/5887]


    def __init__( self, nconv=0 ) :
        """
        Constructor.

        Integral, fwhm and range are dependent on the number of convolutions.

        Parameters
        ----------
        nconv : int
            number of auto-convolutions

        """


        super( Tophat, self ).__init__( integral=self.INTEGRAL[nconv],
                        fwhm=2*self.FWHM[nconv], range=self.RANGE[nconv] )
        self.nconv = nconv

        if nconv == 0 :
            self.conv = Conv0()
        elif nconv == 1 :
            self.conv = Conv1()
        elif nconv == 2 :
            self.conv = Conv2()
        elif nconv == 3 :
            self.conv = Conv3()
        elif nconv == 4 :
            self.conv = Conv4()
        elif nconv == 5 :
            self.conv = Conv5()
        elif nconv == 6 :
            self.conv = Conv6()
        else :
            raise ValueError( "Cannot handle more than 6 convolutions." )

    def result( self, x ):
        return self.conv.result( x )

    def resultsq( self, xsq ):
        return self.result( xsq )                   #  the same

    def partial( self, x ):
        return self.conv.partial( x )

    def isBound( self ):
        return True

    def name( self ):
        return str( "Tophat %d convolved" % self.nconv )

class Conv0( object ) :
    def result( self, x ) :
        ax = numpy.abs( x )
        return numpy.where( ax < 0.5, 1.0, numpy.where( ax == 0.5, 0.5, 0.0 ) )

    def partial( self, x ):
        return numpy.zeros_like( x )

class Conv1( object ) :
    def result( self, x ) :
        ax = numpy.abs( x )
        return numpy.where( ax < 1.0, 1.0 - ax, 0.0 )

    def partial( self, x ):
        ax = numpy.abs( x )
        return numpy.where( ax < 1.0, - numpy.sign( x ), 0 )

class Conv2( object ) :

    def result( self, x ) :
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 1.5 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        x1 = ax[q]
        x2 = x1 * x1

        b = numpy.where( x1 < 0.5 )
        res[q[b]] = ( 3 - 4 * x2[b] ) / 3,
        b = numpy.where( x1 >= 0.5 )
        res[q[b]] = ( 9 - 12 * x1[b] + 4 * x2[b] ) / 6

        return res

    def partial( self, x ):
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 1.5 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        s0 = numpy.sign( x[q] )
        x1 = ax[q]

        b = numpy.where( x1 < 0.5 )
        res[q[b]] = s0[b] * ( -8 * x1[b] ) / 3
        b = numpy.where( x1 >= 0.5 )
        res[q[b]] = s0[b] * ( -12 + 8 * x1[b] ) / 6

        return res

class Conv3( object ) :
    def result( self, x ) :
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 2.0 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        x1 = ax[q]
        x2 = x1 * x1
        x3 = x2 * x1

        b = numpy.where( x1 < 1.0 )
        res[q[b]] = ( 4 - 6 * x2[b] + 3 * x3[b] ) / 4
        b = numpy.where( x1 >= 1.0 )
        res[q[b]] = ( 8 - 12 * x1[b] + 6 * x2[b] - x3[b] ) / 4

        return res

    def partial( self, x ):
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 2.0 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        s0 = numpy.sign( x[q] )
        x1 = ax[q]
        x2 = x1 * x1

        b = numpy.where( x1 < 1.0 )
        res[q[b]] = s0[b] * ( - 12 * x1[b] + 9 * x2[b] ) / 4
        b = numpy.where( x1 >= 1.0 )
        res[q[b]] = s0[b] * ( - 12 + 12 * x1[b] - 3 * x2[b] ) / 4

        return res

class Conv4( object ) :
    def result( self, x ) :
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 2.5 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        x1 = ax[q]
        x2 = x1 * x1
        x3 = x2 * x1
        x4 = x3 * x1

        res[q] = ( 625 - 1000 * x1 + 600 * x2 - 160 * x3 + 16 * x4 ) / 230
        b = numpy.where( x1 < 1.5 )
        res[q[b]] = ( 55 + 20 * x1[b] - 120 * x2[b] + 80 * x3[b] - 16 * x4[b] ) * 2 / 115
        b = numpy.where( x1 < 0.5 )
        res[q[b]] = ( 115 - 120 * x2[b] +  48 * x4[b] ) / 115

        return res

    def partial( self, x ):
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 2.5 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        s0 = numpy.sign( x[q] )
        x1 = ax[q]
        x2 = x1 * x1
        x3 = x2 * x1

        res[q] = s0 * ( - 1000 + 1200 * x1 - 480 * x2 + 64 * x3 ) / 230
        b = numpy.where( x1 < 1.5 )
        res[q[b]] = s0[b] * ( + 20  - 240 * x1[b] + 240 * x2[b] - 64 * x3[b] ) * 2 / 115
        b = numpy.where( x1 < 0.5 )
        res[q[b]] = s0[b] * ( - 240 * x1[b] +  192 * x3[b] ) / 115

        return res

class Conv5( object ) :
    def result( self, x ) :
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 3.0 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        x1 = ax[q]
        x2 = x1 * x1
        x3 = x2 * x1
        x4 = x3 * x1
        x5 = x4 * x1

        res[q] = ( 243 - 405 * x1 + 270 * x2 - 90 * x3 + 15 * x4 - x5 ) / 66
        b = numpy.where( x1 < 2.0 )
        res[q[b]] = ( 51 + 75 * x1[b] - 210 * x2[b] + 150 * x3[b] - 45 * x4[b] + 5 * x5[b] ) / 66
        b = numpy.where( x1 < 1.0 )
        res[q[b]] = ( 33 - 30 * x2[b] + 15 * x4[b] - 5 * x5[b] ) / 33

        return res

    def partial( self, x ):
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 3.0 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        s0 = numpy.sign( x[q] )
        x1 = ax[q]
        x2 = x1 * x1
        x3 = x2 * x1
        x4 = x3 * x1

        res[q] = s0 * ( - 405 + 540 * x1 - 270 * x2 + 60 * x3 - 5 * x4 ) / 66
        b = numpy.where( x1 < 2.0 )
        res[q[b]] = s0[b] * ( + 75 - 420 * x1[b] + 450 * x2[b] - 180 * x3[b] + 25 * x4[b] ) / 66
        b = numpy.where( x1 < 1.0 )
        res[q[b]] = s0[b] * ( - 60 * x1[b] + 60 * x3[b] - 25 * x4[b] ) / 33

        return res

class Conv6( object ) :
    def result( self, x ) :
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 3.5 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        x1 = ax[q]
        x2 = x1 * x1
        x3 = x2 * x1
        x4 = x3 * x1
        x5 = x4 * x1
        x6 = x5 * x1

        b = numpy.where( x1 >= 2.5 )
        res[q[b]] = ( 117649 - 201684 * x1[b] + 144060 * x2[b] - 54880 * x3[b] + 11760 * x4[b] -
                    1344 * x5[b] + 64 * x6[b] ) / 23548
        b = numpy.where( numpy.logical_and( x1 < 2.5, x1 >= 1.5 ) )
        res[q[b]] = ( 4137 + 30408 * x1[b] - 59220 * x2[b] + 42560 * x3[b] - 15120 * x4[b] +
                    2688 * x5[b] - 192 * x6[b] ) / 11774
        b = numpy.where( numpy.logical_and( x1 < 1.5, x1 >= 0.5 ) )
        res[q[b]] = ( 23583 - 420 * x1[b] - 16380 * x2[b] - 5600 * x3[b] + 15120 * x4[b] -
                    6720 * x5[b] + 960 * x6[b] ) / 23548
        b = numpy.where( x1 < 0.5 )
        res[q[b]] = ( 5887 - 4620 * x2[b] + 1680 * x4[b] - 320 * x6[b] ) / 5887

        return res

    def partial( self, x ):
        x = numpy.array( x, copy=False, ndmin=1 )
        ax = numpy.abs( x )
        q = numpy.where( ax < 3.5 )[0]
        res = numpy.zeros_like( x )
        if len( q ) == 0 :
            return res

        s0 = numpy.sign( x[q] )
        x1 = ax[q]
        x2 = x1 * x1
        x3 = x2 * x1
        x4 = x3 * x1
        x5 = x4 * x1

        b = numpy.where( x1 >= 2.5 )
        res[q[b]] = s0[b] * ( - 201684 + 2*144060 * x1[b] - 3*54880 * x2[b] + 4*11760 * x3[b] -
                    5*1344 * x4[b] + 6*64 * x5[b] ) / 23548
        b = numpy.where( numpy.logical_and( x1 < 2.5, x1 >= 1.5 ) )
        res[q[b]] = s0[b] * ( + 30408 - 2*59220 * x1[b] + 3*42560 * x2[b] - 4*15120 * x3[b] +
                    5*2688 * x4[b] - 6*192 * x5[b] ) / 11774
        b = numpy.where( numpy.logical_and( x1 < 1.5, x1 >= 0.5 ) )
        res[q[b]] = s0[b] * ( - 420 - 2*16380 * x1[b] - 3*5600 * x2[b] + 4*15120 * x3[b] -
                    5*6720 * x4[b] + 6*960 * x5[b] ) / 23548
        b = numpy.where( x1 < 0.5 )
        res[q[b]] = s0[b] * ( - 2*4620 * x1[b] + 4*1680 * x3[b] - 6*320 * x5[b] ) / 5887

        return res





