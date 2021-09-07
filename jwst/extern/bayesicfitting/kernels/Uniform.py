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
#  * A JAVA version of this code was part of the Herschel Common
#  * Science System (HCSS), also under GPL3.
#  *
#  *    2010 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2017 Do Kester

class Uniform( Kernel ):
    """
    Uniform is a Kernel function which is constant between [-1,1].

        K( x ) = 1.0        if |x| < 1
                 0.0        elsewhere

    """
    def __init__( self ) :
        """
        Constructor.

        Using
            integral = 2.0
            fwhm = 2.0
        """
        super( Uniform, self ).__init__( integral=2.0, fwhm=2.0 )

    def result( self, x ):
        ax = numpy.abs( x )
        return numpy.where( ax < 1, 1.0, numpy.where( ax == 1, 0.5, 0.0 ) )

    def resultsq( self, xsq ):
        return self.result( xsq )                   #  the same

    def partial( self, x ):
        return numpy.zeros_like( x )

    def isBound( self ):
        return True

    def name( self ):
        return str( "Uniform: 1 if |x| < 1 else 0" )


