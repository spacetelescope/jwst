import numpy as numpy
import math

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
#  *    2013 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2017 Do Kester

class Parabola( Kernel ):
    """
    Parabola is a Kernel function between [-1,1]; it is 0 elsewhere.

        K( x ) = ( 1-x^2 )      if |x| < 1
                 0              elsewhere

    """
    def __init__( self ) :
        """
        Constructor.

        Using
            integral = 4.0/3.0
            fwhm = 2 sqrt( 0.5 )
        """
        super( Parabola, self ).__init__( integral=4.0/3.0, fwhm=2 * math.sqrt(0.5) )

    def result( self, x ):
        return self.resultsq( x * x )

    def resultsq( self, xsq ):
        return numpy.where( xsq <= 1, 1.0 - xsq, 0.0 )

    def partial( self, x ):
        return numpy.where( numpy.abs( x ) <= 1, -2 * x, 0.0 )

    def isBound( self ):
        return True

    def name( self ):
        return str( "Parabola: ( 1-x^2 ) if |x| < 1 else 0" )


