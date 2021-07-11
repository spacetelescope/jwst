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
#  *    2010 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2017 Do Kester

class Gauss( Kernel ):
    """
    Gauss is an unbound Kernel function:

        f( x ) = exp( -0.5 * x * x ).

    """
    def __init__( self ) :
        """
        Constructor.

        Using
            integral = sqrt( 2 &pi; )
            fwhm = sqrt( 2 log( 2 ) )
            range = inf
        """
        integral = numpy.sqrt( 2 * math.pi )
        fwhm = 2 * math.sqrt( 2 * math.log( 2 ) )
        super( Gauss, self ).__init__( integral=integral, fwhm=fwhm, range=math.inf )

    def integral( self ):
        return

    def result( self, x ):
        return self.resultsq( x * x )

    def resultsq( self, xsq ):
        return numpy.exp( -0.5 * xsq )

    def partial( self, x ):
        return -x * self.result( x )

    def isBound( self ):
        return False

    def name( self ):
        return str( "Gauss: exp( -0.5 * x^2 )" )


