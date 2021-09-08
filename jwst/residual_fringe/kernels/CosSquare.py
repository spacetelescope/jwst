import numpy as numpy

from .Kernel import Kernel
from .Cosine import Cosine

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
#  *    2016 - 2017 Do Kester

class CosSquare( Kernel ):
    """
    CosSquare (Cosine Squared) is a Kernel function between [-1,1]; it is 0 elsewhere.

        K( x ) = cos^2( 0.5 &pi; x )    if |x| < 1
                 0                      elsewhere


    """
    def __init__( self ) :
        """
        Constructor.

        Using
            integral = 1.0
            fwhm = 1.0
        """
        super( CosSquare, self ).__init__( integral=1.0, fwhm=1.0 )

    def result( self, x ):
        res = Cosine.result( Cosine, x )
        return res * res

    def resultsq( self, xsq ):
        return self.result( numpy.sqrt( xsq ) )

    def partial( self, x ):
        return 2 * Cosine.result( Cosine, x ) * Cosine.partial( Cosine, x )

    def isBound( self ):
        return True

    def name( self ):
        return str( "CosSquare: cos^2( 0.5*PI*x ) if |x| < 1 else 0" )


