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
#  *    2013 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2017 Do Kester

class Tricube( Kernel ):
    """
    Tricube is a Kernel function between [-1,1]; it is 0 elsewhere.

        K( x ) = ( 1 - |x|^3 )^3        if |x| < 1
                 0                      elsewhere

    """
    def __init__( self ) :
        """
        Constructor.

        Using
            integral = 81.0/70.0
            fwhm = 1.18176021
        """
        super( Tricube, self ).__init__( integral=81.0/70.0, fwhm=1.18176021 )

    def result( self, x ):
        ax = numpy.abs( x )
        res = numpy.where( ax <= 1, 1 - ax * ax * ax, 0.0 )
        return res * res * res

    def partial( self, x ):
        ax = numpy.abs( x )
        u = 1.0 - ax * ax * ax
        return numpy.where( ax <= 1, -9 * u * u * ax * ax * numpy.sign( x ), 0.0 )

    def isBound( self ):
        return True

    def name( self ):
        return str( "Tricube: ( 1 - |x|^3 )^3 if |x| < 1 else 0" )


