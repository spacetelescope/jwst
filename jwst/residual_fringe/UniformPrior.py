import math

from .Prior import Prior

__author__ = "Do Kester"
__year__ = 2020
__license__ = "GPL3"
__version__ = "2.6.2"
__url__ = "https://www.bayesicfitting.nl"
__status__ = "Perpetual Beta"

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
#  *    2003 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2020 Do Kester

class UniformPrior( Prior ):
    """
    Uniform prior distribution, for location parameters.

    A uniform prior is a improper prior ( i.e. its integral is unbound ).
    Because of that it always needs limits, low and high, such that
    -Inf < low < high < +Inf.

        Pr( x ) = 1 / ( high - low )    if low < x < high
                  0                     elsewhere

    domain2Unit: u = ( d - lo ) / range
    unit2Domain: d = u * range + lo

    Examples
    --------
    >>> pr = UniformPrior()                                 # unbound prior
    >>> pr = UniformPrior( limits=[0,10] )                  # limited to the range [0,10]
    >>> pr = UniformPrior( circular=math.pi )               # circular between 0 and pi
    >>> pr = UniformPrior( limits=[2,4], circular=True )    # circular between 2 and 4

    Attributes
    ----------
    _range : float
        highlimit - lowlimit

    Attributes from Prior
    --------------------=
    lowLimit, highLimit, deltaP, _lowDomain, _highDomain


    """

    #  *********CONSTRUCTORS***************************************************
    def __init__( self, limits=None, circular=False, prior=None ):
        """
        Constructor.

        Parameters
        ----------
        limits : None or [float,float]
            None    no limits are set
            2 floats    lowlimit and highlimit
        circular : bool or float
            True : circular with period from limits[0] to limits[1]
            float : period of circularity
        prior : UniformPrior
            to be copied
        """
        super( ).__init__( limits=limits, circular=circular, prior=prior )

    def copy( self ):
        """ Return a (deep) copy of itself. """
        return UniformPrior( prior=self, limits=self.limits, circular=self.circular )

    def getIntegral( self ) :
        """
        Return integral of UniformPrior from lowLimit to highLimit.
        """
        return self._urng

    def domain2Unit( self, dval ):
        """
        Return the dval as uval

        In Prior.limitedDomain2Unit the dval is transformed into a uval

        Parameters
        ----------
        dval : float
            value within the domain of a parameter

        """
        return dval

    def unit2Domain( self, uval ):
        """
        Return the uval as dval

        In Prior.limitedUnit2Domain the uval is transformed into a dval

        Parameters
        ----------
        uval : float
            value within [0,1]

        """
        return uval

    def result( self, x ):
        """
        Return a the result of the distribution function at x.

        Parameters
        ----------
        x : float
            value within the domain of a parameter

        """
        if math.isinf( self._urng ) :
            raise AttributeError( "Limits are needed for UniformPrior" )

        return 0.0 if self.isOutOfLimits( x ) else 1.0 / self._urng


# logResult has no better definition than the default: just take the math.log of result.
# No specialized method here.

    def partialLog( self, p ):
        """
        Return partial derivative of log( Prior ) wrt parameter.

        Parameters
        ----------
        p : float
            the value

        """
        return math.nan if self.isOutOfLimits( p ) else 0

    def isBound( self ):
        """ Return true if the integral over the prior is bound.  """
        return self.hasLowLimit( ) and self.hasHighLimit( )

    def shortName( self ):
        """ Return a string representation of the prior.  """
        return str( "UniformPrior" + ( " unbound." if not self.isBound( ) else "" ) )




