import numpy as numpy
import math as math
import warnings

from .Tools import setAttribute as setatt
from .Tools import printclass

__author__ = "Do Kester"
__year__ = 2021
__license__ = "GPL3"
__version__ = "2.7.0"
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
#  *    2010 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2021 Do Kester

class Prior( object ):
    """
    Base class defining prior distributions.

    Two methods need to be defined in specific priors which map
    the values between [0,1] on to the domain, and vice versa:
    unit2Domain and domain2Unit.

        u = domain2Unit( d )
        d = unit2Domain( u )

    d is a value in the domain of the prior and u is a vlue in [0,1]

    The handling of limits is relegated to this Prior class. Define
        _umin = domain2Unit( lowLimit )
        _urng = domain2Unit( highLimit ) - umin

        u = ( domain2Unit( d ) - umin ) / urange
        d = unit2Domain( u * urange + umin )

    Symmetric priors can be used in a circular variant; i.e.
    the low and high limits are folded on each other, provided
    that the limit values are the same (hence symmetric)

        u = limitedDomain2Unit( d ) + 1 ) / 3
        d = limitedUnit2Domain( ( u * 3 ) % 1 )

    The copy method is also necessary.

    Attributes
    ----------
    lowLimit : float
        low limit on the Prior
    highLimit : float
        high limit on the Prior
    deltaP : float
        width of numerical partial derivative calculation
    circular : bool or float
        whether circular

    Hidden Attributes
    -----------------
    _lowDomain : float
        lower limit of the Priors possible values
    _highDomain : float
        upper limit of the Priors possible values
    _umin : float
        umin lowLimit in unit
    _urng : float
        urange (hi-lo) in unit
    """

    #*********CONSTRUCTORS***************************************************
    def __init__( self, limits=None, circular=False, domain=[-math.inf,math.inf],
                  prior=None ):
        """
        Default constructor.

        Parameters
        ----------
        limits : None or list of 2 floats
            2 limits resp. low and high
        circular : bool or float
            False not circular
            True  circular with period from limits[0] to limits[1]
            period of circularity
        domain : 2 floats
            over which the prior is defined
        prior : Prior
            prior to copy (with new limits if applicable)
        """
        super( object, self ).__init__()

        self.deltaP = 0.0001                # for numerical partials
        # private properties of the Prior
        self._lowDomain  = domain[0]         # Lower limit of the Priors possible domain.
        self._highDomain = domain[1]         # Upper limit of the Priors possible domain

        self.setPriorAttributes( limits, circular )

        if prior is not None :
            self.deltaP = prior.deltaP

    def copy( self ) :
        """ Return a copy """
        return Prior( prior=self, limits=self.limits, circular=self.circular )

    def setLimits( self, limits=None ):
        """
        Set limits.
        It is asserted that lowLimit is smaller than highLimit.

        Parameters
        ----------
        limits : None or list of any combination of [None, float]
            None : no limit (for both or one)
            float : [low,high] limit

        Raises
        ------
        ValueError when low limit is larger than high limit or out of Domain

        """
#        print( "Prior 1 ", limits )
        if limits is None :
            self.lowLimit  = self._lowDomain
            self.highLimit = self._highDomain
            reset = False
        else :
            self.lowLimit  = limits[0] if limits[0] is not None else self._lowDomain
            self.highLimit = limits[1] if limits[1] is not None else self._highDomain
            reset = ( limits[0] is not None ) or ( limits[1] is not None )

#        print( "Prior 2 ", self.lowLimit, self.highLimit )
        if not ( self._lowDomain <= self.lowLimit < self.highLimit <= self._highDomain ) :
            print( self._lowDomain, self.lowLimit, self.highLimit, self._highDomain )
            raise ValueError( "Limits out of order or out of domain" )

        try :
            self._umin = 0.0            ## only for re-setting the limits
            self._urng = 1.0
            umax = self.domain2Unit( self.highLimit )
            self._umin = self.domain2Unit( self.lowLimit )
            self._urng = umax - self._umin
        except :
            self._urng = math.inf

#        print( "Prior 3 ", self._umin, self._urng )

        if math.isinf( self._urng ) :
#            warnings.warn( "%s needs limits" % self.__str__() )
            return

        if hasattr( self, "baseUnit2Domain" ) : return			# aleady reset

        # reset the functions
        if reset :
            setatt( self, "baseUnit2Domain", self.unit2Domain )
            setatt( self, "unit2Domain", self.limitedUnit2Domain )
            setatt( self, "baseDomain2Unit", self.domain2Unit )
            setatt( self, "domain2Unit", self.limitedDomain2Unit )

    def setPriorAttributes( self, limits, circular ) :

        if isinstance( circular, bool ) :
            self.setLimits( limits )
            if circular and limits is not None :
                period = limits[1] - limits[0]
            else :
                return
        else :
            period = circular

            if hasattr( self, "center" ) :
                self.setLimits( limits=[self.center - period / 2, self.center + period / 2] )
            else :
                self.setLimits( limits=[0,period] )

        self.circular = circular

        # reset the functions to circular variant
        setatt( self, "unit2Domain", self.circularUnit2Domain )
        setatt( self, "domain2Unit", self.circularDomain2Unit )


    def isCircular( self ) :
        """
        Whether circular
        """
        return hasattr( self, "circular" ) and self.circular

    def limitedDomain2Unit( self, dval ) :
        """
        Shrink domain to value in [0,1]
        """
        return ( self.baseDomain2Unit( dval ) - self._umin ) / self._urng

    def limitedUnit2Domain( self, uval ) :
        """
        Expand value in [0,1] to domain
        """
        return self.baseUnit2Domain( uval * self._urng + self._umin )

    def circularDomain2Unit( self, dval ) :
        return ( self.limitedDomain2Unit( dval ) + 1 ) / 3

    def circularUnit2Domain( self, uval ) :
        return self.limitedUnit2Domain( ( uval * 3 ) % 1 )

    def __setattr__( self, name, value ) :
        """
        Set attributes: lowLimit, highLimit, deltaP, _lowDomain, _highDomain.

        Also set scale for Priors that need a scale.

        """
        keys = ["lowLimit", "highLimit", "deltaP", "center", "scale",
                "_lowDomain", "_highDomain", "_umin", "_urng", "circular"]

        if name == "scale"  and value <= 0 :
            raise ValueError( "Attribute %s must be positive" % name )

        if name in keys :
            setatt( self, name, value )
        else :
            raise AttributeError( repr( self ) + " object has no attribute " + name )

    def __getattr__( self, name ) :
        """
        Return (virtual) attribute with name

        Parameters
        ----------
        name : str
            name of attribute
        """
        if name == "lowLimit" :
            return self._lowDomain
        if name == "highLimit" :
            return self._highDomain
        if name == "limits" :
            return self.getLimits()
        if name == "circular" :			## when not defined
            return False
        else :
            raise AttributeError( "Prior: Unknown attribute " + name )

    def unsetLimits( self ):
        """ Remove all limits.  """
        self.lowLimit = self._lowDomain
        self.highLimit = self._highDomain

    def setAttributes( self, limits=None, scale=None ) :
        """
        Set possible attributes for a Prior.

        Parameters
        ----------
        limits : float or None
            [low,high] limit
        scale : float or None
            scale factor
        """
        if limits is not None :
            self.setLimits( limits=limits )
        if scale is not None :
            self.scale = scale

    def isOutOfLimits( self, par ):
        """
        True if the parameter is out of limits

        Parameters
        ----------
        par : float
            the parameter to check

        """
        return ( par < self.lowLimit ) or ( par > self.highLimit )

    def checkLimit( self, par ):
        """
        Check whether the parameter is within limits.

        Parameters
        ----------
        par : float
            the parameter to check

        Raises
        ------
            ValueError when outside limits.

        """
        if self.isOutOfLimits( par ):
            raise ValueError( "Parameter outside supplied limits: %8.2f < %8.2f < %8.2f"%
                            (self.lowLimit, par, self.highLimit) )

    def stayInLimits( self, par ):
        """
        Return lower limit or upper limit when parameter is outside.

        Parameters
        ----------
        par : float
            the parameter to check

        """
        if par < self.lowLimit:
            return self.lowLimit
        if par > self.highLimit:
            return self.highLimit
        return par

    def hasLowLimit( self ):
        """ Return true if the prior has its low limits set.  """
        return self.lowLimit > self._lowDomain

    def hasHighLimit( self ):
        """ Return true if the prior has its high limits set.  """
        return self.highLimit < self._highDomain

    def hasLimits( self ):
        """ Return true if it has any limits.  """
        return self.hasLowLimit() and self.hasHighLimit()

    def getLimits( self ):
        """ Return the limits.  """
        return numpy.array( [ self.lowLimit, self.highLimit ] )

    def getIntegral( self ) :
        """
        Return the integral of the prior over the valid range.

        Default: 1.0 (for bound priors)
        """
        return 1.0

    def getRange( self ):
        """ Return the range.  """
        return self.highLimit - self.lowLimit


    def domain2Unit( self, dval ):
        """
        Return a value in [0,1] given a value within the valid domain of
        a parameter for a distribution.

        Parameters
        ----------
        dval : float
            value within the domain of a parameter

        """
        pass

    def unit2Domain( self, uval ):
        """
        Return a value within the valid domain of the parameter given a value
        between [0,1] for a distribution.

        Parameters
        ----------
        uval : float
            value within [0,1]

        """
        pass

    def result( self, p ):
        """
        Return value of the Prior at a given value.

        If result is not defined, fall back to numerical derivative od Domain2Unit.

        Parameters
        ----------
        p : float
            the value

        """
        return self.numPartialDomain2Unit( p )

    def partialDomain2Unit( self, p ):
        """
        Return the derivative of Domain2Unit, aka the result of the distribution at p

        Parameters
        ----------
        p : float
            the value

        """
        return self.result( p )

    def logResult( self, p ) :
        """
        Return the log of the result; -inf when p == 0.

        Parameters
        ----------
        p : float
            the value

        """
        try :
            return math.log( self.result( p ) )
        except :
            return -math.inf

    def numPartialDomain2Unit( self, dval ):
        """
        Return a the numeric derivate of the domain2Unit function to dval.

        Parameters
        ----------
        dval : float
            value within the domain of a parameter

        """
        return ( ( self.domain2Unit( dval + self.deltaP ) -
                   self.domain2Unit( dval - self.deltaP ) ) /
                        ( 2 * self.deltaP ) )

    def partialLog( self, p ):
        """
        Return partial derivative of log( Prior ) wrt parameter.
        default numPartialLog

        Parameters
        ----------
        p : float
            the value

        """
        return self.numPartialLog( p )

    def numPartialLog( self, p ):
        """
        Return the numeric partial derivative of log( Prior ) wrt parameter.
        Parameters
        ----------
        p : float
            the value

        """
        if self.isOutOfLimits( p ) : return math.nan

        rm = self.logResult( p - self.deltaP )
        rp = self.logResult( p + self.deltaP )
        return ( rp - rm ) / ( 2 * self.deltaP )

    def isBound( self ):
        """ Return true if the integral over the prior is bound.  """
        return False

    def __str__( self ):
        """ Return a string representation of the prior.  """

        name = self.shortName()

        if hasattr( self, "center" ) :
            name += str( " at (%.2f %.2f)" % ( self.center, self.scale ) )

        if self.hasLowLimit() or self.hasHighLimit():
            name += str( " circular" if self.isCircular() else " with limits" )
            name += str( " between %.2f and %.2f" % ( self.lowLimit, self.highLimit ) )

        return name

#      * End of Prior


