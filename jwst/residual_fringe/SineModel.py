import numpy as numpy
from astropy import units
import math
from . import Tools
from .Tools import setAttribute as setatt

from .NonLinearModel import NonLinearModel

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
#  *    2003 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2020 Do Kester

class SineModel( NonLinearModel ):
    """
    Sinusoidal Model.

    Two variants are implemented.

    1. By default it is the weighted sum of sine and cosine of the same frequency:

        f( x:p ) = p_1 * cos( 2 * &pi; * p_0 * x ) + p_2 * sin( 2 * &pi; * p_0 * x )

    where
        p_0 = frequency
        p_1 = amplitude cosine and
        p_2 = amplitude sine.
    As always x = input.

    The parameters are initialized at [1.0, 1.0, 1.0]. It is a non-linear model.

    2. If phase == True, the sinusoidal model has an explicit phase:

        f( x:p ) = p_0 * sin( 2 * &pi; * p_1 * x + p_2 )

    where
        p_0 = amplitude
        p_1 = frequency
        p_2 = phase.

    The parameters are initialized as [1.0, 1.0, 0.0].


    Examples
    --------
    >>> sine = SineModel( )
    >>> print( sine.npchain )
    3
    >>> pars = [0.1, 0.0, 1.0]
    >>> sine.parameters = pars
    >>> print( sine( numpy.arange( 11, dtype=float ) ) )    # One sine period
    >>> pars = [0.1, 1.0, 0.0]
    >>> sine.parameters = pars
    >>> print( sine( numpy.arange( 11, dtype=float ) ) )     # One cosine period

    Attributes
    ----------
    phase : bool (False)
        False : original 2 amplitudes model
        True  : phase model

    Attributes from Model
    ---------------------
        npchain, parameters, stdevs, xUnit, yUnit

    Attributes from FixedModel
    --------------------------
        npmax, fixed, parlist, mlist

    Attributes from BaseModel
    --------------------------
        npbase, ndim, priors, posIndex, nonZero, tiny, deltaP, parNames


    """
    TWOPI = 2 * math.pi

    def __init__( self, copy=None, phase=False, **kwargs ):
        """
        Sinusiodal model.

        Number of parameters is 3.

        Parameters
        ----------
        phase : bool
            if True, construct phase variant.
        copy : SineModel
            model to copy
        fixed : dictionary of {int:float}
            int     list if parameters to fix permanently. Default None.
            float   list of values for the fixed parameters.
            Attribute fixed can only be set in the constructor.

        """
        if phase :
            param = [1.0, 1.0, 0.0]
            names = ["amplitude", "frequency", "phase"]
        else :
            param = [1.0, 1.0, 1.0]
            names = ["frequency", "cosamp", "sinamp"]

        super( SineModel, self ).__init__( 3, copy=copy, params=param,
                        names=names, **kwargs )

        setatt( self, "phase", phase )
        if phase :
            setatt( self, "baseResult", self.phaseResult )
            setatt( self, "basePartial", self.phasePartial )
            setatt( self, "baseDerivative", self.phaseDerivative )
            setatt( self, "baseName", self.phaseName )
            setatt( self, "baseParameterUnit", self.phaseParameterUnit )

    def copy( self ):
        """ Copy method.  """
        return SineModel( phase=self.phase, copy=self )


    def baseResult( self, xdata, params ):
        """
        Returns the result of the model function.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.

        """
        x = self.TWOPI * xdata * params[0]
        result = params[1] * numpy.cos( x ) + params[2] * numpy.sin( x )
        return result

    def basePartial( self, xdata, params, parlist=None ):
        """
        Returns the partials at the input value.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.
        parlist : array_like
            list of indices active parameters (or None for all)

        """
        np = self.npbase if parlist is None else len( parlist )
        partial = numpy.ndarray( ( Tools.length( xdata ), np ) )

        #  disregard count
        x = self.TWOPI * xdata
        xf = x * params[0]
        cxf = numpy.cos( xf )
        sxf = numpy.sin( xf )

        parts = { 0 : ( lambda: x * params[2] * cxf - x * params[1] * sxf ),
                  1 : ( lambda: cxf ),
                  2 : ( lambda: sxf ) }

        if parlist is None :
            parlist = range( self.npmax )

        for k,kp in enumerate( parlist ) :
            partial[:,k] = parts[kp]()

        return partial

    def baseDerivative( self, xdata, params ):
        """
        Returns the derivative of f to x (df/dx) at the input values.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.

        """
        x = self.TWOPI * xdata * params[0]
        df = self.TWOPI * params[0] * ( params[2] * numpy.cos( x ) - params[1] * numpy.sin( x ) )
        return df

    def baseName( self ):
        """
        Returns a string representation of the model.

        """
        return str( "Sine: f( x:p ) = p_1 * cos( 2PI * x * p_0 ) + p_2 * sin( 2PI * x * p_0 )" )

    def baseParameterUnit( self, k ):
        """
        Return the unit of a parameter.

        Parameters
        ----------
        k : int
            the kth parameter.

        """
        if k == 0:
            return units.Unit( units.si.rad ) / self.xUnit
        return self.yUnit

    ###### PHASE VARIANT #################################################

    def phaseResult( self, xdata, params ):
        """
        Returns the result of the model function.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.

        """
        x = self.TWOPI * xdata * params[1] + params[2]
        result = params[0] * numpy.sin( x )
        return result


    def phasePartial( self, xdata, params, parlist=None ):
        """
        Returns the partials at the input value.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.
        parlist : array_like
            list of indices active parameters (or None for all)

        """
        np = self.npbase if parlist is None else len( parlist )
        partial = numpy.ndarray( ( Tools.length( xdata ), np ) )

        #  disregard count
        x = self.TWOPI * xdata
        xf = x * params[1] + params[2]
        cxf = params[0] * numpy.cos( xf )

        parts = { 0 : ( lambda: numpy.sin( xf ) ),
                  1 : ( lambda: cxf * x ),
                  2 : ( lambda: cxf ) }

        if parlist is None :
            parlist = range( self.npmax )

        for k,kp in enumerate( parlist ) :
            partial[:,k] = parts[kp]()

        return partial

    def phaseDerivative( self, xdata, params ):
        """
        Returns the derivative of f to x (df/dx) at the input values.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.

        """
        x = self.TWOPI * xdata * params[1] + params[2]
        df = params[0] * numpy.cos( x ) * self.TWOPI * params[1]
        return df

    def phaseName( self ):
        """
        Returns a string representation of the model.

        """
        return str( "Sine: f( x:p ) = p_0 * sin( 2PI * x * p_1 + p_2 )" )

    def phaseParameterUnit( self, k ):
        """
        Return the unit of a parameter.

        Parameters
        ----------
        k : int
            the kth parameter.

        """
        if k == 0:
            return self.yUnit
        if k == 1:
            return units.Unit( units.si.rad ) / self.xUnit
        return units.Unit( units.si.rad )


