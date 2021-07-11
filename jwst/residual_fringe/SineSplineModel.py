import numpy as numpy
import math
from . import Tools
from .Tools import setAttribute as setatt
from .LinearModel import LinearModel
from .SplinesModel import SplinesModel

__author__ = "Do Kester"
__year__ = 2020
__license__ = "GPL3"
__version__ = "2.5.3"
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
#  *    2017 - 2020 Do Kester

class SineSplineModel( LinearModel ):
    """
    Sine of fixed frequency with splineslike amplitudes/phases.

        f( x:p ) = SM0 \cos( 2 &pi; &omega; x ) + SM1 sin( 2 &pi; &omega; x )

    Where SM0 and SM1 are splines models with defined knots and order.

    It is a linear model with 2 * ( len(knots) + order - 1 ) papameters.

    Examples
    --------
    >>> knots = [3.0*k for k in range( 11 )]
    >>> sine = SineSplineModel( 150, knots )        # fixed frequency of 150 Hz
    >>> print( sine.npbase )                        # number of parameters
    26

    Attributes
    ----------
    frequency : float
        (fixed) frequency of the sine
    knots : array_like
        positions of the spline knots
    order : int
        order of the spline. default: 3
    cm : SplinesModel
        amplitude of the cosine
    sm : SplinesModel
        amplitude of the sine

    Attributes from Model
    ---------------------
        npchain, parameters, stdevs, xUnit, yUnit

    Attributes from FixedModel
    --------------------------
        npmax, fixed, parlist, mlist

    Attributes from BaseModel
    --------------------------
        npbase, ndim, priors, posIndex, nonZero, tiny, deltaP, parNames


    Alternate
    ---------
    The model

    >>> model = SineSplineModel( frequency, knots )

    is equivalent to :

    >>> cm = SplinesModel( knots )
    >>> sm = SplinesModel( knots )
    >>> fxd = {0:cm, 1:sm}
    >>> model = SineAmpModel( frequency, fixed=fxd )


    """

    def __init__( self, frequency, knots, order=3, copy=None, fixed=None, **kwargs ):
        """
        Sine model of a fixed frequency with a splineslike changing amplitude/phase.

        Number of parameters is 2 * ( len(knots) + order - 1 ).

        Parameters
        ----------
        frequency : float
            the frequency
        copy : SineSplineModel
            model to be copied
        fixed : dict
            If not None raise AttributeError.

        Raises
        ------
        AttributeError
            when fixed is not None

        """
        if fixed is not None :
            raise AttributeError( "SineSplineModel cannot have fixed parameters" )

        np = 2 * ( len( knots ) + order - 1 )
        super( SineSplineModel, self ).__init__( np, copy=copy, **kwargs )

        self.frequency = frequency
        self.knots = knots
        self.order = order
        if copy is not None :
            setatt( self, "cm", copy.cm.copy() )
            setatt( self, "sm", copy.sm.copy() )
        else :
            setatt( self, "cm", SplinesModel( knots, order=order ) )
            setatt( self, "sm", SplinesModel( knots, order=order ) )

    def copy( self ):
        """ Copy method.  """
        return SineSplineModel( self.frequency, self.knots, order=self.order, copy=self )

    def __setattr__( self, name, value ) :
        if name == "frequency" :
            setatt( self, name, value, type=float )
        elif name == "knots" :
            setatt( self, name, value, type=float, islist=True )
        elif name == "order" :
            setatt( self, name, value, type=int )
        elif name == "cm" or name == "sm" :
            raise AttributeError( "Attributes cm or sm can not be set" )
        else :
            super( SineSplineModel, self ).__setattr__( name, value )

    def basePartial( self, xdata, params, parlist=None ):
        """
        Returns the partials at the input value.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the partials
        params : array_like
            parameters of the model (ignored in LinearModels)
        parlist : array_like
            list of indices active parameters (or None for all)

        """
        nxdata = Tools.length( xdata )
        part = numpy.zeros( ( nxdata, self.npbase ), dtype=float )
        nh = self.npbase // 2
        cx = numpy.cos( 2 * math.pi * self.frequency * xdata )
        sx = numpy.sin( 2 * math.pi * self.frequency * xdata )

        part[:,:nh] = ( cx * self.cm.basePartial( xdata, params ).transpose() ).transpose()
        part[:,nh:] = ( sx * self.sm.basePartial( xdata, params ).transpose() ).transpose()
        return part

    def baseDerivative( self, xdata, params ):
        """
        Returns the derivative of f to x (df/dx) at the input value.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the partials
        params : array_like
            parameters of the model

        """
        tpf = 2 * math.pi * self.frequency
        tx = tpf * xdata
        cx = numpy.cos( tx )
        sx = numpy.sin( tx )
        amps = self.getAmplitudes( xdata, params )
        cadx = self.cm.baseDerivative( xdata, params )
        sadx = self.sm.baseDerivative( xdata, params )
        return tpf * ( cadx * cx - amps[0] * sx + sadx * sx + amps[1] * cx )

    def getAmplitudes( self, xdata, params ) :
        """
        Return the amplitudes if cosine and sine, resp.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the partials
        params : array_like
            parameters of the model

        """
        nh = self.npbase // 2
        return ( self.cm.result( xdata, params[:nh] ),
                 self.sm.result( xdata, params[nh:] ) )

    def baseName( self ):
        """
        Returns a string representation of the model.

        """
        return ( "SineSpline: f( x:p ) = spline_0 * cos( 2 pi * x * f ) + " +
                 "spline_1 * sin( 2 pi * x * f ); f = %f"%self.frequency )

    def baseParameterUnit( self, k ):
        """
        Return the name of a parameter.
        Parameters
        ----------
        k : int
            the kth parameter.

        """
        k = k % ( self.npbase / 2 )
        if k > self.order :
            k = self.order
        return self.yUnit / ( self.xUnit ** k )



