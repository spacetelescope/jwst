import numpy as numpy
from . import Tools
from .Tools import setAttribute as setatt
from .LinearModel import LinearModel

#__author__ = "Do Kester"
#__year__ = 2020
#__license__ = "GPL3"
#__version__ = "2.5.3"
#__url__ = "https://www.bayesicfitting.nl"
#__status__ = "Perpetual Beta"

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
#  *    2003 - 2014 Do Kester, SRON (JAVA code)
#  *    2016 - 2020 Do Kester

class SplinesModel( LinearModel ):
    """
    General splines model of arbitrary order and with arbitrary knot settings.
    It is a linear model.

    order   behaviour between knots     continuity at knots
      0     piecewise constant          not continuous at all
      1     piecewise linear            lines are continuous (connected)
      2     parabolic pieces            1st derivatives are also continuous
      3     cubic pieces                2nd derivatives are also continuous
     n>3    n-th order polynomials      (n-1)-th derivatives are also continuous

    The user lays out a number ( << datapoints ) of knots on the x-axis at
    arbitrary position, generally more knots where the curvature is higher.
    The knots need to be monotonuously increasing in x.
    Alternatively one can ask this class to do the lay-out which is then
    equidistant in x over the user-provided range.
    Through these knots a splines function is obtained which best
    fits the datapoints. One needs at least 2 knots, one smaller and one
    larger than the x-values in the dataset.

    If the end knots are put in between the x-values in the dataset, a kind of
    extrapoling spline is obtained. It still works more or less. Dont push it.

    This model is NOT for (cubic) spline interpolation.

    Examples
    --------
    >>> knots = numpy.arange( 17, dtype=float ) * 10    # make equidistant knots from 0 to 160
    >>> csm = SplinesModel( knots=knots, order=2 )
    >>> print csm.getNumberOfParameters( )
    18
    # or alternatively:
    >>> csm = SplinesModel( nrknots=17, order=2, min=0, max=160 )    # automatic layout of knots
    >>> print csm.getNumberOfParameters( )
    18
    # or alternatively:
    >>> npt = 161                                               # to include both 0 and 160.
    >>> x = numpy.arange( npt, dtype=float )                    # x-values
    >>> csm = SplinesModel( nrknots=17, order=2, xrange=x )     # automatic layout of knots
    >>> print csm.getNumberOfParameters( )
    18

    Attributes
    ----------
    knots : array_like
        positions of the spline knots
    order : int
        order of the spline. default: 3

    Attributes from Model
    ---------------------
        npchain, parameters, stdevs, xUnit, yUnit

    Attributes from FixedModel
    --------------------------
        npmax, fixed, parlist, mlist

    Attributes from BaseModel
    --------------------------
        npbase, ndim, priors, posIndex, nonZero, tiny, deltaP, parNames


    Limitations
    -----------
    Dont construct the knots so closely spaced, that there are no datapoints in between.

    """
    def __init__( self, knots=None, order=3, nrknots=None, min=None, max=None, xrange=None,
                        copy=None, **kwargs ):
        """
        Splines on a given set of knots and a given order.

        The number of parameters is ( length( knots ) + order - 1 )

        Parameters
        ----------
        knots : array_like
            a array of arbitrarily positioned knots
        order : int
            order of the spline. Default 3 (cubic splines)
        nrknots : int
            number of knots, equidistantly posited over xrange or [min,max]
        min : float
            minimum of the knot range
        max : float
            maximum of the knot range
        xrange : array_like
            range of the xdata
        copy : SplinesModel
            model to be copied.
        fixed : None or dictionary of {int:float|Model}
            int         index of parameter to fix permanently.
            float|Model values for the fixed parameters.
            Attribute fixed can only be set in the constructor.
            See: @FixedModel

        Raises
        ------
        ValueError : At least either (`knots`) or (`nrknots`, `min`, `max`) or
                (`nrknots`, `xrange`) must be provided to define a valid model.

        Notes
        -----
        The SplinesModel is only strictly valid inside the domain defined by the
        minmax of knots. It deteriorates fastly going outside the domain.

        """

        if copy is not None :
            knots = copy.knots
            order = copy.order
        if knots is not None : nrknots = len( knots )
        elif nrknots is None :
            raise ValueError( "Need either knots or (nrknots,min,max) or (nrknots,xrange)" )

        if 'nparams' in kwargs :
            npar = kwargs['nparams']
            del kwargs['nparams']
        else :
            npar = order + nrknots - 1

        super( SplinesModel, self ).__init__( npar, copy=copy, **kwargs )
        self.order = order
        if knots is None :
            if xrange is not None :
                min = numpy.min( xrange )
                max = numpy.max( xrange )
            knots = numpy.linspace( min, max, nrknots, dtype=float )
        self.knots = knots

    def copy( self ):
        return SplinesModel( copy=self )

    def __setattr__( self, name, value ):
        """
        Set attributes: knots, order

        """
        if name == "knots" :
            setatt( self, name, value, type=float, islist=True )
        elif name == "order" :
            setatt( self, name, value, type=int )
        else :
            super( SplinesModel, self ).__setattr__( name, value )

    def basePartial( self, xdata, params, parlist=None ):
        """
        Returns the partials at the input value.

        The partials are the powers of x (input) from 0 to degree.

        Parameters
        ----------
        xdata : array_like
            value at which to calculate the partials
        params : array_like
            parameters to the model (ignored in LinearModels)
        parlist : array_like
            list of indices active parameters (or None for all)

        """
        np = self.npmax
        ni = Tools.length( xdata )
        partial = numpy.zeros( ( ni, np), dtype=float )
        x = numpy.ones( ni )
        for i in range( self.order ):
            partial[:,i] = x
            x *= xdata

        i = self.order
        ko = i - 1

        ki = numpy.arange( ni )
        while i < np - 1 :
            partial[ki,i] = x[ki]

            ki = numpy.where( xdata >= self.knots[i - ko] )[0]
            if ki.size == 0 : break

            d = xdata[ki] - self.knots[i-ko]
            t = numpy.ones( len( d ) )
            for k in range( self.order ):
                t *= d
            partial[ki,i] -= t

            x[ki] = t
            i += 1
        partial[ki,i] = x[ki]

        if parlist is None or len( parlist ) == np :
            return partial

        return partial[:,parlist]


    def baseDerivative( self, xdata, params ) :
        """
        Return the derivative df/dx at each xdata (=x).

        Parameters
        ----------
        xdata : array_like
            value at which to calculate the partials
        params : array_like
            parameters to the model

        """
        ff = numpy.append( numpy.arange( self.order + 1, dtype=float ),
                           numpy.full( self.npbase - self.order - 1, 3, dtype=float ) )
        pars = params * ff
        pars = pars[1:]
        return SplinesModel( knots=self.knots, order=self.order-1 ).result( xdata, pars )


    def baseName( self ):
        """ Returns a string representation of the model. """
        return "Splines of order %d with %d knots."%( self.order, len( self.knots) )

    def baseParameterUnit( self, k ):
        """
        Return the name of the parameter.

        Parameters
        ----------
        k : int
            index of the parameter.
        """
        if k > self.order :
            k = self.order
        return self.yUnit / ( self.xUnit ** k )


