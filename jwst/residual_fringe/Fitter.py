import numpy as numpy
from .BaseFitter import BaseFitter

from .Formatter import formatter as fmt

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
#  * A JAVA version of this code was part of the Herschel Common
#  * Science System (HCSS), also under GPL3.
#  *
#  *    2003 - 2014 Do Kester, SRON (JAVA code)
#  *    2016 - 2020 Do Kester

class Fitter( BaseFitter ):
    """
    Fitter for linear models.

    The Fitter class is to be used in conjunction with Model classes.

    The Fitter class and its descendants fit data to a model. Fitter itself
    is the variant for linear models, ie. models linear in its parameters.

    Examples
    --------
    # assume x and y are numpy.asarray data arrays:
    >>> x = numpy.arange( 100 )
    >>> y = numpy.arange( 100 ) // 4        # digitization noise
    >>> poly = PolynomialModel( 1 )         # line
    >>> fitter = Fitter( x, poly )
    >>> param = fitter.fit( y )
    >>> stdev = fitter.stdevs               # stdevs on the parameters
    >>> chisq = fitter.chisq
    >>> scale = fitter.scale                # noise scale
    >>> yfit  = fitter.getResult( )         # fitted values
    >>> yfit  = poly( x )                   # same as previous
    >>> yband = fitter.monteCarloError( )        # 1 sigma confidence region


    Limitations
    -----------
    1. The Fitter does not work with limits.
    2. The calculation of the evidence is an Gaussian approximation which is
       only exact for linear models with a fixed scale.

    Author  Do Kester

    """

    def __init__( self, xdata, model, map=False, keep=None, fixedScale=None ):
        """
        Create a new Fitter, providing xdatas and model.

        A Fitter class is defined by its model and the input vector (the
        independent variable). When a fit to another model and/or another
        input vector is needed a new object should be created.

        Parameters
        ----------
        xdata : array_like
            array of independent input values
        model : Model
            the model function to be fitted
        map : bool (False)
            When true, the xdata should be interpreted as a map.
            The fitting is done on the pixel indices of the map,
            using ImageAssistant
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)
            The values of keep will be used by the Fitter as long as the Fitter exists.
            See also `fit( ..., keep=dict )`
        fixedScale : float
            the fixed noise scale

        """
        super( Fitter, self ).__init__( xdata, model, map=map, keep=keep, fixedScale=fixedScale )

    def fit( self, ydata, weights=None, keep=None, plot=False ):
        """
        Return model parameters fitted to the data, including weights.

        For Linear models the matrix equation

            H * p = &beta;

        is solved for p. H is the Hessian matrix ( D * w * D^T )
        and &beta; is the inproduct of the data with the D, design matrix.

            &beta; = y * w * D^T

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted
        weights : array_like
            weights pertaining to the data ( = 1.0 / sigma^2 )
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)
            The values will override those at initialization.
            They are only used in this call of fit.
        plot : bool
            Plot the results

        Raises
        ------
            ValueError when ydata or weights contain a NaN

        """
        fitIndex, ydata, weights = self.fitprolog( ydata, weights=weights, keep=keep )

        if self.model.isNullModel() :
            self.chiSquared( ydata, weights )
            return numpy.asarray( 0 )

        hessian = self.getHessian( weights=weights, index=fitIndex )
        ydatacopy = ydata.copy( )
        # subtract influence of fixed parameters on the data
        if fitIndex is not None :
            fxpar = numpy.copy( self.model.parameters )
            fxpar[fitIndex] = 0.0
            ydatacopy = numpy.subtract( ydatacopy, self.model.result( self.xdata, fxpar) )

        if weights is not None :
            ydatacopy *= weights
        if hasattr( self, "normdfdp" ) :
            ydatacopy = numpy.append( ydatacopy, self.normdata * self.normweight )

        vector = self.getVector( ydatacopy, index=fitIndex )
#        print( fmt( hessian ) )
        params = numpy.linalg.solve( hessian, vector )

        params = self.insertParameters( params, index=fitIndex )
        self.model.parameters = params
        self.chiSquared( ydata, params=params, weights=weights)

        self.fitpostscript( ydata, plot=plot )

        return params

    def __str__( self ):
        """ Return the name of the fitter. """
        return "Fitter"


