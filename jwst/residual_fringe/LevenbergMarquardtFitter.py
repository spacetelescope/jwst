import numpy as numpy
from astropy import units
import math
from . import Tools

from .IterativeFitter import IterativeFitter
from .ConvergenceError import ConvergenceError
from .Formatter import formatter as fmt

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
#  *    2017 - 2020 Do Kester


class LevenbergMarquardtFitter( IterativeFitter ):
    """
    Non-linear fitter using the Levenberg-Marquardt method.

    Implementation of the Levenberg-Marquardt algorithm to fit the parameters
    of a non-linear model. It is a gradient fitter which uses partial
    derivatives to find the downhill gradient. Consequently it ends in the
    first minimum it finds.

    The original C-version stems from Numerical Recipes with some additions of my own.
    This might be the third or fourth transcription of it.

    Author       Do Kester.

    Examples
    --------
    # assume x and y are Double1d data arrays:
    >>> x = numpy.arange( 100, dtype=float ) / 10
    >>> y = numpy.arange( 100, dtype=float ) / 122            # make slope
    >>> rg = RandomGauss( seed=12345L )            # Gaussian random number generator
    >>> y += rg( numpy.asarray( 100, dtype=float ) ) * 0.2            # add noise
    >>> y[Range( 9,12 )] += numpy.asarray( [5,10,7], dtype=float )         # make some peak
    # define a model: GaussModel + background polynomial
    >>> gauss = GaussModel( )                            # Gaussian
    >>> gauss += PolynomialModel( 1 )                    # add linear background
    >>> gauss.setParameters( numpy.asarray( [1,1,0.1,0,0], dtype=float ) )    # initial parameter guess
    >>> print gauss.getNumberOfParameters( )                # 5 ( = 3 for Gauss + 2 for line )
    >>> gauss.keepFixed( numpy.asarray( [2] ), numpy.asarray( [0.1], dtype=float ) )    # keep width fixed at 0.1
    >>> lmfit = LevenbergMarquardtFitter( x, gauss )
    >>> param = lmfit.fit( y )
    >>> print param.length( )                             # 4 ( = 5 - 1 fixed )
    >>> stdev = lmfit.getStandardDeviation( )             # stdevs on the parameters
    >>> chisq = lmfit.getChiSquared( )
    >>> scale = lmfit.getScale( )                         # noise scale
    >>> yfit  = lmfit.getResult( )                        # fitted values
    >>> yband = lmfit.monteCarloError( )                       # 1 sigma confidence region
    # for diagnostics ( or just for fun )
    >>> lmfit = LevenbergMarquardtFitter( x, gauss )
    >>> lmfit.setVerbose( 2 )                             # report every 100th iteration
    >>> plotter = IterationPlotter( )                     # from BayesicFitting
    >>> lmfit.setPlotter( plotter, 20 )                   # make a plot every 20th iteration
    >>> param = lmfit.fit( y )

    Notes
    -----
    In case of problems look at the "Troubles" page in the documentation area.


    Limitations
    -----------
    1. LMF is <b>not</b> guaranteed to find the global minimum.
    2. The calculation of the evidence is an Gaussian approximation which is
    only exact for linear models with a fixed scale.

    Attributes
    ----------
    xdata : array_like
        vector of numbers as input for model
    model : Model
        the model to be fitted
    lamda : float
        to balance the curvature matrix (see Numerical Recipes)

    """
    #  *************************************************************************
    def __init__( self, xdata, model, **kwargs ):
        """
        Create a class, providing xdata and model.

        Parameters
        ----------
        xdata : array_like
            vector of independent input values
        model : Model
            a model function to be fitted

        kwargs : dict
            Possibly includes keywords from
                IterativeFitter :       maxIter, tolerance, verbose
                BaseFitter :            map, keep, fixedScale


        """
        super( LevenbergMarquardtFitter, self ).__init__( xdata, model, **kwargs )

        self.lamda = 0.001
        self.converged = False
        self.first = True

    #  *************************************************************************
    def fit( self, data, weights=None, par0=None, keep=None, limits=None,
                maxiter=None, tolerance=None, verbose=None, plot=False,
                callback=None ):
        """
        Return Model fitted to the data arrays.

        It will calculate the hessian matrix and chisq.

        Parameters
        ----------
        data  : array_like
            the data vector to be fitted
        weights : array_like
            weights pertaining to the data
        par0 : array_like
            initial values for the parameters of the model
            default: from model
        keep : dict of {int : float}
            dictionary of indices (int) of parameters to be kept at fixed value (float)
            The values of `keep` are only valid for *this* fit
            see also `LevenbergMarquardtFitter( ..., keep=dict )
        limits : None or list of 2 floats or list of 2 array_like
            None : no limits applied
            [lo,hi] : low and high limits for all values of the parameters
            [la,ha] :  arrays of low and high limits for all values of the parameters
        maxiter : int
            max number of iterations. default=1000,
        tolerance : float
            absolute and relative tolrance. default=0.0001,
        verbose : int
            0 : silent
            >0 : more output
            default=1
        plot : bool
            plot the results
        callback : callable
            is called each iteration as
            `val = callback( val )`
            where val is the parameter list.

        Raises
        ------
        ConvergenceError if it stops when the tolerance has not yet been reached.

        """
        if maxiter is None : maxiter = self.maxIter
        if tolerance is None : tolerance = self.tolerance
        if verbose is None : verbose = self.verbose

        fitIndex, data, weights = self.fitprolog( data, weights=weights, keep=keep )


        trypar = self.model.parameters if par0 is None else par0

#        if fitIndex is not None and len( fitIndex ) < len( par0 ) :
#            par0 = par0[fitIndex]

        self.chi = self.chiSquaredExtra( data, trypar, weights=weights ) + 1

        self.lamda = 0.001

        if verbose > 1 and self.first :
            print( "  Iter  loglambda    chisq   parameters" )
            self.first = False

        while self.iter < maxiter :

            trypar, trychi = self.trialfit( trypar, fitIndex, data, weights, verbose, maxiter )
            self.model.parameters = trypar

            tol = tolerance if self.chi < 1 else tolerance * self.chi
            if abs( self.chi - trychi ) < tol :
                self.chi = trychi
                self.report( verbose, trypar, trychi, more=math.log10( self.lamda ), force=True )
                self.converged = True

                self.fitpostscript( data, plot=plot )
                return trypar                   # we are done.

            self.chi = trychi
            if self.lamda > 1e-100 :
                self.lamda *= 0.1

        raise ConvergenceError( "LevenbergMarquardtFitter.fit( ): " +
                                "Too many iterations: ", self.iter )


    def chiSquaredExtra( self, data, params, weights=None ) :
        """
        Add normalizing data to chisq.
        """
        chisq = self.chiSquared( data, params=params, weights=weights )
        if hasattr( self, "normdfdp" ) :
            res = numpy.inner( self.normdfdp, params ) - self.normdata
            res *= res * self.normweight
            chisq += numpy.sum( res )

        return chisq


    #  *************************************************************************
    def trialfit( self, params, fi, data, weights, verbose, maxiter ):

        hessian = self.getHessian( params=params, weights=weights, index=fi )

        residu = data - self.model.result( self.xdata, params )
        if weights is not None :
            residu *= weights
        if hasattr( self, "normdfdp" ) :
            nres = ( self.normdata - numpy.inner( self.normdfdp, params ) ) * self.normweight
            residu = numpy.append( residu, nres )

        vector = self.getVector( residu, index=fi )

        nfit = len( fi )
        fitpar = params[fi]

        while self.iter < maxiter :

            for k in range( nfit ) :
                hessian[k,k] *= ( 1 + self.lamda )

            newpar = fitpar + 0.5 * numpy.linalg.solve( hessian, vector )

            onEdge, edgePar, edgeInd = self.checkLimits( newpar, fi )

            if onEdge :
                newpar = edgePar

            trypar = self.insertParameters( newpar, index=fi, into=params )

            trychi = self.chiSquaredExtra( data, trypar, weights=weights )

            self.ntrans += 1
            self.iter += 1

#            print( "trialfit", self.chi, trychi, math.log10( self.lamda ) )
            self.report( verbose, trypar, trychi, more=math.log10( self.lamda ),
                         force=(verbose >= 3) )

            if trychi <= self.chi:

                if onEdge :                        # further convergence on the edge plane(s)
                    self.chi = trychi
                    self.model.parameters = trypar
                    trypar, trychi = self.trialfit( trypar, edgeInd, data, weights,
                            verbose, maxiter )


                return ( trypar, trychi )          #  succesfull step

            if self.lamda < 1e20 :
                self.lamda *= 10
            self.fitpar = trypar                   #  keep to report back

        raise ConvergenceError( "LevenbergMarquardtFitter. Too many iterations: ", self.iter )


    ### TBD to __getattr__  ???
    def getParameters( self ):
        """
        Return status of the fitter: parameters ( for debugging ).

        Only for debugging; use Model.getParameters( ) otherwise

        """
        return self.model.parameters if self.converged else self.fitpar

    def __str__( self ):
        """ Return name of the fitter.  """
        return "LevenbergMarquardtFitter"

    def checkLimits(  self, fitpar, fitindex ):
        """

        """
        onedge = False
#        if self.model.priors is not None and self.model.priors.hasLimits( ):
        if self.model.priors is not None :
            fitin = []
            for i,k in enumerate( fitindex ) :
                pr = self.model.getPrior( k )
                if pr.hasLowLimit( ) and fitpar[i] < pr.lowLimit :
                    fitpar[i] = pr.lowLimit
                    onedge = True
                elif pr.hasHighLimit( ) and fitpar[i] > pr.highLimit :
                    fitpar[i] = pr.highLimit
                    onedge = True
                else:
                    fitin += [k]

            return ( onedge, fitpar, fitin )
        else:
            return ( onedge, fitpar, fitindex )




