import numpy as numpy
import math

from .ConvergenceError import ConvergenceError
from .BaseFitter import BaseFitter
from .IterationPlotter import IterationPlotter
from .GaussErrorDistribution import GaussErrorDistribution
from .LaplaceErrorDistribution import LaplaceErrorDistribution
from .CauchyErrorDistribution import CauchyErrorDistribution
from .PoissonErrorDistribution import PoissonErrorDistribution
from .ExponentialErrorDistribution import ExponentialErrorDistribution
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
#  *    2016 - 2020 Do Kester

class IterativeFitter( BaseFitter ):
    """
    Base class with methods common to all iterative fitters.

    Author:      Do Kester.

    Attributes
    ----------
    tolerance : float
        When absolute and relative steps in subsequent chisq steps are less than
        tolerance, the fitter stops. Default = 0.0001
    maxIter : int
        When the number of iterations gets larger than maxiter the fitter
        stops with a ConvergenceError Default = 1000 * nparams
    iter : int (read only)
        iteration counter
    ntrans : int (read only)
        number of transforms
    verbose : int
        information per iteration.
        0 : silent
        1 : base information (default)
        2 : report about every 100th iteration
        3 : report about every ietration
    tooLarge : int
        When the length parameter array is too large to make a Hessian.
        To avert OutOfMemory. Default = 100

    plotter : Plotter
        Iteration plotter class. Default = IterationPlotter
    plotIter : int
        Produce a plot for every plotIter-th iteration.
        Default = 0 (no plotting)

    Raises
    ------
    ConvergenceError    Something went wrong during the convergence if the fit.

    """
    def __init__( self, xdata, model, maxIter=None, tolerance=0.0001, verbose=1, **kwargs ) :

        """
        Create a new iterative fitter, providing xdatas and model.

        This is a base class. It collects stuff common to all iterative fitters.
        It does not work by itself.

        Parameters
        ----------
        xdata : array_like
            array of independent input values
        model : Model
            the model function to be fitted

        tolerance : float
            When absolute and relative steps in subsequent chisq steps are less than
            tolerance, the fitter stops. Default = 0.01
        maxIter : None or int
            When the number of iterations gets larger than maxiter the fitter
            stops with a ConvergenceError Default = 1000 * nparams
        verbose : int
            0 : silent
            1 : report result
            2 : report every 100th iteration
            3 : report every iteration
        kwargs for @BaseFitter
            map, keep, fixedScale

        """
        super( IterativeFitter, self ).__init__( xdata, model, **kwargs )

        self.tolerance = tolerance
        if maxIter is None :
            maxIter = 1000 * model.npchain
        self.maxIter = maxIter
        self.verbose = verbose
        self.tooLarge = 100
        self.plotter = IterationPlotter()

        self.plotfreq = 0

    #  *************************************************************************
    def setParameters( self, params ):
        """
        Initialize the parameters of the model
        A little superfluous: see {@link Model#setParameters}

        Parameters
        ----------
        params : array_like
            initial parameters

        """
        self.model.parameters = params

    def doPlot( self, param, force=False ):
        """
        Plot intermediate result.

        Parameters
        ----------
        param : array_like
            of the model
        force : bool
            do the plot

        """
        if ( self.plotfreq > 0 and self.iter % self.plotfreq == 0 ) or force :
            y = self.model.result( self.xdata, param )
            self.plotter.plotResult( self.xdata, y, self._iter )

    def fitprolog( self, ydata, weights=None, keep=None ) :
        """
        Prolog for all iterative Fitters.

        1. Sets up plotting (if requested)
        2. Sets self.iter and self.ntrans to 0
        3. Checks data/weighs for Nans
        4. Makes fitIndex.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted
        weights : array_like
            weights pertaining to the data
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)

        Returns
        -------
        fitIndex : ndarray of int
            Indices of the parameters that need fitting

        """
        if self.plotfreq > 0 :
            self.plotter.plotData( self.xdata, ydata, self.__str__( ) )

        self.iter = 0
        self.ntrans = 0

        return super( IterativeFitter, self ).fitprolog( ydata, weights=weights, keep=keep )

    #  *************************************************************************
    def fit( self, ydata, weights=None, keep=None, **kwargs ):
        """
        Return model parameters fitted to the data.

        It will calculate the hessian matrix and chisq.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted
        weights : array_like
            weights pertaining to the data
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)
        kwargs :
            passed to the fitter

        Raises
        ------
        ConvergenceError if it stops when the tolerance has not yet been reached.

        """
        raise ConvergenceError( "IterativeFitter is a base class, not suitable itself to perform fits." )

    def report( self, verbose, param, chi, more=None, force=False ) :
        """
        Report on intermediate results.
        """

        if verbose > 1 and ( self.iter % 100 == 0 or force ) :
            mr = "" if more is None else fmt( more, format='   %6.1f ' )
            mx = 5 if verbose < 4 else None if verbose == 4 else verbose
            print( fmt( self.iter, format='%6d' ), mr, fmt( chi, format="%8.1f " ),
                   fmt( param, max=mx ) )



    def __str__( self ):
        return "IterativeFitter"


