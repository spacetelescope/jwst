import numpy as numpy
import math
import warnings

from .ErrorDistribution import ErrorDistribution
from .LogFactorial import logFactorial

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
#  *    2010 - 2014 Do Kester, SRON (Java code)
#  *    2017 - 2020 Do Kester


class PoissonErrorDistribution( ErrorDistribution ):
    """
    To calculate a Poisson likelihood.

    For one observation with n counts it holds

        f( n,x ) = x^n / ( e^x * n! )

    where x is the expected counts

    The function is mostly used to calculate the likelihood L, or easier
    to use log likelihood, logL.

        logL = &sum;( n * log( x ) - x - log( n! ) )

    Weights are not accepted in this ErrorDistribution; they are silently ignored.


    Author       Do Kester.

    """

    PARNAMES = []

    #  *********CONSTRUCTORS***************************************************
    def __init__( self, copy=None ):
        """
        Constructor.

        Parameters
        ----------
        copy : PoissonErrorDistribution
            distribution to be copied.

        """
        super( PoissonErrorDistribution, self ).__init__( copy=copy )

    def copy( self ):
        """ Return copy of this.  """
        return PoissonErrorDistribution( copy=self )

    def acceptWeight( self ):
        """
        True if the distribution accepts weights.
        Always false for this distribution.
        """
        return False


    def getScale( self, problem, allpars=None ) :
        """
        Return the noise scale.

        *** Gaussian approximation ***

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            None take parameters from problem.model
            list of all parameters in the problem
        """
        return self.getGaussianScale( problem, allpars=allpars )

    #  *********LIKELIHOODS***************************************************
    def logLikelihood_alt( self, problem, allpars ):
        """
        Return the log( likelihood ) for a Poisson distribution.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            list of all parameters in the problem

        """
        self.ncalls += 1

        mock = problem.result( allpars )
        if numpy.any(  mock <= 0.0 ) :
            return -math.inf

        lfdata = logFactorial( problem.ydata )

        logl = numpy.sum( problem.ydata * numpy.log( mock ) - mock - lfdata )

        if math.isnan( logl ) :
            return -math.inf

        return logl

    def logLdata( self, problem, allpars, mockdata=None ) :
        """
        Return the log( likelihood ) for each residual

        logL = sum( logLdata )

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            list of all parameters in the problem
        mockdata : array_like
            as calculated by the model

        """
        if mockdata is None :
            mockdata = problem.result( allpars )
        lfdata = logFactorial( problem.ydata )

        with warnings.catch_warnings():
            warnings.simplefilter( "ignore", category=RuntimeWarning )
            lld = problem.ydata * numpy.log( mockdata ) - mockdata - lfdata

        lld = numpy.where( numpy.isfinite( lld ), lld, -math.inf )

#        lld = numpy.where( mockdata <= 0, -math.inf,
#                problem.ydata * numpy.log( mockdata ) - mockdata - lfdata )


        return lld

    def partialLogL_alt( self, problem, allpars, fitIndex ):
        """
        Return the partial derivative of log( likelihood ) to the parameters.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            list of all parameters in the problem
        fitIndex : array_like
            indices of parameters to be fitted
        """
        self.nparts += 1
        mock = problem.result( allpars )
        dM = problem.partial( allpars )
        dL = numpy.zeros( len( fitIndex ), dtype=float )

        i = 0
        for k in fitIndex :
            dL[i] = numpy.sum( ( problem.ydata / mock - 1 ) * dM[:,k] )
            i += 1

        return dL

    def nextPartialData( self, problem, allpars, fitIndex, mockdata=None ):
        """
        Return the partial derivative of log( likelihood ) to the parameters.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            list of all parameters in the problem
        fitIndex : array_like
            indices of parameters to be fitted
        mockdata : array_like
            as calculated by the model
        """
        if mockdata is None :
            mockdata = problem.result( allpars )
        dM = problem.partial( allpars )
##      TBD import mockdata into partial
#        dM = problem.partial( allpars, mockdata=mockdata )

        for k in fitIndex :
            yield ( problem.ydata / mockdata - 1 ) * dM[:,k]

        return

    def __str__( self ) :
        return "Poisson error distribution"

