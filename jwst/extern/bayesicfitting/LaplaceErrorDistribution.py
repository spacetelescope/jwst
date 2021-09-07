import numpy as numpy
import math

from .Formatter import formatter as fmt
from .ScaledErrorDistribution import ScaledErrorDistribution

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


class LaplaceErrorDistribution( ScaledErrorDistribution ):
    """
    To calculate a Laplace likelihood.

    For one residual, x, it holds

        f( x ) = 1 / ( 2 s ) exp( - |x| / s )

    where s is the scale.
    s is a hyperparameter, which might be estimated from the data.

    The variance of this function is &sigma;^2 = 2 s ^ 2.
    See: toSigma()

    The function is mostly used to calculate the likelihood L over N
    residuals, or easier using log likelihood, logL.

        logL = log( N / ( 2 s ) ) - &sum;( |x| / s  )

    Using weights this becomes:

        logL = log( &sum;( w ) / ( 2 s ) ) - &sum;( w |x| / s  )

    Using this error distribution results in median-like solutions.

    Author       Do Kester.

    """
    SQRT2 = math.sqrt( 2 )
    LGSQ2 = math.log( SQRT2 )
    LOG2 = math.log( 2.0 )


    #  *********CONSTRUCTORS***************************************************
    def __init__( self, scale=1.0, limits=None, copy=None ) :
        """
        Constructor of Laplace Distribution.

        Parameters
        ----------
        scale : float
            noise scale
        limits : None or list of 2 floats [low,high]
            None : no limits implying fixed scale
            low     low limit on scale (needs to be >0)
            high    high limit on scale
            when limits are set, the scale is *not* fixed.

        copy : LaplaceErrorDistribution
            distribution to be copied.
        """
        super( LaplaceErrorDistribution, self ).__init__( scale=scale,
                limits=limits, copy=copy )

    def copy( self ):
        """ Return copy of this.  """
        return LaplaceErrorDistribution( copy=self )

    #  *********DATA & WEIGHT***************************************************
    def acceptWeight( self ):
        """
        True if the distribution accepts weights.
        Always true for this distribution.
        """
        return True

    def toSigma( self, scale ) :
        """
        Return sigma, the squareroot of the variance.
        Parameter
        --------
        scale : float
            the scale of this Laplace distribution.
        """
        return scale * math.sqrt( 2.0 )

    def getScale( self, problem, allpars=None ) :
        """
        Return the noise scale

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            None take parameters from problem.model
            list of all parameters in the problem
        """
        sumres = self.getSumRes( problem, allpars=allpars )
        return sumres / problem.sumweight


    def getSumRes( self, problem, allpars=None ):
        """
        Return the sum of the absolute values of the residuals.

        ..math ::
            \sum ( | res | )

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            None take parameters from problem.model
            list of all parameters in the problem
        """
        res = self.getResiduals( problem, allpars=allpars )

        if problem.weights is not None :
            res *= problem.weights
        return numpy.sum( numpy.abs( res ) )


    #  *********LIKELIHOODS***************************************************
    def logLikelihood_alt( self, problem, allpars ) :
        """
        Return the log( likelihood ) for a Gaussian distribution.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            parameters of the problem

        """
        self.ncalls += 1

        scale = allpars[-1]
        sumres = self.getSumRes( problem, allpars ) / scale
        return - problem.sumweight * ( self.LOG2 + math.log( scale ) ) - sumres

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
        np = problem.npars
        res = problem.residuals( allpars[:np], mockdata=mockdata )

        scale = allpars[-1]
        res = - numpy.abs( res ) / scale - ( self.LOG2 + math.log( scale ) )
        if problem.weights is not None :
            res = res * problem.weights
        return res

    def partialLogL_alt( self, problem, allpars, fitIndex ) :
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
        scale = allpars[-1]

        dM = problem.partial( allpars[:-1] )
        res = problem.residuals( allpars[:-1] )
        wgt = numpy.ones_like( res, dtype=float ) if problem.weights is None else problem.weights
        wgt = numpy.copysign( wgt, res )

        dL = numpy.zeros( len( fitIndex ), dtype=float )
        i = 0
        for k in fitIndex :
            if k >= 0 :
                dL[i] = numpy.sum( wgt * dM[:,k] )
                i += 1
            else :
                dL[-1] = self.getSumRes( problem, allpars ) / scale - problem.sumweight
        return dL / scale

    def nextPartialData( self, problem, allpars, fitIndex, mockdata=None ) :
        """
        Return the partial derivative of elements of the log( likelihood )
        to the parameters.

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
        param = allpars[:-1]
        scale = allpars[-1]
        res = problem.residuals( param, mockdata=mockdata )

        dM = problem.partial( param )
##      TBD import mockdata into partial
#        dM = model.partial( self.xdata, param, mockdata=mockdata )

        wgt = numpy.ones_like( res, dtype=float ) if problem.weights is None else problem.weights
        swgt = numpy.copysign( wgt, res )

        res *= swgt / scale             ## make all residuals >= 0

        for k in fitIndex :
            if k >= 0 :
                yield ( swgt * dM[:,k] ) / scale
            else :
                yield ( res - wgt ) / scale
        return

    def __str__( self ) :
        return "Laplace error distribution"

