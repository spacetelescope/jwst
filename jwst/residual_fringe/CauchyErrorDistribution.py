import numpy as numpy
import scipy
import math

from .ScaledErrorDistribution import ScaledErrorDistribution
from .NoiseScale import NoiseScale

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


class CauchyErrorDistribution( ScaledErrorDistribution ):
    """
    To calculate a Cauchy or Lorentz likelihood.

        f( x ) = s / ( &pi; * ( s^2 + x^2 ) )

    where x = residual and s = scale

    The function is mostly used to calculate the likelihood L, or easier
    to use log likelihood, logL.

        logL = N ( log( s ) - log( &pi; ) ) - sum( log( x^2 + s^2 ) )

    Weights are not possible in this error distribution. They are silently ignored.

    s is a hyperparameter, which might be estimated from the data.

    Author       Do Kester.

    """

    LOGPI = math.log( math.pi )

    #  *********CONSTRUCTORS***************************************************
    def __init__( self, scale=1.0, limits=None, copy=None ):
        """
        Constructor.

        Parameters
        ----------
        scale : float
            noise scale
        limits : None or list of 2 floats [low,high]
            None : no limits implying fixed scale
            low     low limit on scale (needs to be >0)
            high    high limit on scale
            when limits are set, the scale is not fixed.
         copy : CauchyErrorDistribution
            distribution to be copied.

        """

        super( CauchyErrorDistribution, self ).__init__( scale=scale, limits=limits,
            copy=copy )

    def copy( self ):
        """ Return copy of this.  """
        return CauchyErrorDistribution( copy=self )

    def __setattr__( self, name, value ):
        """
        Set attributes. Needed for getScale()

        """
        if name == "res2" or name == "sumweight":
            object.__setattr__( self, name, value )
        else :
            super( CauchyErrorDistribution, self ).__setattr__( name, value )

    def acceptWeight( self ):
        """
        True if the distribution accepts weights.
        False for this distribution.
        """
        return False

    def getScale( self, problem, allpars=None ) :
        """
        Return the noise scale as calculated from the residuals.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            None take parameters from problem.model
            list of all parameters in the problem
        """
        res = self.getResiduals( problem, allpars=allpars )
        self.res2 = res * res
        self.sumweight = problem.sumweight

        scale = scipy.optimize.bisect( self.funct, 0.001, 100.0 )
        return scale

    def funct( self, scale ) :
        """
        Internal use, only.
        """
        return ( numpy.sum( numpy.log( self.res2 + scale * scale ) ) -
                self.sumweight * ( 2 * math.log( scale ) + math.sqrt( 2.0 ) ) )

    #  *********LIKELIHOODS***************************************************
    def logLikelihood_alt( self, problem, allpars ):
        """
        Return the log( likelihood ) for a Cauchy distribution.
        Cauchy distr : f( x ) = s / ( pi * ( s^2 + x^2 ) )

        where x = residual and s = scale

        Alternate calculation

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            list of all parameters in the problem

        """
        self.ncalls += 1

        scale = allpars[-1]
        res = problem.residuals( allpars[:-1] )
        res2 = numpy.square( res )
        return ( problem.ndata * ( math.log( scale ) - self.LOGPI ) -
                 numpy.sum( numpy.log( res2 + scale * scale ) ) )

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
        res = problem.residuals( allpars[:-1], mockdata=mockdata )

        scale = allpars[-1]
        s2 = scale * scale
        res2 = res * res
        return math.log( scale ) - self.LOGPI - numpy.log( res2 + s2 )

    def partialLogL_alt( self, problem, allpars, fitIndex ) :
        """
        Return the partial derivative of log( likelihood ) to the parameters
        in fitIndex.

        Alternate calculation

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            parameters of the problem
        fitIndex : array_like
            indices of parameters to be fitted

        """
        self.nparts += 1

        scale = allpars[-1]
        res = problem.residuals( allpars[:-1] )
        r2s = res * res + scale * scale
        dM = problem.partial( allpars[:-1] )

        dL = numpy.zeros( len( fitIndex ), dtype=float )
        i = 0
        for k in fitIndex :
            if k >= 0 :
                dL[i] = 2 * numpy.sum( res * dM[:,k] / r2s )
                i += 1
            else :
                dL[-1] = self.ndata / scale - numpy.sum( 2 * scale / r2s )

        return dL

    def nextPartialData( self, problem, allpars, fitIndex, mockdata=None ) :
        """
        Return the partial derivative of log( likelihood ) to the parameters
        in fitIndex.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            parameters of the problem
        fitIndex : array_like
            indices of parameters to be fitted
        mockdata : array_like
            as calculated by the model

        """
        res = problem.residuals( allpars[:-1], mockdata=mockdata )
        scale = allpars[-1]
        r2s = res * res + scale * scale

        dM = problem.partial( allpars[:-1] )
##      TBD import mockdata into partial
#        dM = model.partial( self.xdata, param, mockdata=mockdata )

        for k in fitIndex :
            if k >= 0 :
                yield 2 * res * dM[:,k] / r2s
            else :
                yield 1.0 / scale - 2 * scale / r2s

        return

    def __str__( self ) :
        return "Cauchy error distribution"

