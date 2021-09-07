import numpy as numpy
from scipy import special
import math

from .ScaledErrorDistribution import ScaledErrorDistribution
from .HyperParameter import HyperParameter
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
#  *    2017 - 2020 Do Kester


class ExponentialErrorDistribution( ScaledErrorDistribution ):
    """
    Also known as generalized gaussian errordistribution.

    To calculate an Exponential likelihood.

    For one residual, x, it holds

        f( x ) = p / ( 2 s &Gamma;( 1 / p ) ) exp( - ( |x| / s ) ^ p )

    where s is the scale and p is the power.
    s and p are hyperparameters, which might be estimated from the data.

    The variance of this function is

        &sigma; ^ 2 = s ^ 2 &Gamma;( 3 / p ) / &Gamma;( 1 / p )

    See toSigma()

    The function is mostly used to calculate the likelihood L over N residuals,
    or easier to use log( L )

        logL = log( N p / ( 2 s &Gamma;( 1 / p ) ) ) - &sum;( ( |x| / s ) ^ p )

    Using weights this becomes:

        logL = log( &sum;( w ) p / ( 2 s &Gamma;( 1 / p ) ) ) - &sum;( w ( |x| / s ) ^ p )

    Note
    ----
    The scale s in Exponential is NOT the same as the scale in Gaussian or in Laplace.

    Attributes from ErrorDistibution
    --------------------------------
    hyperpar, deltaP, ncalls, nparts, sumweight, ndata, hypar, nphypar


    Author       Do Kester.

    """
    LOG2PI = math.log( 2 * math.pi )
    PARNAMES = ["scale", "power"]

    #  *********CONSTRUCTORS***************************************************
    def __init__( self, scale=1.0, power=2.0, limits=None, copy=None ):
        """
        Default Constructor.

        Parameters
        ----------
        scale : float
            noise scale
        power : float
            power of the distribution
        limits : None or [low,high] or [[low],[high]]
            None : no limits implying fixed scale
            low     low limit on scale (needs to be >0)
            high    high limit on scale
            [low]   low limit on [scale,power] (need to be >0)
            [high]  high limit on [scale,power]
            when limits are set, the scale cq. power are *not* fixed.
        copy : ExponentialErrorDistribution
            distribution to be copied.
        """
        super( ExponentialErrorDistribution, self ).__init__( limits=None, copy=copy )

        plim = None
        if limits is None :
            slim = None
        else :
            lo = limits[0]
            hi = limits[1]
            try :
                slim = [lo[0],hi[0]]
            except :
                slim = [lo,hi]
            try :
                plim = [lo[1],hi[1]]
            except :
                pass

        if copy is None :
            self.hyperpar = [NoiseScale( scale=scale, limits=slim ),
                             HyperParameter( hypar=power, limits=plim )]
        else :
            self.hyperpar = copy.hyperpar

    def copy( self ):
        """ Return copy of this.  """
        return ExponentialErrorDistribution( copy=self )

    def acceptWeight( self ):
        """
        True if the distribution accepts weights.
        Always true for this distribution.
        """
        return True

    def toSigma( self, hypar ) :
        """
        Return sigma, the squareroot of the variance.
        Parameter
        --------
        hypar : array_like (2 floats)
            the [scale,power] of this Exponential distribution.
        """
        p = hypar[1]
        return hypar[0] * math.sqrt( special.gamma( 3.0 / p ) / special.gamma( 1.0 / p ) )

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

        scale = allpars[-2]
        power = allpars[-1]

        chipow = self.getChipow( problem, allpars ) / math.pow( scale, power )
        norm = math.log( 0.5 * power / scale ) - special.gammaln( 1.0 / power )

        return self.sumweight * norm - chipow

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
        scale = allpars[-2]
        power = allpars[-1]

        res = problem.residuals( allpars[:-2], mockdata=mockdata )

        lld = - numpy.power( numpy.abs( res / scale ), power )
        norm = math.log( power / ( 2 * scale ) ) - special.gammaln( 1.0 / power )
        lld += norm
        if problem.weights is not None :
            lld *= problem.weights
        return lld

    def getChipow( self, problem, allpars=None ) :
        """
        Return chisq.

        return Sum over the (weighted) powered residuals

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            None take parameters from problem.model
            list of all parameters in the problem

        """
        res = problem.residuals( allpars[:-2] )
        power = allpars[-1]

        ares = numpy.power( numpy.abs( res ), power )
        if problem.weights is not None :
            ares = ares * problem.weights
        return numpy.sum( ares )

    def getScale( self, problem, allpars=None ) :
        """
        Return the noise scale calculated from the residuals.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            None take parameters from problem.model
            list of all parameters in the problem
        """
        power = allpars[-1]
        chi = self.getChipow( problem, allpars=allpars )
        return math.pow( chi / problem.sumweight, 1.0 / power )


    def partialLogL_alt( self, problem, allpars, fitIndex ) :
        """
        Return the partial derivative of log( likelihood ) to the parameters.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            parameters of the problem
        fitIndex : array_like
            indices of parameters to be fitted

        """
        self.ncalls += 1

        scale = allpars[-2]
        power = allpars[-1]
        res = problem.residuals( allpars[:-2] )

        ars = numpy.abs( res / scale )
        rsp = numpy.power( ars, power )
        if problem.weights is not None :
            rsp = rsp * problem.weights

        dLdm = power * rsp / res
        dM = problem.partial( allpars[:-2] )

        dL = numpy.zeros( len( fitIndex ), dtype=float )
        i = 0
        for  k in fitIndex :
            if k >= 0 :
                dL[i] = numpy.sum( dLdm * dM[:,k] )
                i += 1
            elif k == -2 :
                dL[-2] = - problem.sumweight / scale + power * numpy.sum( rsp ) / scale
            else :
                # special.psi( x ) is the same as special.polygamma( 1, x )
                dldp = problem.sumweight * ( power + special.psi( 1.0 / power ) )
                dldp /= ( power * power )
                dldp -= ( numpy.sum( rsp * numpy.log( ars ) ) )
                dL[-1] = dldp

        return dL

    def nextPartialData( self, problem, allpars, fitIndex, mockdata=None ) :
        """
        Return the partial derivative of all elements of the log( likelihood )
        to the parameters.

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
        param = allpars[:-2]
        res = problem.residuals( param, mockdata=mockdata )
        scale = allpars[-2]
        power = allpars[-1]

        ars = numpy.abs( res / scale )
        rsp = numpy.power( ars, power )
        if problem.weights is not None :
            rsp = rsp * problem.weights
            wgt = problem.weights
        else :
            wgt = 1.0

        dLdm = power * rsp / res
        dM = problem.partial( param )
##      TBD import mockdata into partial
#        dM = problem.partial( param, mockdata=mockdata )

        # special.psi( x ) is the same as special.polygamma( 1, x )
        dlp = wgt * ( power + special.psi( 1.0 / power ) ) / ( power * power )

        for  k in fitIndex :
            if k >= 0 :
                yield ( dLdm * dM[:,k] )
            elif k == -2 :
                yield ( power * rsp - wgt ) / scale
            else :
                yield dlp - rsp * numpy.log( ars )

        return

    def __str__( self ) :
        return "Exponential error distribution"

