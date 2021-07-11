import numpy as numpy
import math

from .HyperParameter import HyperParameter
from . import Tools
from .Tools import setAttribute as setatt

__author__ = "Do Kester"
__year__ = 2021
__license__ = "GPL3"
__version__ = "2.7.0"
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
#  *    2017 - 2021 Do Kester


class ErrorDistribution( object ):
    """
    ErrorDistribution defines general methods for a error distribution.

    Error distributions are used to calculate the likelihoods.

    Author       Do Kester.

    Attributes
    ----------
    hyperpar : HyperParameter
        hyperparameter for the error distribution
    deltaP : float
        delta for calculating numerical derivatives
    ncalls : int
        number of calls to the logLikelihood
    nparts : int
        number of calls to the partial of the logLikelihood
    sumweight : float
        sum over the weights or ndata
    ndata : int
        number of points in data
    hypar : [float]
        list of values for the hyperparameters
    nphypar : int
        number of hyper parameters in this error distribution
    constrain : None or callable
        None:     Use logLikelihood as is
        callable: logL = func( logL, problem, allpars, lowLhood )
                  returning a (modified) value of the logLikelihood.
    """

    PARNAMES = ["hypar"]

    #  *********CONSTRUCTORS***************************************************
    def __init__( self, fixed=None, constrain=None, copy=None ):
        """
        Constructor.

        Parameters
        ----------
        fixed : dictionary of {int:float}
            int     list if parameters to fix permanently. Default None.
            float   list of values for the fixed parameters.
        constrain : None or callable
            function as: func( logL, problem, allpars )
            returning a (modified) value of the logLikelihood.
        copy : ErrorDistribution
            distribution to be copied.

        Raise
        -----
        ValueError when constrain is not a callable method.

        """
        super( ErrorDistribution, self ).__init__()
        self.ncalls = 0
        self.nparts = 0
        if copy is None :
            self.hyperpar = []
            self.deltaP = 0.000001
            self.fixed = self.keepFixed( fixed )
            self.constrain = constrain
        else :                          ### TBC do we need copies here ????
            self.hyperpar = copy.hyperpar
            self.deltaP = copy.deltaP
            self.fixed = self.keepFixed( copy.fixed )
            self.constrain = copy.constrain

    def copy( self ):
        """ Return copy of this.  """
        return ErrorDistribution( copy=self )

    def __setattr__( self, name, value ):
        """
        Set attributes.

        """
        if name == "constrain" :
            if hasattr( self, "logLdata" ) :
                ## Choose one of 2 likelihood functions in "logLikelihood"
                if value is None :
                    setatt( self, "logLikelihood", self.logLhood )
                elif callable( value ) :
                    setatt( self, "logLikelihood", self.logCLhood )
                else :
                    raise ValueError( "Constrain is not a callable function" )

            setatt( self, "constrain", value )
            return

        key0 = ["hyperpar", "fixed"]
        keys = {"hyperpar": HyperParameter}
        key1 = {"deltaP": float, "lowLhood": float, "ncalls": int, "nparts": int, "fixed": dict }
        if ( Tools.setNoneAttributes( self, name, value, key0 ) or
             Tools.setListOfAttributes( self, name, value, keys ) or
             Tools.setSingleAttributes( self, name, value, key1 ) ) :
            pass
        else :
            raise AttributeError(
                "Object has no attribute " + name + " of type " + str( value.__class__ ) )

    def __getattr__( self, name ) :
        """
        Return value belonging to attribute with name.

        Parameters
        ----------
        name : string
            name of the attribute
        """
        if name == 'nphypar' :
            return len( self.hyperpar )
        elif name == "hypar" :
            return numpy.asarray( [s.hypar for s in self.hyperpar], dtype=float )
        else :
            raise AttributeError( "Unknown attribute " + name )

        return None

    def getGaussianScale( self, problem, allpars=None ) :
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
        chi = self.getChisq( problem, allpars=allpars )
        return math.sqrt( chi / problem.sumweight )


    def getResiduals( self, problem, allpars=None ):
        """
        Return residuals: ydata - model.result

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            None take parameters from problem.model
            list of all parameters in the problem
        """
        if allpars is None :
            param = problem.model.parameters
        else :
            nh = self.nphypar
            param = allpars[:-nh]

        return problem.residuals( param )


    def getChisq( self, problem, allpars=None ):
        """
        Return chisq

        *** Gaussian approximation ***

        Sum over the (weighted) squared residuals

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            None take parameters from problem.model
            list of all parameters in the problem
        """
        if allpars is None :
            param = problem.model.parameters
        else :
            nh = self.nphypar
            param = allpars[:-nh]

        res2 = problem.weightedResSq( param )

        return numpy.sum( res2 )

    def toSigma( self, scale ) :
        """
        Return sigma, the squareroot of the variance.

        Parameter
        --------
        scale : float
            the scale of this distribution.

        Return default value : scale
        """
        return scale

    def isBound( self ) :
        """
        True when all priors of its (hyper)parameters are bound
        """
        bound = True
        if self.hyperpar is None :
            return bound

        for sp in self.hyperpar :
            bound = bound and sp.prior.isBound()

        return bound

    def acceptWeight( self ):
        """ True if the distribution accepts weights.  """
        return False

    def keepFixed( self, fixed=None ) :
        """
        Keeps (hyper)parameters fixed at the provided values.

        1. Repeated calls start from scratch.<br>
        2. Reset with keepFixed( fixed=None )

        Parameters
        ----------
        fixed : dictionary of {int:float}
            int     list if parameters to fix permanently. Default None.
            float   list of values for the fixed parameters.
        """
        if self.hyperpar is None :
            return None

        # Reset all fixed to False
        for s in self.hyperpar :
            s.fixed = False

        if fixed is None :
            return None

        for f in fixed.keys() :
            self.hyperpar[f].fixed = True
            self.hyperpar[f].hypar = fixed[f]

        return fixed


    def setPriors( self, priors ) :
        """
        Set priors on the hyper parameter(s).

        Parameters
        ----------
        priors  : (list of) Prior
            prior distribution for the hyperparameters
        """
        for spp,pr in zip( self.hyperpar, priors ) :
            spp.setPrior( pr )

    def setLimits( self, limits ) :
        """
        Set limits on the hyper parameter(s).

        Parameters
        ----------
        limits : [low,high]
            low : float or array_like
                low limits
            high : float or array_like
                high limits
        """
        if self.nphypar == 1 :
            self.hyperpar[0].setLimits( limits )
            return
        low = limits[0]
        high = limits[1]
        for spp,lo,hi in zip( self.hyperpar, low, high ) :
            spp.setLimits( [lo,hi] )

    def domain2Unit( self, dval, ks ) :
        """
        Return value in [0,1] for the selected parameter.
        Parameters
        ----------
        dval : float
            hyper parameter value in domain
        ks : int
            selecting index
        """
        try :
            return self.hyperpar[ks].domain2Unit( dval )
        except :
            return 0.0

    def unit2Domain( self, uval, ks ) :
        """
        Return domain value for the selected parameter.
        Parameters
        ----------
        uval : float
            unit value of hyper parameter
        ks : int
            selecting index
        """
        try :
            return self.hyperpar[ks].unit2Domain( uval )
        except :
            return 0.0

    #  *********LIKELIHOODS***************************************************

    def logCLhood( self, problem, allpars ):
        """
        Return the constrained log( likelihood ).

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            parameters of the problem
        """
        logL = self.logLhood( problem, allpars )

        return self.constrain( logL, problem, allpars, self.lowLhood )

    def logLhood( self, problem, allpars ):
        """
        Return the log( likelihood ).

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            parameters of the problem
        """
        self.ncalls += 1

        return numpy.sum( self.logLdata( problem, allpars ) )


    def partialLogL( self, problem, allpars, fitIndex ) :
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
        self.nparts += 1                ## counts calls tp partialLogL
        mock = problem.result( allpars[:problem.model.npars] )

        try :
            pg = self.nextPartialData( problem, allpars, fitIndex, mockdata=mock )
            np = len( fitIndex )
            dL = numpy.zeros( np, dtype=float )
            for k in range( np ) :
                dL[k] = numpy.sum( next( pg ) )

            try :
                pg.close()
            except :
                pass

            return dL
        except :
            print( "Using numeric partialLogL." )
            return self.numPartialLogL( problem, allpars, fitIndex )


    def partialLogL_alt( self, problem, allpars, fitIndex ) :
        """
        Return the partial derivative of log( likelihood ) to the parameters.

        Alternative calculation.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            parameters of the problem
        fitIndex : array_like
            indices of parameters to be fitted

        """
        return self.numPartialLogL( problem, allpars, fitIndex )

    def numPartialLogL( self, problem, allpars, fitIndex ) :
        """
        Return d log( likelihood ) / dp, numerically calculated.

        Parameters
        ----------
        problem : Problem
            to be solved
        allpars : array_like
            parameters of the problem
        fitIndex : array_like
            indices of parameters to be fitted

        """
        dL = numpy.zeros( len( fitIndex ), dtype=float )
        p = allpars.copy( )
        i = 0
        for k in fitIndex :
            p[k] = allpars[k] - self.deltaP
            lm = self.logLikelihood( problem, p )
            p[k] = allpars[k] + self.deltaP
            lp = self.logLikelihood( problem, p )
            dL[i] = 0.5 * ( lp - lm ) / self.deltaP
            i += 1
            p[k]= allpars[k]
        return dL


    def updateLogL( self, problem, allpars, parval=None ):
        """
        Return a update of the log( likelihood ) given a change in a few parameter.

        This method provides the opportunity to optimize the logL calculation.
        Providing this one, automatically provides the previous one.
        For now it just refers to logLikelihood() itself.

        Parameters
        ----------
        problem : Problem
            to be solved
        param : array_like
            parameters of the model
        parval : dict of {int : float}
            int index of a parameter
            float (old) value of the parameter
        """
        return self.logLikelihood( problem, allpars )

    def setResult( self ):
        pass

    def __str__( self ) :
        """
        Return a string representation of this class.
        """
        return "Error distribution"

    def hyparname( self, k ) :
        """
        Return name of the hyperparameter

        Parameters
        ----------
        k : int
            index of the hyperparameter
        """
        return self.PARNAMES[k]



