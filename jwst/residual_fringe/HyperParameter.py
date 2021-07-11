import numpy as numpy
from astropy import units
import math
from . import Tools
from .Prior import Prior

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
#  *    2016 - 2020 Do Kester

class HyperParameter( object ):
    """
    Values and priors for the parameter(s) of an ErrorDistribution.

    Hyperparameters are not directly related to the model, they are
    parameters of the error distribution.

    Information about the scale of the noise is stored in a derived class,
    noiseScale.

    The full use of priors is reserved for Bayesian calculations as
    in NestedSampler

    Attributes
    ----------
    hypar : float
        the value of the hyperparameter.  Default: 1.0
    stdev : float
        the standard deviation of the hyperparameter.  Default: None
    prior : Prior
        the prior for the hyperparameter.
    isFixed : boolean
        keep the hyperparameter fixed at the value given by hypar.
        default: True

    """

    #  *********CONSTRUCTORS***************************************************
    def __init__( self, hypar=1, isFixed=True, prior=None, limits=None,
                    copy=None ):
        """
        Constructor.

        Parameters
        ----------
        hypar : float
            value of the hyperparameter
        isFixed : bool
            True:   Consider the hyperparameter as fixed
            False:  Optimize the parameter too (when relevant)
                    It might need a prior and/or limits to be set
        prior : None or Prior
            None : no prior set
            prior probability on the hyperparameter
        limits : None or list of 2 floats [lo,hi]
            low limit and high limit on hypar.
        copy : HyperParameter
            HyperParameter to copy

        """
        super( HyperParameter, self ).__init__( )
        self.hypar = hypar
        self.isFixed = isFixed
        self.stdev = None
        self.prior = prior
        if prior is not None :
            self.prior.setLimits( limits )
            self.isFixed = False

        if copy is not None :
            self.hypar = copy.hypar
            self.stdev = copy.stdev
            if copy.prior is not None :
                self.prior = copy.prior.copy()
            self.isFixed = copy.isFixed

    def copy( self ):
        """ Return a copy.  """
        return HyperParameter( copy=self )

    def checkPrior( self ) :
        """
        Raises
        ------
        ValueError when no prior has been set.
        """
        if self.prior is None :
            raise ValueError( "Need a prior set to %s" % self.__str__() )

    def setLimits( self, limits ):
        """
        Set the limits on the scale within the prior.

        Parameters
        ----------
        limits : list of 2 float
            the [low,high] limits.

        """
        self.checkPrior()
        self.prior.setLimits( limits )
        self.isFixed = False

    def getLimits( self ):
        """ Return the limits on the scale. """
        return self.prior.getLimits( )

    def isBound( self ) :
        """ Return true is the itergral over the prior is bound. """
        self.checkPrior()
        return self.prior.isBound()

    def domain2Unit( self, dval ):
        """
        Return a value in [0,1] given a value within the valid domain of
        a parameter for the prior distribution.

        Parameters
        ----------
        dval : float
            value within the domain of a parameter

        """
        self.checkPrior()
        return self.prior.domain2Unit( dval )

    def unit2Domain( self, uval ):
        """
        Return a value within the valid domain of the parameter given a value
        between [0,1] for the prior distribution.

        Parameters
        ----------
        uval : float
            value within [0,1]

        """
        self.checkPrior()
        return self.prior.unit2Domain( uval )

    def partialDomain2Unit( self, dval ):
        """
        Return a the derivate of the domain2Unit function to dval.

        Parameters
        ----------
        dval : float
            value within the domain of a parameter

        """
        self.checkPrior()
        return self.prior.partialDomain2Unit( dval )

    def __str__( self ) :
        return str( "HyperParameter. value = %f" % self.hypar )
