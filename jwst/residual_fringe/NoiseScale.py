import numpy as numpy
from astropy import units
import math
from . import Tools

from .HyperParameter import HyperParameter
from .JeffreysPrior import JeffreysPrior

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
#  *    2011 - 2014 Do Kester, SRON (JAVA code)
#  *    2016 - 2020 Do Kester

class NoiseScale( HyperParameter ):
    """
    Hyperparameter for the scale of a ScaledErrorDistribution

    it is a measure of the noise.

    Information about the scale of the noise is stored in his class.
    It is either in the form of a fixed number, when the noise scale
    is known or in the form of a Prior with limits.
    By default this prior is a JeffreysPrior..

    The full use of priors is reserved for Bayesian calculations as
    in NestedSampler

    Attributes
    ----------
    scale : float
        the value of the noiseScale.  Default: 1.0
    stdev : float
        the standard deviation of the noise scale.  Default: None
    prior : Prior
        the prior for the noiseScale.  Default: JeffreysPrior
    fixed : boolean
        keep the noise scale fixed at the value given by scale.
        default: True
    minimum : boolean
        automatic noise scaling with a minimum. default: False

    """

    #  *********CONSTRUCTORS***************************************************
    def __init__( self, scale=1.0, isFixed=True, prior=None, limits=None,
                    copy=None ):
        """
        Constructor.

        Parameters
        ----------
        scale : float
            float   value of the noise scale
        isFixed : bool
            True:   Consider the hyperparameter as fixed
            False:  Optimize the parameter too (when relevant)
                    It might need a prior and/or limits to be set
                    The default prior is JeffreysPrior
        prior : None or Prior
            None : no prior set
            Prior : the prior probability on scale
        limits : None or list of 2 floats
            None : no limits set
            [lo,hi] : limits to be passed to the Prior.
            If limits are set, the default for Prior is JeffreysPrior
        copy : NoiseScale
            NoiseScale to copy

        """
        if limits is not None and prior is None :
            prior = JeffreysPrior()
        super( NoiseScale, self ).__init__( hypar=scale, isFixed=isFixed,
                        prior=prior, limits=limits )
        self.minimum = False

        if copy is not None :
            self.minimum = copy.minimum

    def copy( self ):
        """ Return a copy.  """
        return NoiseScale( scale=self.scale, copy=self )

    def __setattr__( self, name, value ) :
        """
        Rename scale to hypar and stdevScale to stdev.
        """
        if name == "scale" :
            self.hypar = value
        elif name == "stdevScale" :
            self.stdev = value
        else :
            object.__setattr__( self, name, value )


    def __getattr__( self, name ) :
        if name == "scale" :
            return self.hypar
        elif name == "stdevScale" :
            return self.stdev
        else :
            raise AttributeError( "Unknown attribute " + name )

        return None

    def minimumScale( self, scale=None ) :
        """
        Fit the noise scale with a minimum value.

        Parameters
        ----------
        scale : float
            the value of the noise scale. Default: noiseScale.scale

        """
        if scale is not None : self.hypar = scale
        self.minimum = True

    def __str__( self ) :
        return str( "Noise scale. value = %f" % self.hypar )


