import numpy as numpy
import math

from .ErrorDistribution import ErrorDistribution
from .NoiseScale import NoiseScale
from .JeffreysPrior import JeffreysPrior
from . import Tools

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


class ScaledErrorDistribution( ErrorDistribution ):
    """
    Base class that defines methods common to error distributions with a scale.

    GaussErrorDistribution
    LaplaceErrorDistribution
    CauchyErrorDistribution
    ExponentialErrorDistribution
    UniformErrorDistribution

    Author       Do Kester.

    """
    PARNAMES = ["scale"]

    #  *********CONSTRUCTORS***************************************************
    def __init__( self, scale=1.0, limits=None, fixed=None, copy=None ):
        """
        Default Constructor.

        Parameters
        ----------
        scale : float
            noise scale
        limits : None or list of 2 floats [low,high]
            None : no limits implying fixed scale
            low     low limit on scale (needs to be >0)
            high    high limit on scale
            when limits are set, the scale is to be fitted
        fixed : dictionary of {int:float}
            int     list if parameters to fix permanently. Default None.
            float   list of values for the fixed parameters.

        copy : ScaledErrorDistribution
            distribution to be copied.

        """
        super( ScaledErrorDistribution, self ).__init__( fixed=fixed, copy=copy )

        if copy is None :
            self.hyperpar = NoiseScale( scale=scale, limits=limits )
        else :
            self.hyperpar = copy.hyperpar                   ## TBC copy ???

        if limits is None :
            self.fixed = self.keepFixed( {0: scale} )

    def copy( self ):
        """ Return copy of this.  """
        return ScaledErrorDistribution( copy=self )

    def setLimits( self, limits ) :
        """
        Set limits for scale.

        Parameters
        ----------
        limits : [low,high]
            low : float or array_like
                low limits
            high : float or array_like
                high limits
        """
        if self.hyperpar[0].prior is None :
            self.hyperpar[0].prior = JeffreysPrior()

        super( ScaledErrorDistribution, self ).setLimits( limits )


    def __setattr__( self, name, value ):
        """
        Set attributes.

        """
        if name == "scale" and Tools.isInstance( value, float ) :
            self.hyperpar[0].hypar = value
        else :
            super( ScaledErrorDistribution, self ).__setattr__( name, value )

    def __getattr__( self, name ) :
        """
        Return value belonging to attribute with name.

        Parameters
        ----------
        name : string
            name of the attribute
        """
        if name == "scale" :
            return self.hyperpar[0].hypar
        else :
            return super( ScaledErrorDistribution, self ).__getattr__( name )


