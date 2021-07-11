import numpy as numpy
from . import Tools
import math

from .Formatter import formatter as fmt

__author__ = "Do Kester"
__year__ = 2020
__license__ = "GPL3"
__version__ = "2.6.1"
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
#  *    2003 - 2014 Do Kester, SRON (JAVA code)
#  *    2016 - 2020 Do Kester

class MonteCarlo( object ):
    """
    Helper class to calculate the confidence region of a fitted model.

    MonteCarlo for models.

    The MonteCarlo class is to be used in conjunction with Model classes.

    Author:      Do Kester

    Attributes
    ----------
    xdata : array_like
        array of independent input values
    model : Model
        the model function to be fitted
    mcycles : int
        Sets number of cycles in the MonteCarlo procedure to estimate
        error bars. Default = 25

    Hidden Attributes
    -----------------
    _eigenvectors : array_like (read only)
        from eigenvalue decomposition of covariance matrix
    _eigenvalues : array_like (read only)
        from eigenvalue decomposition of covariance matrix
    _random : random
        random number generator

    """

    #  *****CONSTRUCTORS********************************************************
    def __init__( self, xdata, model, covariance, index=None, seed=12345, mcycles=25 ):
        """
        Create a new MonteCarlo, providing inputs and model.

        A MonteCarlo object is defined by its model and the input vector (the
        independent variable). When a fit to another model and/or another
        input vector is needed a new object should be created.

        Parameters
        ----------
        xdata : array_like
            array of independent input values
        model : Model
            the model function to be fitted
        covariance : matrix
            the covariance matrix of the problem. Default from the Model.
        index : list of int
            indices of parameters to fit
        seed : int
            seed for random number generator
        mcycles : int
            number of cycles in the MonteCarlo procedure to estimate error bars.

        Raises
        ------
        ValueError when model and input have different dimensions

        """

        self.xdata = xdata
        self.nxdata = Tools.length( xdata )
        self.model = model
        self.index = index
        self.mcycles = mcycles
        self._random = numpy.random
        self._random.seed( seed )
        self.decompose( covariance )

    def decompose( self, covariance ):
        val, vec = numpy.linalg.eigh( covariance )
        self._eigenvectors = vec
        self._eigenvalues = numpy.sqrt( val )

    #  *****MONTE CARLO ERROR***************************************************
    def getError( self, xdata=None ):
        """
        Calculates 1 &sigma;-confidence regions on the model given some inputs.

        From the full covariance matrix ( = inverse of the Hessian ) random
        samples are drawn, which are added to the parameters. With this new
        set of parameters the model is calculated. This procedure is done
        by default, 25 times.
        The standard deviation of the models is returned as the error bar.

        Parameters
        ----------
        xdata ; array_like
            input data over which to calculate the error bars. default provided xdata

        """
        if xdata is None :
            xdata = self.xdata
        ninp = Tools.length( xdata )
        sm1 = numpy.zeros( ninp )
        sm2 = numpy.zeros( ninp )

        for k in range( self.mcycles ) :
            model = self.randomVariant( xdata )
            sm1 += model
            sm2 += numpy.square( model )
        sm1 /= self.mcycles
        sm2 /= self.mcycles
        return numpy.sqrt( sm2  - sm1 * sm1 )

    def randomVariant( self, xdata ):
        """
        Return a random variant of the model result.
        Taking into account the stdev of the parameters and their covariance.

        Parameters
        ----------
        xdata : array_like
            input data at these indpendent points

        """
        nfit = self.model.npchain if self.index is None else len( self.index )
        err = self._random.standard_normal( nfit )
        err = numpy.inner( self._eigenvectors, self._eigenvalues * err )
        par = self.model.parameters.copy()
        if self.index is None :
            par += err
        else :
            par[self.index] += err

        return self.model.result( xdata, par )

    #  *************************************************************************
    def __str__( self ):
        """ Return name of the class.  """
        return "MonteCarlo with seed %d and %d cycles" % (self._random.seed, self.mcycles)




