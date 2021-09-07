import numpy as numpy
import inspect

from .IterativeFitter import IterativeFitter
from .kernels.Kernel import Kernel
from .kernels.Biweight import Biweight

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
#  *    2008 - 2014 Do Kester, SRON (Java code)
#  *    2017 - 2020 Do Kester


class RobustShell( IterativeFitter ):
    """
    RobustShell tries to make a fit more robust in the presence of outliers.

    It is a shell around another fitter. Technically it is in itself a
    a Fitter, but with limited possiblities.

    A RobustShell tries to make a fit more robust in the presence of outliers.
    It does it by manipulating the weights: outliers are
    downweighted. "Normal" points keep their weights, more or less.

    Apart from methods specific to the robustification, RobustShell has a fit method
    and little else from the Fitter family. Methods to get the :math:`\chi^2`,
    the covariance matrix, the evidence, the noise scale etc. should be taken from
    the embedded fitter.

    The theory behind robust fitting can be found in Wikipedia: Robust Statistics,
    and the references therein.

    To determine which points are "normal" and which are "outliers" we
    need a previous fit. The difference to this earlier fit determines
    the normalcy. Its amount is used to adjust the weights.

    The fact that we need a previous iteration makes robust estimation a
    iterative procedure.
    Typically the procedure should stabilize in half a dozen steps maximally.

    There are several schemes to adjust the weights. But before we go into that
    we need two values in each scheme.
    Firstly the amount of noise present in the data. By default the noise
    is taken from the previous fit via Fitter.scale.
    The second value needed is the size of the influence domain in terms
    of the noise scale.

    For all schemes the deviant is calculated as the difference
    between data and fit divided by noisescale and domainsize:

    d = ( data - fit ) / ( noisescale * domainsize )

    The domainsize is defined such that deviants upto 3 times the noisescale fall
    within the fwhm.

    A number of weighting schemes are provided.

    With bound support and smooth edges:
       Kernel       Name        domain      comment
       Biweight     Tukey         5.54      Default kernel
       CosSquare                  6.00
       Tricube                    5.08
       Triweight                  6.60

    With bound support and hard edges:
       Uniform      Clip          3.00      Ignore all outside 3 sigma
       Cosine                     4.50
       Triangle                   6.00
       Parabola                   4.50

    with unbound support:
       Huber        Median        1.50      Inside domain mean; outside domain median
       Gauss                      2.12
       Lorentz                    3.00

    Other schemes can be written by making another Kernel or writing a function
        wgts = func( d )
    where d is the deviant as above.


    Notes
    -----
    Robust fitting is even more dangerous than ordinary fitting.
    *Never trust what you get without thorough checking.*

    Example
    -------
    >>> model = PolynomialModel( 1 )                # some model
    >>> x = numpy.arange( 100, dtype=float ) / 100  # some x values
    >>> y = numpy.arange( 100, dtype=float ) / 4    # digitization noise
    >>> y[1,11,35,67] += 10                         # create outliers
    >>> ftr = Fitter( x, model )                    # a fitter, a model and a x
    >>> rob = RobustShell( ftr, domain=7 )          # robust fitter using someFtr
    >>> par = rob.fit( y )                          # solution
    >>> print( rob )                                #
    >>> print( rob.weights )                        # print final weights
    >>> print( ftr.chisq )                          # get from someFtr

    Author       Do Kester.

    """
    def __init__( self, fitter, kernel=Biweight, domain=None, onesided=None, **kwargs ):
        """
        Create a new class, providing the fitter to be used.

        Parameters
        ----------
        fitter : BaseFitter
             to be used
        kernel : Kernel or callable
            All Kernels have a method `result( d )` which is applied to the deviants.
            where d = ( data - model ) / ( domain * scale )
            If kernel is a callable method it is assumed to be a similar result mathod.
        domain : None or float
            Width of the kernel.
            None : automatic calculation of domain according to table in class doc.
            float : overrides autocalculation.
        onesided : None or "positive" or "p" or "negative" or "n"
            None : apply robust weights to positive and negative residuals
            "positive" : apply robust weights to positive residuals only
            "negative" : apply robust weights to negative residuals only

        """
        super( RobustShell, self ).__init__( fitter.xdata, fitter.model, **kwargs )
        self.fitter = fitter
        self.kernelfunc = self.setKernel( kernel )
        if domain is None :
            domain = 6.0 / self.kernel.fwhm
        self.domain = domain
        self.onesided = self.setOneSided( onesided )

    def setKernel( self, kernel ) :
        """
        Set the robust kernel to be used.

        Parameters
        ----------
        kernel : Kernel or callable
            All Kernels have a method `result( d )` which is applied to the deviants.
            where d = ( data - model ) / ( domain * scale )
            If kernel is a callable method it is assumed to be a similar result mathod.

        Raises
        ------
        ValueError when kernel is not recognized.

        """
        if inspect.isclass( kernel ) :
            kernel = kernel()
        if isinstance( kernel, Kernel ) :
            self.kernel = kernel
            kernelfunc = self.kernel.result
        elif callable( kernel ) :
            self.kernel = None
            kernelfunc = kernel
        else :
            raise ValueError( "Could not interpret the value for kernel: %s" % str( kernel ) )

        return kernelfunc

    def setOneSided( self, onesided ) :
        """
        set self.onesided to either 0 or +1 or -1.

        Parameters
        ----------
        onesided : None or "positive" or "negative"
            None : apply robust weights to positive and negative residuals
            "positive" : apply robust weights to positive residuals only
            "negative" : apply robust weights to negative residuals only

        Raises
        ------
        ValueError when onesided could not be interpreted.

        """
        if onesided is None :
            return 0
        elif isinstance( onesided, str ) :
            if onesided[0] is "p" or onesided[0] is "P" :
                return +1
            elif onesided[0] is "n" or onesided[0] is "N" :
                return -1
            else :
                raise ValueError( "Unknown string for onesided: %s" % onesided )
        else :
            raise ValueError( "Could not interpret the value for onesided: %s" % str( onesided ) )


    def fit( self, data, weights=None, kernel=None, domain=None, onesided=None, **kwargs ) :
        """
        Perform a robustification step.

        Parameters
        ----------
        data : array_like
            the data as they go into a fitter
        kwargs : dict
            keyword args to be passed to fitter.fit()
        """

        kernelfunc = self.kernelfunc if kernel is None else self.setKernel( kernel )
        onesided = self.onesided if onesided is None else self.setOneSided( onesided )
        domain = self.domain if domain is None else domain

        self.iter = 0

        param = self.fitter.fit( data, weights, **kwargs )
        self.npfit = self.fitter.npfit
        chi = self.fitter.chisq

        while self.iter < self.maxIter :

            residuals = data - self.model.result( self.xdata, param )
            scale = self.fitter.scale
            residuals /= ( scale * domain )
            robwgt = kernelfunc( residuals )
            robwgt = self.getOneSidedWeights( robwgt, residuals, onesided )

            self.weights = robwgt
            if weights is not None :
                self.weights *= weights
            param = self.fitter.fit( data, self.weights, **kwargs )

            trychi = self.fitter.chisq
            tol = self.tolerance if trychi < 1 else self.tolerance * trychi

            if abs( chi - trychi ) < tol :
                self.report( self.verbose, param, trychi, more=scale, force=True )
                break

            self.report( self.verbose, param, trychi, more=scale, force=(self.verbose>=3) )
            chi = trychi
            self.iter += 1

        self.hessian = self.fitter.getHessian()
        self.chisq = trychi

        return param

    def getOneSidedWeights( self, wgt, res, onesided ) :
        if onesided == 0 :
            return wgt
        if onesided == +1 :
            return numpy.where( res < 0, 1.0, wgt )
        if onesided == -1 :
            return numpy.where( res > 0, 1.0, wgt )
        return wgt

    def __getattr__( self, name ) :
        return super( RobustShell, self ).__getattr__( name )

    def __str__( self ):
        """
        Return the name and weight of the fitter.
        """
        name = "RobustShell( " + str( self.fitter ) + " ) using a "
        if self.onesided == 0 :
            return name + self.kernel.name()
        if self.onesided == +1 :
            return name + " +1 sided " + self.kernel.name( )
        if self.onesided == -1 :
            return name + " -1 sided " + self.kernel.name( )
        return name + self.kernel.name()

