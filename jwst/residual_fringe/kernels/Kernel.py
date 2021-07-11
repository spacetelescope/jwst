
__author__ = "Do Kester"
__year__ = 2017
__license__ = "GPL3"
__version__ = "0.9"
__maintainer__ = "Do"
__status__ = "Development"

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
#  * 2013 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2017 Do Kester


class Kernel( object ):
    """
    A kernel is an even real-valued integrable function.

    It is satisfying the following two requirements:

    1. The integral over [-Inf,+Inf] exists ( < Inf ).
    2. It is an even function: K( x ) = K( -x ).

    A kernel is called bound when it is 0 everywhere
    except when |x| < range < inf.

    All kernels are scaled such that its value at x=0 is 1.0.


    Several kernel functions, K( x ) are defined in this package:

    Name        Definition          Integral    FWHM    range   comment
    Biweight    ( 1-x^2 )^2            16/15    1.08      1.0   aka Tukey
    CosSquare   cos^2( 0.5*PI*x )        1.0    1.00      1.0
    Cosine      cos( 0.5*PI*x )         4/PI    1.33      1.0
    Gauss       exp( -0.5*x*x )   sqrt(2*PI)    1.22      inf
    Huber       min( 1, 1/|x| )          inf    4.00      inf   improper kernel; aka Median
    Lorentz     1 / ( 1 + x*x )           PI    2.00      inf
    Parabola    1 - x*x                  4/3    1.41      1.0
    Sinc        sin(x) / x               1.0    1.21      1.0
    Triangle    1 - |x|                  1.0    1.00      1.0
    Tricube     ( 1 - |x|^3 )^3        81/70    1.18      1.0
    Triweight   ( 1 - x^2 )^3          32/35    0.91      1.0
    Uniform     1.0                      2.0    2.00      1.0   aka Clip

    Tophat 0    1.0                      1.0    1.00      0.5   like Uniform
           1    1 - |x|                  1.0    1.00      1.0   like Triangle
           2    2nd order polynomials    1.0    1.26      1.5
           3    3rd order polynomials    1.0    1.44      2.0
           4    4th order polynomials    1.0    1.60      2.5
           5    5th order polynomials    1.0    1.73      3.0
           6    6th order polynomials    1.0    1.86      3.5

    For all bound Kernels the definition in the table is true for |x| < range;
    elsewhere it is 0.

    Huber is not a proper Kernel as the integral is inf. However it is important
    in robust fitting (RobustShell) to get a madian-like solution for the outliers.

    Attributes
    ----------
    integral : float
        the integral over the valid range
    fwhm : float
        the full width at half maximum
    range : float
        the region [-range..+range] where the kernel is non-zero.

    Author:      Do Kester

    Category:    mathematics/Fitting

    """
    def __init__( self, integral=1.0, fwhm=1.0, range=1.0  ) :
        """
        Constructor

        Parameters
        ----------
        integral : float
            over [-inf, +inf]
        fwhm : float
            full width at half maximum
        range : float
            the region [-range,+range] where the kernel is nonzero
        """

        self.integral = integral
        self.fwhm = fwhm
        self.range = range

    def result( self, x ):
        """
        Return the result for one input value.

        Parameters
        ----------
        x : array-like
            the input value

        """
        pass

    def resultsq( self, xsq ):
        """
        Return the result for one input value.

        Parameters
        ----------
        x : array-like
            the input value
        Parameters: x the square of the input value

        """
        pass

    def partial( self, x ):
        """
        Return the partial derivative wrt input value.

        Parameters
        ----------
        x : array-like
            the input value

        """
        pass

    def isBound( self ):
        """
        Return true when the kernel is bound.
        All non-zero values are between -1 and +1

        """
        pass

    def name( self ):
        """
        Return the name of the kernel.
        """
        return "Kernel of unknown function"


    def __str__( self ):
        """ Return a name.  """
        return self.name()



