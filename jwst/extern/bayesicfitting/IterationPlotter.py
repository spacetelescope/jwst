import time
import numpy as numpy
from astropy import units
import math
from . import Tools
import matplotlib.pyplot as pyplot

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
#  *    2016 - 2020 Do Kester

class IterationPlotter( object ):
    """
    The IterationPlotter plots intermediate results from a iterative fitter.

    Author:      Do Kester

    """
    def plotData( self, x, y, title ):
        """
        Plot the data.

        Parameters
        ----------
        x : array_like
            x-axis values of the data
        y : array_like
            y-axis values of the data
        title : string
            the title of the plot
        """

        self.p = pyplot.gca()
        self.p.plot( x, y, 'k.' )
        pyplot.show( block=False )
        time.sleep( 1 )

    def plotResult( self, x, r, iter ):
        """
        Plot the ( intermediate ) result.

        Parameters
        ----------
        x : array_like
            x-axis values of the data
        r : array_like
            model result
        iter : int
            iteration number
        """
        self.p.plot( x, r, 'r-' )
        pyplot.show( block=False )
        time.sleep( 1 )

    def plotProgress( self, percent ):
        """
        Plot (estimated) progress upto now.

        """
        pass


