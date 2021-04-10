import numpy as numpy
import math as math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.testing import assert_array_almost_equal as assertAAE
from . import Tools
from .Formatter import formatter as fmt

__author__ = "Do Kester"
__year__ = 2020
__license__ = "GPL3"
__version__ = "2.5.3"
__url__ = "https://www.bayesicfitting.nl"
__status__ = "Perpetual Beta"

# This file is part of the BayesicFitting package.
#
# BayesicFitting is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or ( at your option ) any later version.
#
# BayesicFitting is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# The GPL3 license can be found at <http://www.gnu.org/licenses/>.
#
#   2017 - 2020 Do Kester

def plotFit( x, data=None, yfit=None, model=None, fitter=None, show=True,
             residuals=False ) :
    """
    Plot the data of a fit.

    Parameters
    ----------
    x : array_like
        xdata of the problem
    data : array_like
        ydata of the problem
    yfit : array_like
        fit of the data to the model
    model : Model
        the model the data are fitted to at x
    fitter : BaseFitter
        the fitter being used
        If set it displays a confidence region for the fit
    show : bool
        display the plot
    residuals : bool
        plot the residuals in a separate panel
    """

    minx = numpy.min( x )
    maxx = numpy.max( x )
    plt.figure( "Fitter Results" )
    ax0 = plt

    if yfit is None and model is not None :
        yfit = model( x )

    if residuals :
        plt.subplots_adjust( hspace=0.001 )
        gs = gridspec.GridSpec( 2, 1, height_ratios=[4, 1])

        ax1 = plt.subplot( gs[1] )
        res = data - yfit
        nd = int( math.log10( len( data ) ) )
        mrksz = ( 5 - nd ) if nd < 4 else 1
        ax1.plot( x, res, 'k.', markersize=mrksz )
        ax1.margins( 0.05, 0.05 )
        xt = plt.ylim()
        xtk = Tools.nicenumber( ( xt[1] - xt[0] ) / 4 )

        plt.yticks( [-xtk, 0.0, xtk] )
        plt.ylabel( "residual" )
        plt.xlabel( "xdata" )

        ax0 = plt.subplot( gs[0] )
        xticklabels = ax0.get_xticklabels()
        plt.setp( xticklabels, visible=False )

    if data is not None :
        nd = int( math.log10( len( data ) ) )
        mrksz = ( 5 - nd ) if nd < 4 else 1
        ax0.plot( x, data, 'k.', markersize=mrksz )

    if model is not None :
        xx = numpy.linspace( minx, maxx, 10 * len( x ) )
        yy = model( xx )
        if fitter is not None :
            err = fitter.monteCarloError( xx )
            ax0.plot( xx, yy - err, 'g-' )
            ax0.plot( xx, yy + err, 'g-' )
        ax0.plot( xx, yy, 'r-' )
    elif yfit is not None :
        ax0.plot( x, yfit, 'r-' )

    ax0.margins( 0.05, 0.05 )
    plt.ylabel( "ydata" )
    if not residuals :
        plt.xlabel( "xdata" )

    if show :
        plt.show()

