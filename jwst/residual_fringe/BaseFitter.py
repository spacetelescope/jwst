import numpy as numpy
import math
import warnings
import matplotlib.pyplot as plt
from astropy.table import Table

from .ImageAssistant import ImageAssistant
from .MonteCarlo import MonteCarlo
from .ConvergenceError import ConvergenceError
from . import Tools
from . import Plotter
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

class BaseFitter( object ):
    """
    Base class for all Fitters.

    The Fitter class is to be used in conjunction with *Model classes.

    The Fitter class and its descendants fit data to a model. Fitter itself
    is the variant for linear models, ie. models linear in its parameters.

    For both linear and nonlinear models it holds that once the optimal
    estimate of the parameters is found, a variety of calculations is exactly
    the same: standard deviations, noise scale, evidence and model errors.
    They all derive more or less from the inverse Hessian matrix ( aka the
    covariance matrix ). All these calculations are in this Fitter class.
    Other Fitter classes relegate their calculation in these issues to this one.

    Examples
    --------
    # It is not possible to use this class. User Fitter, CurveFitter etc. in stead

    Note Also
    ---------
    1. The calculation of the evidence is an Gaussian approximation which is
       only exact for linear models with a fixed scale.
    2. Attributes labelled as read only should not be set by a user.

    Author: Do Kester

    Attributes
    ----------
    model : Model
        the model to be fitted
    xdata : array_like
        independent variable(s)
    nxdata : int (read only)
        length of the xdata vector(s)
    ndim : int (read only)
        number of xdata vectors
    weights : array_like
        the weights on the data from the last fit
    imageAssistant : ImageAssistant
        to convert images to pixels indices, needed for a fit
    keep : dict of {int : float}
        to keep the indexed (int) parameter at the provided value (float)
    fitIndex : list of int (or None)
        list of parameter indices to fit (None is all)
    npfit : int
        the number of parameters fitted in the last fit.
    fixedScale : float
        the fixed noise scale.
        The presence of `fixedScale` has consequences for the definitions of `chisq`,
        `(co)variance`, `stdevs` and `evidence`

    minimumScale : float
        introduce a minimum value for the noise scale
    design : matrix (read only)
        the design matrix (partial of model to parameters)
        returns self.getDesign()

    Attributes (available after a call to fit())
    ----------
    yfit : array_like
        The model result at the optimal value for the parameters.
        If map is true, a map is returned.
    chisq : float (read only)
        chisquared of the fit
    parameters : ndarray
        parameters fitted to the model
        returns self.model.parameters
    stdevs, standardDeviations : array_like (read only)
        the standard deviations on the parameters
        returns self.getStandardDeviations()
    hessian : matrix (read only)
        the hessian matrix
    covariance : matrix (read only)
        the covariance matrix
        returns self.getCovarianceMatrix()
    scale : float
        the noise scale
        returns self.getScale()
    sumwgt : float (read only)
        sum of the weights
    logZ : float (read only)
        the e-log of the evidence
        returns self.getLogZ()
    evidence : float (read only)
        the 10log of the evidence (logZ / log(10))
        returns self.getEvidence()

    Attributes (available after a call to getLogZ() or getEvidence())
    ----------
    logOccam : float (read only)
        Occam factor
    logLikelihood : float (read only)
        log of the likelihood

    """

    #  *****CONSTRUCTORS********************************************************
    def __init__( self, xdata, model, map=False, keep=None, fixedScale=None ):
        """
        Create a new Fitter, providing inputs and model.

        A Fitter class is defined by its model and the input vectors ( the
        independent variable ). When a fit to another model and/or another
        input vector is needed a new object should be created.

        Parameters
        ----------
        xdata : array_like
            independent input variable(s)
        model : Model
            the model function to be fitted
        map : bool (False)
            When true, the xdata should be interpreted as a map.
            The fitting is done on the pixel indices of the map,
            using ImageAssistant
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)
            The values of keep will be used by the Fitter as long as the Fitter exists.
            See also fit( ..., keep=dict )
        fixedScale : None or float
            None : the noise scale is not fixed
            float: value of the fixed noise scale
            The value of fixedScale only influences the evidence calculation

        Raises
        ------
        ValueError when one of the following is true
            1. Dimensionality of model and input does not match.
            2. Nans in input stream.
            3. Model is not the head of a compound model chain.

        """
        if model != model._head:
            raise ValueError( "Model is not the head of a compound model chain" )

        if map :
            self.imageAssistant = ImageAssistant()
            self.xdata = self.imageAssistant.getIndices( xdata )
        else :
            self.imageAssistant = None
            self.xdata = numpy.asarray( xdata )

        if isinstance( xdata, Table ) :
            ndim = len( xdata.columns )
        else :
            if numpy.any( numpy.isnan( xdata ) ) :
                raise ValueError( "NaNs in xdata array" )

            ndim = 1 if self.xdata.ndim <= 1 else self.xdata.shape[1]

        ninp = Tools.length( xdata )
        self.nxdata = ninp
        self.ndim = ndim
        self.model = model
        self.keep = keep
        self.fitIndex = self.keepFixed( keep )
        self.fixedScale = fixedScale

        if self.ndim != model.ndim:
            raise ValueError( "Model (%d) and xdata (%d) must be of the same dimensionality."
                                % (model.ndim, self.ndim) )

    def setMinimumScale( self, scale=0 ) :
        """
        Introduce a minimum in scale calculation and consequently in chisq.
            chi^2 >= sumwgt * scale^2

        Parameters
        ----------
        scale : float
            minimum scale
        """
        self.minimumScale = scale


    def fitprolog( self, ydata, weights=None, keep=None ) :
        """
        Prolog for all Fitters.

        1. Checks data/weighs for Nans
        2. Makes fitIndex.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted
        weights : array_like
            weights pertaining to the data
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)

        Returns
        -------
        fitIndex : ndarray of int
            Indices of the parameters that need fitting

        """
        self.checkNan( ydata, weights=weights )

        if self.imageAssistant is not None :
            ydata = self.imageAssistant.getydata( ydata )
            if weights is not None :
                weights = self.imageAssistant.getydata( weights )
        else :
            ydata = numpy.asarray( ydata )
            if weights is not None :
                weights = numpy.asarray( weights )

        self.weights = weights

        if keep is not None :
            return ( self.keepFixed( keep ), ydata, weights )

        if self.fitIndex is None :
            self.npfit = self.model.npchain
            return ( numpy.arange( self.model.npchain, dtype=int ), ydata, weights )

        self.npfit = len( self.fitIndex )
        return ( self.fitIndex, ydata, weights )

    def fitpostscript( self, ydata, plot=False ) :
        """
        Produce a plot of the results.
        """
        if plot :
            self.plotResult( xdata=self.xdata, ydata=ydata, model=self.model,
                    residuals=True, confidence=False, show=True )

# *****KEEP**************************************************************
    def keepFixed( self, keep=None ) :
        """
        Keeps parameters fixed at the provided values.

        1. The model will act exactly as if it were a model with less
           parameters, although slightly less efficient.
        2. Repeated calls start from scratch.
        3. Reset with keepFixed( None )

        Parameters
        ----------
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)

        Returns
        -------
        fitIndex : list of int (or None)
            list of parameter indices to be kept
        """
        if keep is None :
            self.npfit = self.model.npchain
            return None
        else :
#            print( "Keep  ", keep )
            self.npfit = len( keep )
            fitIndex = numpy.arange( self.model.npchain )      # start from scratch
            fitIndex = numpy.setxor1d( fitIndex, list( keep.keys() ) )
            self.model.parameters[list(keep.keys())] = list( keep.values() )
            return fitIndex

    def insertParameters( self, fitpar, index=None, into=None ) :
        """
        Insert fitparameters into the parameters when fitIndex is present.

        Parameters
        ----------
        fitpar : list of float
            (fitted) parameters
        index : list of int
            list of parameter indices to be kept
        into : list of float
            array into which the fitpar need to be inserted.

        """
        if index is None :
            index = self.fitIndex

        fitpar = numpy.array( fitpar, dtype=float, copy=False, ndmin=1 )

        if index is None or len( fitpar ) == len( self.model.parameters ) :
            return fitpar
        pars = self.model.parameters.copy() if into is None else into
        pars[index] = fitpar
        return pars

    #  *****FIT*****************************************************************
    def modelFit( self, ydata, weights=None, keep=None ):
        """
        Return model fitted to the data.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        weights : None or array_like
            weights to be used
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)
            The values will override those at initialization.
            They are only used in this call of fit.

        """
        self.model.parameters = self.fit( ydata, weights=weights, keep=keep )
        return self.yfit


    def limitsFit( self, fitmethod, ydata, weights=None, keep=None ) :
        """
        Fit the data to the model.
        When a parameter(s) transgresses the limits, it set and fixed at that limit
        and the fit is done again, excluding the parameter(s)
        When the chisq landscape is largely monomodal (no local minima) this is OK.

        Parameter
        ---------
        fitmethod : callable fitmethod( ydata, weights=weights )
            A fit method from the BaseFitter family
        ydata : array_like
            data that the model needs to be fit to
        weights : array_like
            weights partaining to the data.
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)
            The values will override those at initialization.
            They are only used in this call of fit.

        Returns
        -------
        pars : array_like
            the parameters of the fit

        Raises
        ------
        Warning when parameters have been reset at the limits.

        """
        pars = fitmethod( ydata, weights=weights )          # perform the fit
        if self.model.priors is None :                           # no priors -> no limits
            return pars

        sfix = self.keep                                    # save original fixed setting

        npchain = self.model.npchain
        params = self.model.parameters
        index = self.fitIndex
        if index is None :
            index = range( npchain )
            fix = []
        else :
            fix = sfix.copy()

        ool = 0
        for k in index :                                    # check limits for params
            if params[k] < self.model.getPrior( k ).lowLimit :
                fix = fix + [k]                             # add to original fixed
                params[k] = self.model.getPrior( k ).lowLimit
                ool += 1
            elif params[k] > self.model.getPrior( k ).highLimit :
                fix = fix + [k]
                params[k] = self.model.getPrior( k ).highLimit
                ool += 1

        if ool > 0 :                                        # some transgressions
            self.keepFixed( fix )                           # fix them
            pars = fitmethod( ydata, weights=weights )      # run fit again
            self.keepFixed( sfix )                          # reset original fitIndex
            if sfix is not None :                           # select parameters
                pars = [self.model.parameters[k] for k in self.fitIndex]
                fix = list( set( fix ) - set( sfix ) )
            warnings.warn( "Parameters ", fix, " exceeded limits." )

        return pars                                         # return parameters


    def fit( self, ydata, weights=None, keep=None ) :
        """
        Return model parameters fitted to the data.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        weights : array_like
            weights to be used
        keep : dict of {int:float}
            dictionary of indices (int) to be kept at a fixed value (float)
            The values will override those at initialization.
            They are only used in this call of fit.

        Raises
        ------
        NotImplementedError. BaseFitter cannot perform fits by itself.

        """
        raise NotImplementedError( "BaseFitter is a base class, not suitable itself to perform fits." )

    def checkNan( self, ydata, weights=None ):
        """
        Check there are no Nans or Infs in ydata or weights.

        Parameters
        ----------
        ydata : array_like
            data to be fitted.
        weights : array_like
            weights pertaining to ydata

        Raises
        ------
        ValueError.
        """
        if not numpy.all( numpy.isfinite( ydata ) ) :
            raise ValueError( "Fitter: NaNs or Infs in ydata" )
        if  weights is not None and not numpy.all( numpy.isfinite( weights ) ) :
            raise ValueError( "Fitter: NaNs or Infs in weights" )

    def __getattr__( self, name ) :
        """
        Return value belonging to attribute with name.

        Parameters
        ----------
        name : string
            name of the attribute
        """
        if name in ['logOccam', 'logLikelihood', 'chisq'] :
            raise AttributeError( str( self ) + ": " + name + " is not yet available." )
        elif name == 'parameters' :
            return self.model.parameters
        elif name == 'weights' :            ## not present return None
            return None
        elif name == 'sumwgt' :             ## not present return nxdata
            return self.nxdata
        elif name == 'yfit' :
            yfit = self.model.result( self.xdata )
            return yfit if self.imageAssistant is None else self.imageAssistant.resizeData( yfit )
        elif name == 'design' :
            return self.getDesign()
        elif name == 'hessian' :
            return self.getHessian()
        elif name == 'covariance' :
            return self.getCovarianceMatrix()
        elif name == 'stdevs' or name == 'standardDeviations' :
            return self.getStandardDeviations()
        elif name == 'scale' :
            return self.getScale()
        elif name == "stdevScale" :
            return self.scale / math.sqrt( 2 * self.sumwgt )
#            return self.scale / math.sqrt( 2 * self.nxdata )
        elif name == 'evidence' :
            return self.getEvidence()
        elif name == 'logZ' :
            return self.getLogZ()
        else :
            raise AttributeError( str( self ) + ": Unknown attribute " + name )

        return None


    #  *****VECTOR**************************************************************
    def getVector( self, ydata, index=None ):
        """
        Return the &beta;-vector.

        It includes "normalized" data if present. See normalize().

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted. When such is appliccable, it should be
            multiplied by weights and/or appended by normdata.
        index : list of int
            index of parameters to be fixed

        """
        if self.model.isNullModel() :
            return numpy.asarray( 0 )
        design = self.getDesign( index=index )

        return numpy.inner( design.transpose(), ydata )

    #  *****HESSIAN**************************************************************
    def getHessian( self, params=None, weights=None, index=None ):
        """
        Calculates the hessian matrix for a given set of model parameters.

        It includes "normalized" data if present. See normalize()

        Parameters
        ----------
        params : array_like
            the model parameters to be considered
        weights : array_like
            weights to be used
        index : list of int
            index of parameters to be fixed

        """
        if params is None : params = self.model.parameters
        if weights is None : weights = self.weights

        self.sumwgt = self.nxdata if weights is None else numpy.sum( weights )

        if self.model.isNullModel() :
            return

        design = self.getDesign( xdata=self.xdata, params=params, index=index )

        if hasattr( self, "normweight" ) :
            if weights is None :
                weights = numpy.ones( self.nxdata, dtype=float )
            weights = numpy.append( weights, self.normweight )

        design = design.transpose()

        if weights is not None :
            self.hessian = numpy.inner( design, design * weights )
        else :
            self.hessian = numpy.inner( design, design )

        return self.hessian

#      * TBD Condition number see Wikipedia: Condition Number and Matrix Norm

    #  *************************************************************************
    def getInverseHessian( self, params=None, weights=None, index=None ):
        """
        Return the inverse of the Hessian Matrix, H.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        weights : array_like
            weights to be used
        index : list of int
            index of parameters to be fixed

        """
        hes = self.getHessian( params, weights, index )
        inh = numpy.linalg.inv( hes )
        return inh

    def getCovarianceMatrix( self ):
        """
        Returns the inverse hessian matrix over the fitted parameters,
                multiplied by the variance.

        Stdevs are found from this as np.sqrt( np.diag( covarianceMatrix ) )

        """
        cov =  self.getInverseHessian( index=self.fitIndex ) * self.makeVariance()
        return cov

    #  *************************************************************************
    def makeVariance( self, scale=None ):
        """
        Return the (calculated) variance of the remaining noise. I.e.
            var = chisq / dof
        when automatic noise scaling is requested or
            var = scale^2
        when we have a fixed scale.

        Parameters
        ----------
        scale : float
            noise scale to be used
        """
        if scale is not None :
            return scale * scale

        scale = self.getScale()
        var = scale * scale
        if hasattr( self, "minimumScale" ) :       # add minimum to scale when requested
            var += self.minimumScale * self.minimumScale
        return var

    def normalize( self, normdfdp, normdata, weight=1.0 ) :
        """
        If for some reason the model is degenerate, e.g when two parameters measure
        essentially the same thing, This method can disambiguate these parameters.

        It is like adding a dummy measurement of one (or more) parameter to the data.

        Parameters
        ----------
        normdfdp : array_like
            for each parameter to sum to value (same length as self.parameters)
        normdata : float
            simulated data value
        weight : float
            weight of this measurement
        """
        normdfdp = numpy.asarray( normdfdp, dtype=float )
        if normdfdp.ndim == 1 :
            normdfdp = normdfdp.reshape( 1, len( normdfdp ) )
            normdata = numpy.asarray( [normdata], dtype=float )
            weight   = numpy.asarray( [weight], dtype=float )

        if not hasattr( self, "normdfdp" ) :
            self.normdfdp = normdfdp
            self.normdata = normdata
            self.normweight = weight
        else :
            self.normdfdp = numpy.append( self.normdfdp, normdfdp, axis=0 )
            self.normdata = numpy.append( self.normdata, normdata )
            self.normweight = numpy.append( self.normweight, weight )


    #  *****DESIGN**************************************************************
    def getDesign( self, params=None, xdata=None, index=None ):
        """
        Return the design matrix, D.
        The design matrix is also known as the Jacobian Matrix.

        Parameters
        ----------
        xdata : array_like
            the independent input data
        params : array_like
            parameters of the model
        index : list of int
            index of parameters to be fixed

        """
        if params is None : params = self.model.parameters
        if xdata is None :  xdata = self.xdata

        design = self.model.partial( xdata, params )
        if hasattr( self, "normdfdp" ) :
            design = numpy.append( design, self.normdfdp, axis=0 )

        if index is not None :
            design = design[:,index]

        return design

    #  *****CHI-SQUARED*********************************************************
    def chiSquared( self, ydata, params=None, weights=None ):
        """
        Calculates Chi-Squared for data and weights.

        It is the (weighted) sum of the squared residuals.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        params : array_like
            parameters for the model
        weights : array_like
            weights to be used

        Raises
        ------
        ValueError when chisq <= 0.

        """
        res2 = numpy.square( ydata - self.model.result( self.xdata, params ) )
        if weights is not None:
            res2 *= weights
            self.sumwgt = numpy.sum( weights )
        else:
            self.sumwgt = self.nxdata
        self.chisq = numpy.sum( res2 )
        if self.chisq <= 0 :
            raise ValueError( str( self ) + ": chisq <= 0" )

        return self.chisq

    #  *****STANDARD DEVIATION**************************************************
    def getStandardDeviations( self ):
        """
        Calculates of standard deviations pertaining to the parameters.

            &sigma;_i = s * sqrt( C[i,i] )

        where C is the Covariance matrix, the inverse of the Hessian Matrix and
        s is the noiseScale.

        Standard deviation are calculated for the fitted parameters only.

        Note that the stdev will decrease with sqrt( N ) of the number of
        datapoints while the noise scale, s, does not.

        """
        cov = self.covariance

        stdevs = numpy.sqrt( cov.diagonal() )
        if self.keep is None :
            self.model.stdevs = stdevs
        else :
            intostd = numpy.zeros( self.model.npchain, dtype=float )
            self.model.stdevs = self.insertParameters( stdevs, into=intostd )

        return self.model.stdevs

    #  *****MONTE CARLO ERROR***************************************************
    def monteCarloError( self, xdata=None, monteCarlo=None):
        """
        Calculates &sigma;-confidence regions on the model given some inputs.

        From the full covariance matrix (inverse of the Hessian) random
        samples are drawn, which are added to the parameters. With this new
        set of parameters the model is calculated. This procedure is done
        by default, 25 times.
        The standard deviation of the models is returned as the error bar.

        The calculation of the confidence region is delegated to the class
        MonteCarlo. For tweaking of that class can be done outside BaseFitter.

        Parameters
        ----------
        xdata : array_like
            input data over which to calculate the error bars.
        monteCarlo : MonteCarlo
            a ready-made MonteCarlo class.

        """
        if xdata is None : xdata = self.xdata
        if monteCarlo is None :
            monteCarlo = MonteCarlo( xdata, self.model, self.covariance,
                                          index=self.fitIndex )

        return monteCarlo.getError( xdata )

    #  *************************************************************************
    def getScale( self ):
        """
        Return
        ------
        float : the noise scale
            scale = sqrt( chisq / DOF )

        Raise
        -----
        RuntimeError when DoF <= 0. The number of (weighted) datapoints is too small.

        """
#        dof = self.nxdata - self.npfit
        dof = self.sumwgt - self.npfit
        if dof <= 0 :
            raise RuntimeError( "More parameters than (weighted) data points" )

        return math.sqrt( self.chisq / dof )

    #  *****EVIDENCE************************************************************
    def getEvidence( self, limits=None, noiseLimits=None ):
        """
        Calculation of the evidence, log10( Z ), for the model given the data.

            E = log10( P( Model | data ) )

        The calculation of the evidence uses a Gaussion approximation of the Posterior
        probability.
        It needs to know the limits of the parameters (and the noise scale if applicable),
        either from the priors in the model or from keywords "limits/noiseLimits".


        Parameters
        ----------
        limits : list of 2 floats/array_likes
            possible range of the parameters. ( [low,high] )
        noiseLimits : list of 2 floats
            possible range on noise scale ( [low,high] )

        Raises
        ------
        ValueError when no Prior is available

        """
        return self.getLogZ( limits, noiseLimits ) / math.log( 10.0 )

    def getLogLikelihood( self, autoscale=False, var=1.0 ) :
        """
        Return the log likelihood.

        It is implementing eq 19/20 last parts (Kester 2002) term by term

        Parameters
        ----------
        autoscale : bool
            whether the noise scale is optimized too
        var : float
            variance
        """
        return -0.5 * ( self.sumwgt * math.log( 2 * math.pi * var ) +
                        self.chisq / var )

    def getLogZ( self, limits=None, noiseLimits=None ):
        """
        Calculation of the evidence, log( Z ), for the model given the data.

            logZ = log( P( Model | data ) )

        The calculation of the evidence uses a Gaussion approximation of the Posterior
        probability.
        It needs to know the limits of the parameters (and the noise scale if applicable),
        either from the priors in the model or from keywords "limits/noiseLimits".


        Parameters
        ----------
        limits : list of 2 floats/array_likes
            possible range of the parameters. ( [low,high] )
        noiseLimits : list of 2 floats
            possible range on noise scale ( [low,high] )

        Raises
        ------
        ValueError when no Prior is available

        RuntimeError when DoF <= 0. The number of (weighted) datapoints is too small.

        """
#        dof = self.nxdata - self.npfit
        dof = self.sumwgt - self.npfit
        if dof <= 0 :
            raise RuntimeError( "More parameters than (weighted) data points" )

        if noiseLimits is not None :
            scalerange = math.log( noiseLimits[1] ) - math.log( noiseLimits[0] )
            autoScale = True
            s2 = self.makeVariance( )
            self.logOccam = 0.5 * math.log( math.pi * dof ) - scalerange
        else :
            autoScale = False
            scale = 1.0 if self.fixedScale is None else self.fixedScale
            s2 = scale * scale
            self.logOccam = 0.0

        # add integral over the scale prior, if scale is to be fitted.
        spr = math.log( scalerange ) if autoScale else 0.0

        # obtain the loglikelihood.
        self.logLikelihood = self.getLogLikelihood( autoscale=autoScale, var=s2 )

        if self.npfit == 0 :
            return self.logLikelihood + self.logOccam

        priors = self.model.priors
        if limits is not None :
            prirange = Tools.toArray( limits[1] ) - Tools.toArray( limits[0] )
        elif self.model.hasLimits() :
            prirange = [p.getIntegral() for p in priors]        # convert to list
        elif self.model.npchain > 0 :
            raise ValueError( "No limits provided on the parameters. Cannot calculate evidence." )

        priorlength = Tools.length( prirange )                         # maybe less than npar
        prirange = numpy.log( numpy.asarray( prirange ) )

        if priorlength == 1 :
            spr += prirange[0]
        elif limits is None and self.fitIndex is not None :
            fi = self.fitIndex
            q = numpy.where( fi < priorlength )
            fq = fi[q]
            priorlength = len( fq )
            spr += numpy.sum( prirange[fq] )
        else :
            spr += numpy.sum( prirange )

        if priorlength < self.npfit :                          # add enough times the last one
            spr += ( self.npfit - priorlength ) * prirange[-1]

        lidet = math.log( numpy.linalg.det( self.hessian ) )

        # implementing eq 18 (Kester 2002) term by term
        self.logOccam += -spr + 0.5 * ( self.npfit *
                         math.log( 2 * math.pi * s2 ) - lidet )

        return self.logLikelihood + self.logOccam

    def __str__( self ):
        """ Return name of the fitter.  """
        return "BaseFitter"


    def plotResult( self, xdata=None, ydata=None, model=None, residuals=True,
                    confidence=False, show=True ) :
        """
        Plot the results of the fit.

        Parameters
        ----------
        xdata : array_like
            xdata of the problem
        ydata : array_like
            ydata of the problem
        model : Model
            the model the ydata are fitted to at xdata.
        residuals : bool
            plot the residuals in a separate panel.
        confidence : bool
            plot confidence region
        show : bool
            display the plot.
        """
        if xdata is None :
            xdata = self.xdata
        if model is None :
            model = self.model
        fitter = self if confidence else None

        Plotter.plotFit( xdata, data=ydata, model=model, show=show,
                fitter=fitter, residuals=residuals )


