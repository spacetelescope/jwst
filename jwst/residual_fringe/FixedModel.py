import numpy as numpy
from astropy import units
import re
import string
import warnings
from . import Tools

from .BaseModel import BaseModel
from .Tools import setAttribute as setatt

__author__ = "Do Kester"
__year__ = 2021
__license__ = "GPL3"
__version__ = "2.7.0"
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
#  *   2017 - 2021 Do Kester

class FixedModel( BaseModel ):
    """
    A FixedModel is a BaseModel where some parameters are permanently fixed.

    A parameter can be fixed either with a constant float value or more
    dynamically, with another Model. The parameters of this latter model
    also appear as parameters of the FixedModel.

    The methods result (f(x:p)) and partial (df/dp) are calculated
    in here.
    Unfortunately the methods derivative (df/dx) is model dependent.
    It is reset to numDerivative.

    Examples
    --------
    >>> m1 = PolynomialModel( 1 )
    >>> m1 += SineModel()
    >>> print( m1.npchain )         # number of params: 2 + 3
    5
    >>> fixed = { 0: 1.0, 1: m1 }
    >>> em = EtalonModel( fixed=fixed )
    >>> print( em.npbase, em.npmax, em.npchain )          # ( 4 - 2 ) + 5
    7 9 7
    >>> print( em )

    Attributes
    ----------
    npmax : int
        maximum number of parameters of the simple (not-fixed) model
    fixed : dictionary of {int:float|Model}
        int     index of parameter to fix permanently. Default None.
        float   value for the fixed parameter.
        Model   model to replace the parameter
        Attribute fixed can only be set in the constructor.
    parlist : array_like or None
        list of active (not-fixed) indices. None is all.
    mlist : list of Model
        list of parameter indices which are replaced by a Model in fixed.

    Attributes from BaseModel
    --------------------------
        npbase, ndim, priors, posIndex, nonZero, tiny, deltaP, parNames

    Author       Do Kester

    """

    #  *************************************************************************
    def __init__( self, nparams=0, ndim=1, copy=None, fixed=None,
                  names=None, **kwargs ) :
        """
        FixedModel Constructor.

        Parameters
        ----------
        nparams : int
            Number of parameters in the model (default: 0)
        ndim : int
            Number of dimensions of the input (default: 1)
        copy : BaseModel
            to be copied
        fixed : dictionary of {int:float|Model}
            int     index of parameter to fix permanently. Default None.
            float   value for the fixed parameter.
            Model   model to replace the parameter
            Attribute fixed can only be set in the constructor.
        names : list of string
            names for the parameters
        kwargs for @BaseModel :
            posIndex, nonZero

        """

        mlist = []
        npar = nparams
        if fixed is not None :
            for k in fixed.keys() :
                if k < 0 :                  # convert negative indices into positive ones
                    val = fixed[k]
                    del( fixed[k] )
                    k += npar
                    fixed[k] = val
            for k in fixed.keys() :
                if isinstance( fixed[k], BaseModel ) :
                    mlist += [k]
                    npar += fixed[k].npchain

        super( FixedModel, self ).__init__( nparams=npar, ndim=ndim,
                copy=copy, **kwargs )


        if copy is None :
            setatt( self, "npmax", nparams )
            if fixed is not None :
                setatt( self, "npbase", self.npbase - len( fixed ) )
                setatt( self, "_npb", npar )
            setatt( self, "mlist", mlist )
            self._setFixed( fixed )
            if names is None :
                parNames = ["parameter_%d"%k for k in range( self.npbase )]
            else :
                parNames = self.selectNames( names )
        else :
            setatt( self, "npmax", copy.npmax )
            setatt( self, "mlist", copy.mlist )
            if copy.fixed is not None :
                fixed = {}
                setatt( self, "_npb", copy._npb )
                for k in copy.fixed.keys() :
                    if k in self.mlist :
                        fixed[k] = copy.fixed[k].copy()
                    else :
                        fixed[k] = copy.fixed[k]
            else :
                fixed = None
            parNames = copy.parNames[:]
        self._setFixed( fixed )
        setatt( self, "parNames", parNames )

    def copy( self ) :
        """ Return a copy. """
        return FixedModel( copy=self )

    def __setattr__( self, name, value ):
        """
        Set attributes.

        """
        if name == "fixed" :
            raise AttributeError( "Attribute fixed can only be set in the constructor" )
        else :
            super( FixedModel, self ).__setattr__( name, value )

    def _setFixed( self, fixed ) :
        setatt( self, "fixed", fixed )
        if fixed is None :
            setatt( self, "parlist", None )
            return

        if not isinstance( fixed, dict ) :
            raise AttributeError( "Attribute fixed needs to be a dictionary" )

        parlist = numpy.arange( self._npb )
        parlist = numpy.setxor1d( parlist, list( fixed.keys() ) )
        setatt( self, "parlist", parlist )

    def select( self, params ) :
        """
        Select the relevant parameters and store them.

        Parameters
        ----------
        params : array of float
            parameters of the head model
        """

        if self.fixed is None :
            return params

        # only fixed floats; np models
        if len( self.mlist ) == 0 :
            return params[self.parlist]

        pl = numpy.intersect1d( self.parlist, numpy.arange( self.npmax ) )
        pars = params[pl]
        for m in self.mlist :
            pars = numpy.append( pars, self.fixed[m].parameters )
        return pars

    def selectNames( self, names ) :
        """
        Select the relevant parameter names and store them.

        Parameters
        ----------
        names : list of string
            parameter names of the head model
        """
        if self.fixed is None :
            return names
        if len( self.mlist ) == 0 :
            return [names[k] for k in self.parlist]

        pl = numpy.intersect1d( self.parlist, numpy.arange( self.npmax ) )
        nms = [names[k] for k in pl]

        for m in self.mlist :
            mdl = self.fixed[m]
            while mdl is not None:
                nms += mdl.parNames
                mdl = mdl._next
        return nms

    #  *****RESULT**************************************************************
    def result( self, xdata, param ):
        """
        Returns the result calculated at the xdatas.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.

        """
        expparam = self.expand( xdata, param )
        return super( FixedModel, self ).result( xdata, expparam )

    def expand( self, xdata, param ) :
        """
        Returns a complete list of parameters, where the fixed parameters
        have been replaced by either a constant value or by the results of
        the fixed function.
        """
        if self.fixed is None :
            return param

        pars = numpy.zeros( self._npb, dtype=float )

        pars[self.parlist] = param
        if len( self.mlist ) == 0 :
            fl =  list( self.fixed.keys() )
            pars[fl] = list( self.fixed.values() )
            return pars

        par = pars.tolist()
        ns = self.npmax - len( self.fixed )
        for k in list( self.fixed.keys() ) :
            if k in self.mlist :
                nf = ns + self.fixed[k].npchain
                par[k] = self.fixed[k].result( xdata, param[ns:nf] )
                ns = nf
            else :
                par[k] = self.fixed[k]

        ## par is now a heterogeneous list of floats and [floats]
        ## It can only be cast into an array of object

        return numpy.asarray( par, dtype=object )

    #  *****PARTIAL*************************************************************
    def partial( self, xdata, param ):
        """
        Returns the partial derivatives calculated at the inputs.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.

        """
        expparam = self.expand( xdata, param )
        self.checkParameter( expparam )
        if len( self.mlist ) == 0 :
            return super(FixedModel, self).partial( xdata, expparam, parlist=self.parlist )

        np = len( param )
        ns = len( expparam )
        nf = self.npmax - len( self.fixed )
        plist = numpy.append( self.parlist[:nf], self.mlist )
        partial = numpy.zeros( ( Tools.length( xdata ), np ), dtype=float )

        bp = self.basePartial( xdata, expparam[:self.npmax], parlist=plist )
        partial[:,:nf] = bp[:,:nf]
        ns = nf
        nb = self.npmax
        for k in self.mlist :
            ne = self.fixed[k].npchain
            pa = expparam[nb:nb+ne]
            pf = self.fixed[k].partial( xdata, pa )
            ppp = bp[:,nf] * pf.transpose()
            partial[:,ns:ns+ne] = ppp.transpose()
            ns += ne
            nf += 1
            nb += ne
        return partial

    def numPartial( self, xdata, params, parlist=None ) :
        """
        Returns numerical partial derivatives of the model to params.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters; default is self.parameters
        parlist : array_like
            list of indices active parameters (or None for all)

        """
        nxdata = Tools.length( xdata )
        partial = numpy.zeros( ( nxdata, self.npbase ), dtype=float  )
        dpg = Tools.makeNext( self.deltaP, 0 )

        par = params.copy()
        if parlist is None :
            parlist = numpy.arange( self.npbase )

        i = 0

        for k in parlist :
            dp = next( dpg )
            par[k] += 0.5 * dp
            r1 = FixedModel.result( self, xdata, par )
            par[k] -= dp
            r2 = FixedModel.result( self, xdata, par )
            partial[:,i] = ( r1 - r2 ) / dp
#            print( k, i, params[k], dp, par[k], r1, r2, partial[:,i] )
            par[k] = params[k]
            i += 1
        try :
            dpg.close()
        except :
            pass
        return partial

    def basePartial( self, xdata, param, parlist=None ) :
        """
        Replacement for models that dont define a partial.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        param : array_like
            values for the parameters.
        parlist : array_like
            list of indices active parameters (or None for all)

        """
        warnings.warn( self.shortName() + " has no partials defined; " +
                    "using numeric partials." )

        if parlist is None :
            return FixedModel.numPartial( self, xdata, param, parlist=parlist )

        par = param[parlist]
        pl = numpy.arange( len( par ), dtype=int )
        return FixedModel.numPartial( self, xdata, par, parlist=pl )

    def derivative( self, xdata, param ) :
        """
        Returns the derivative of the model to xdata.

        It is a numeric derivative as the analytic derivative is not present
        in the model.

        If `fixed` contains a Model, the derivative cannot be constructed
        from the constituent models. Use `numDerivative` instaed.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the derivative
        param : array_like
            values for the parameters. (default: self.parameters)

        """
        if len( self.mlist ) == 0 :
            params = self.expand( xdata, param )
            return self.baseDerivative( xdata, params )
        else :
            return FixedModel.numDerivative( self, xdata, param )

    def numDerivative( self, xdata, param ) :
        """
        Returns the numeric derivative of the model to input

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the derivative
        param : array_like
            values for the parameters. (default: self.parameters)

        Raises
        ------
        ValueError when the number of xdata dimensions > 1.

        """
        if self.ndim > 1 :
            raise ValueError( "Not implemented for dimensions > 1" )

        dx = self.deltaP[0]
        xd = xdata + dx
        r1 = FixedModel.result( self, xd, param )
        xd = xdata - dx
        r2 = FixedModel.result( self, xd, param )

        return ( r1 - r2 ) / ( 2 * dx )


    #  *****TOSTRING***********************************************************
    def _toString( self, npars=0 ) :
        basename = self.baseName()

#        print( "FM %d  " % npars, basename )

        if self.fixed is not None :
            for k in self.fixed.keys() :
                par = r"p_%d([^0-9]|$)"%k
                if k in self.mlist :
                    val = r"(%s) " % self.fixed[k].shortName()
                else :
                    val = r"(%.1f) "%self.fixed[k]
                basename = re.sub( par, val, basename )
#            print( "FM1   ", basename )
#            print( "ParL  ", self.parlist )
            i = 0
            for k in self.parlist :
                old = r"p_%d([^0-9]|$)"%k
                par = r"q_%d\1"%i
                basename = re.sub( old, par, basename )
                i += 1
#                print( "FM    ", basename )
            basename = re.sub( r"q_", r"p_", basename )

#        print( "FM2   ", basename )


        if npars == 0 : return basename

#        print( "FM11  ", basename )

        i = 0
        while True :
            mo = re.search( r'p_([0-9]*)', basename, flags=re.MULTILINE )
            if mo is None : break

            span = mo.span()
            k = int( basename[span[0]+2:span[1]] )
#            print( "FM    ", span, k, npars )
            basename = re.sub( r'p_%d'%k, r'q_%d'%(k+npars), basename, flags=re.MULTILINE )


#        print( "FM12  ", basename )

        basename = basename.replace( "q_", "p_" )

#        print( "FM13  ", basename )

        return basename



