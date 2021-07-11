from __future__ import print_function

import numpy as numpy
import random
from astropy import units
import warnings

from .FixedModel import FixedModel
from . import Tools
from .Formatter import formatter as fmt
from .Tools import setAttribute as setatt

__author__ = "Do Kester"
__year__ = 2020
__license__ = "GPL3"
__version__ = "2.6.0"
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
#  *    2003 - 2011 Do Kester, SRON (JAVA code)
#  *    2016 - 2020 Do Kester


class Model( FixedModel ):
    """
    Model implements the common parts of (compound) models.
    It is the last common anchestor of all Models.

    Models can be handled by the Fitter classes.

    A model consists of one or more instantiations of (base) models which
    are concatenated in a chain of models using various operations
    (+-*/). A special operation is the pipe (|). It works like a unix pipe,
    i.e. the output of the left-hand process in used as input of the
    right-hand process.

    Methods defined in BaseModel as eg. baseResult() are recursively called
    here as result(). They are the ones used in the fitters.

    The Model is the place where model-related items are kept, like parameters,
    stdevs.

    Model also implements a numerical derivation of partial to be
    used when partial is not given in the model definition itself. This same
    numerical derivation of partial is used in testPartial to indeed test
    whether the partial has been implemented properly.

    Example:
    --------
    >>> x = numpy.arange( 10 )
    >>> poly = PolynomialModel( 2 )             # quadratic model
    >>> poly.parameters = [3,2,1]               # set the parameters for the model
    >>> y = poly( x )                           # evaluate the model at x
    >>> p0 = poly[0]                            # 3: the first parameter
    >>>
    >>> # To make a compound model consisting of a gaussian and a constant background
    >>>
    >>> gauss = GaussModel( )                   # gaussian model
    >>> gauss += PolynomialModel( 0 )           # gaussian on a constant background
    >>> print( gauss.getNumberOfParameters( ) )
    >>> 4
    >>>
    >>> # Set limits to this model
    >>>
    >>> lolim = [0,-10,0,-5]                    # lower limits for the parameters
    >>> hilim = [10,10,2, 5]                    # high limits for parameters
    >>> gauss.setLimits( lolim, hilim )         # set limits. Does not work with all Fitters
    >>>
    >>> # Pipe a model; The order of operation matters.
    >>> # m5 = ( m1 | m2 ) + m3
    >>>
    >>> m1 = PolynomialModel( 1 )               # m1( x, p )
    >>> m2 = SineModel()                        # m2( x, q )
    >>> m3 = PolynomialModel( 0 )               # m3( x, r )
    >>> m4 = m1 | m2                            # m2( m1( x, p ), q )
    >>> m5 = m4 + m3                            # m2( m1( x, p ), q ) + m3( x, r )
    >>> print( m5.parameters )                  # [p, q, r]
    >>>
    >>> # Implicit brackets
    >>> # m5 = m1 | ( m2 + m3 )
    >>>
    >>> m1 = PolynomialModel( 1 )               # m1( x, p )
    >>> m2 = SineModel()                        # m2( x, q )
    >>> m3 = PolynomialModel( 0 )               # m3( x, r )
    >>> m4 = m2 + m3                            # m2( x, q ) + m3( x, r )
    >>> m5 = m1 | m4                            # m2( m1( x, p ), q ) + m3( m1( x, p ), r )
    >>> print( m5.parameters )                  # [p, q, r]

    Attributes
    ----------
    parameters : array_like
        parameters of the model
    stdevs : None or array_like
        standard deviations after a fit to the data
    xUnit : astropy.units or list of
        unit of the x-values (list of in case of more dimensions)
    yUnit : astropy.units
        unit of the y-values
    npars : int (read only)
        number of parameters in this model
    npchain : int (read only)
        identical to npars

    Attributes from FixedModel
    --------------------------
        npmax, fixed, parlist, mlist

    Attributes from BaseModel
    --------------------------
        npbase, ndim, priors, posIndex, nonZero, tiny, deltaP, parNames

    Author       Do Kester

    """

    NOP = 0
    ADD = 1
    SUB = 2
    MUL = 3
    DIV = 4
    PIP = 5

    #  *****CONSTRUCTOR*********************************************************
    def __init__( self, nparams=0, ndim=1, copy=None, params=None, **kwargs ):
        """
        Initializes the Model with all attributes set to None, except for
        the parammeters which are all initialized to 0.

        Parameters
        ----------
        nparams : int
            the number of parameters in this model
        ndim : int
            the dimensionality of the xdatas (default: 1)
        copy : Model
            model to be copied (default: None)
        params : array_like
            initial parameters of the model
        fixed : None or dictionary of {int:float|Model}
            int         index of parameter to fix permanently.
            float|Model values for the fixed parameters.
            Attribute fixed can only be set in the constructor.
            See: @FixedModel

        """
        super( Model, self ).__init__( nparams=nparams, ndim=ndim, copy=copy,
                            **kwargs )

        setatt( self, "_next", None )
        setatt( self, "_head", self )

        if params is None :
            params = numpy.zeros( nparams, dtype=float )
        else :
            params = Tools.toArray( params, dtype=float )
        setatt( self, "parameters", self.select( params ), type=float, islist=True )

        nparams = self.npbase                       # accounting len( fixed )
        setatt( self, "_npchain", nparams )
        setatt( self, "stdevs", None )

        # xUnit is by default a (list[ndim] of) scalars, unitless
        xUnit = units.Unit( 1.0 ) if ndim == 1 else [units.Unit( 1.0 )]*ndim
        setatt( self, "xUnit", xUnit )
        setatt( self, "yUnit", units.Unit( 1.0 ) )                  # scalar
        setatt( self, "_operation", self.NOP )

        if copy is None : return

        ## the head has been copied; append now the rest of the chain
        last = copy._next
        while last is not None:
            self.appendModel( last.isolateModel( 0 ), last._operation )
            last = last._next

        setatt( self, "_npchain", copy._npchain )
        if copy.parameters is not None :
            setatt( self, "parameters", copy.parameters.copy() )
        if copy.stdevs is not None :
            setatt( self, "stdevs", copy.stdevs.copy() )

        setatt( self, "xUnit", copy.xUnit )
        setatt( self, "yUnit", copy.yUnit )

    def copy( self ):
        """ Return a copy.  """
        return Model( copy=self )

    def __setattr__( self, name, value ):
        """
        Set attributes.

        Parameters
        ----------
        name :  string
            name of the attribute
        value :
            value of the attribute

        """

        if name in ['parameters', 'stdevs'] :
            if value is not None and Tools.length( value ) != self.npchain :
                raise ValueError( "%s: Nr of %s does not comply. expect %d; got %d" %
                                ( self.shortName(), name, self.npchain, Tools.length( value ) ) )
            setatt( self, name, value, type=float, islist=True, isnone=True )

        elif name in ['xUnit', 'yUnit'] :
            isl = ( name == 'xUnit' and self.ndim == 2 )
            setatt( self, name, value, type=units.core.UnitBase, islist=isl, isnone=True )

        else :
            super( Model, self ).__setattr__( name, value )

    def __getattr__( self, name ) :
        """
        Return value belonging to attribute with name.

        Parameters
        ----------
        name : string
            name of the attribute
        """
        if name == 'npars' or name == 'npchain' :
            return self._head._npchain

        return super( Model, self ).__getattr__( name )


    def chainLength( self ):
        """ Return length of the chain.  """
        last = self
        i = 0
        while last != None:
            last = last._next
            i += 1
        return i

    def isNullModel( self ) :
        """
        Return True if the model has no parameters (a NullModel).
        """
        return len( self.parameters ) == 0

    def isolateModel( self, k ):
        """
        Return a ( isolated ) copy of the k-th model in the chain.
        Fixed parameters and priors which might be present in the compound model
        will be lost in the isolated model.

        Parameters
        ----------
        k : int
            the model number ( head = 0 )

        Raises
        ------
        IndexError when the chain is shorter than k.

        """
        last = self
        i = 0
        np = 0
        while last is not None:
            if i == k :
                next = last._next               #  save

                setatt( last, "_next", None )   #  isolate
                head = last._head               #  save
                setatt( last, "_head", last )   #  isolate
                mdl = last.copy()               #  copy
                setatt( last, "_next", next )   #  restore
                setatt( last, "_head", head )
                setatt( mdl, "_operation", self.NOP )

                n2 = np + last.npbase
                # hard overwrite of the base parameters (over the chain parameters)
                setatt( mdl, "parameters", self._head.parameters[np:n2] )

                if self.stdevs is not None :
                    setatt( mdl, "stdevs", self._head.stdevs[np:n2] )
                else :
                    setatt( mdl, "stdevs", None )
                setatt( mdl, "_npchain", mdl.npbase )      # only extract one model
                return mdl

            np += last.npbase
            last = last._next
            i += 1
        raise IndexError( "There are only " + i + " models in this compound model" )

    #  *************************************************************************
    def addModel( self, model ):
        """
        Make a compound model by concatinating/adding model to this.

        The final result is the sum of the individual results.

        The compound model is implemented as a chain of Models.
        Each of these base models contain the attributes ( parameters, limits etc. )
        and when needed these attributes are taken from there, or stored there.

        The operation (addition in this case) is always with the total result of the
        existing chain. For the use of "brackets" in a chain use BracketModel.

        Parameters
        ----------
        model : Model
            model to be added to

        """
        self.appendModel( model, self.ADD )

    def subtractModel( self, model ):
        """
        Make a compound model by concatinating/subtracting a model from this.

        The final result is the difference of the models.

        Parameters
        ----------
        model : Model
            model to be subtracted from

        """
        self.appendModel( model, self.SUB )

    def multiplyModel( self, model ):
        """
        Make a compound model by concatinating/multiplying a model with this.

        The final result is the product of the models.

        Parameters
        ----------
        model : Model
            model to be multiplied by

        """
        self.appendModel( model, self.MUL )

    def divideModel( self, model ):
        """
        Make a compound model by concatinating/dividing by a model.

        The final result is the division of the models.

        Parameters
        ----------
        model : Model
            model to be divided by

        """
        self.appendModel( model, self.DIV )

    def pipeModel( self, model ):
        """
        Make a compound model by piping the result into the next.

        Parameters
        ----------
        model : Model
            model to pipe into

        """
        self.appendModel( model, self.PIP )

    def appendModel( self, model, operation ):
        """
        Append a model to the present chain using a operation.

        Parameters
        ----------
        model : Model
            the model to be appended
        operation : int
            operation index

        Raises
        ------
        ValueError when a model of a different dimensionality is offered.
        """
        if self.isDynamic() and model.isDynamic():
            raise ValueError( "Only one Dynamic model in a chain" )
        if self.ndim != model.ndim:
            raise ValueError( "Trying to add incompatible models, with dimensions: %d and %d"%
                                (self.ndim, model.ndim) )

        if model._next is not None :
            model = Brackets( model )     # provide brackets if model is a chain

        last = self
        while last._next != None:
            last = last._next
        setatt( last, "_next",  model )
        setatt( model, "_operation", operation )
        while last._next != None:
            last = last._next
            setatt( last, "_head", self._head )

        setatt( self, "_npchain", len( self.parameters ) + len( model.parameters ) )

        setatt( self, "parameters", self._optAppend( self.parameters, model.parameters ) )

        # Erase the model's attributes; not needed anymore
        setatt( model, "parameters", None )

        return

    def _optAppend( self, x, y ) :
        if y is not None :
            if x is None : return y
            else : return numpy.append( x, y )
        else : return x

    #  *****CHECK**************************************************************
    def correctParameters( self, params ):
        """
        Check parameters for non-zero and positivity

        Parameters
        ----------
        params : array_like
            parameters for the model.

        """
        newpar = numpy.asarray( [], dtype=float )

        pp = self._recursiveCorrect( params, newpar )
        return pp

    def _recursiveCorrect( self, params, newpar ) :

        np = self.npbase
        newpar = numpy.append( newpar,  super( Model, self ).checkParameter( params[:np] ) )
        model = self._next
        if model is None :
            return newpar

        newpar = model._recursiveCorrect( params[np:], newpar )
        return newpar

    #  *****RESULT**************************************************************
    def result( self, xdata, param=None ):
        """
        Return the result of the model as applied to an array of input data.

        Parameters
        ----------
        xdata : array_like
            input data
        param : array_like
            parameters for the model. Default parameters from the Model

        """
        if param is None :
            param = self.parameters

        res = None
        xdata = Tools.toArray( xdata )
        return self._recursiveResult( xdata, param, res )

    def _recursiveResult( self, xdata, param, res ) :

        np = self.npbase

        if not self._operation == self.PIP :
            nextres = super( Model, self ).result( xdata, param[:np] )
        else :
            nextres = None
        res = self.operate( res, param[:np], nextres )
        model = self._next
        if model is None :
            return res

        res = model._recursiveResult( xdata, param[np:], res )
        return res

    def operate( self, res, pars, next ):
        if res is None or self._operation == self.NOP: # first one
            res = next
        elif self._operation == self.ADD:                 # NOP & ADD
            res = numpy.add( res, next )
        elif self._operation == self.SUB:
            res = numpy.subtract( res, next )
        elif self._operation == self.MUL:
            res = numpy.multiply( res, next )
        elif self._operation == self.DIV:
            res = numpy.divide( res, next )
        elif self._operation == self.PIP:
            res = super( Model, self ).result( res, pars )

        return res

    #  *****DERIVATIVE*********************************************************
    def derivative( self, xdata, param, useNum=False ):
        """
        Return the derivatives (df/dx) of the model at the inputs

        Parameters
        ----------
        xdata : array_like
            an input vector or array
        param : array_like
            parameters for the model
        useNum : bool
            if true, numeric derivatives are used.

        """
        result = None
        df = numpy.zeros_like( xdata )
        xdata = Tools.toArray( xdata )
        df = self._recursiveDerivative( xdata, param, result, df,
                    useNum=useNum )
        return df

    def _recursiveDerivative( self, xdata, param, result, df, useNum=False ):
        """
        Workhorse for derivative.

        Implements the optional concatenation of models via recursion

        """
        np = self.npbase
        par = param[:np]
        nextres = None

#        print( self.shortName(), self._operation )

        if self.npmax > 0:            #  the base model has no parameters at all: skip
            xd = xdata if self._operation is not self.PIP else result
            if useNum:
                nextdf = super( Model, self ).numDerivative( xd, par )
            else:
                nextdf = super( Model, self ).derivative( xd, par )

        else :
            nextdf = numpy.zeros_like( xdata, dtype=float )

        if self._operation <= self.ADD :
            df += nextdf

        elif self._operation == self.SUB :
            df -= nextdf

        elif self._operation == self.MUL :
            nextres = super( Model, self ).result( xdata, par )
            df = df * nextres + nextdf * result

        elif self._operation == self.DIV :
            nextres = super( Model, self ).result( xdata, par )
            df = ( df * nextres - nextdf * result ) / ( nextres * nextres )

        elif self._operation == self.PIP :
            nextres = 'dummy'           ## not needed, but needs something
            df *= nextdf
        else :
            raise ValueError( "Unknown operation: %d" & self.PIP )

        if self._next is None :
            return df
        if nextres is None:
            nextres = super( Model, self ).result( xdata, par )
        result = self.operate( result, par, nextres )
        #  append the dfs of the _next model

        model = self._next
        return model._recursiveDerivative( xdata, param[np:], result, df,
                    useNum=useNum )

    #  *****PARTIAL*************************************************************
    def partial( self, xdata, param, useNum=False ):
        """
        Return the partial derivatives of the model at the inputs

        Parameters
        ----------
        xdata : array_like
            an input vector or array
        param : array_like
            parameters for the model
        useNum : bool
            if true, numeric partials are used.

        """
        result = None
        partial = None
        xdata = Tools.toArray( xdata )
        partial = self._recursivePartial( xdata, param, 0, result,
                            partial, useNum=useNum )
        return partial

    def _recursivePartial( self, xdata, param, at, result, partial, useNum=False ):
        """
        Workhorse for partial.

        Implements the optional concatenation of models via recursion

        """
        np = self.npbase
        par = param[at:at+np]
        nextres = None

        if np > 0 :
            xd = xdata if self._operation is not self.PIP else result
            if useNum:
                nextpartial = super( Model, self ).numPartial( xd, par )
            else:
                nextpartial = super( Model, self ).partial( xd, par )

        else :
            inlen = Tools.length( xdata )
            nextpartial = numpy.ndarray( (inlen,0), dtype=float )

        if self._operation == self.SUB :
            nextpartial = numpy.negative( nextpartial )

        elif self._operation == self.MUL :
            nextres = super( Model, self ).result( xdata, par )
            partial = numpy.multiply( partial.transpose(), nextres ).transpose()
            nextpartial = numpy.multiply( nextpartial.transpose(), result ).transpose()

        elif self._operation == self.DIV :
            nextres = super( Model, self ).result( xdata, par )
            partial = numpy.divide( partial.transpose(), nextres ).transpose()
            invres = - result / ( nextres * nextres )
            nextpartial = numpy.multiply( nextpartial.transpose(), invres ).transpose()

        elif self._operation == self.PIP :
            nextres = 'dummy'       ## not needed, but needs something

            if useNum :
                dfdx = super( Model, self ).numDerivative( result, par )
            else :
                dfdx = super( Model, self ).derivative( result, par )
            partial = numpy.multiply( partial.transpose(), dfdx ).transpose()

        partial = ( nextpartial if partial is None
                    else numpy.append( partial, nextpartial, axis=1 ) )

        model = self._next
        if model is None:
            return partial
        if nextres is None:
            nextres = super( Model, self ).result( xdata, par )

        result = self.operate( result, param[at:], nextres )
        #  append the partials of the _next model
        at += np
        return model._recursivePartial( xdata, param, at, result, partial,
                    useNum=useNum )

    #  *****TOSTRING***********************************************************
    def __str__( self ):
        """ Returns a string representation of the model.  """
        return self._toString( "" )

    def _toString( self, indent, npars=0 ) :
        opname = [" null\n", " +\n", " -\n", " *\n", " /\n", " |\n" ]
        np = self.npbase

        if self._next is None :
            return super( Model, self )._toString( npars=npars )

        return ( super( Model, self )._toString( npars=npars ) +
                 opname[self._next._operation] +
                 indent + self._next._toString( indent, npars=npars+np ) )

    def shortName( self ) :
        """
        Return a short version the string representation: upto first non-letter.
        """
        opname = [" null\n", " + ", " - ", " * ", " / ", " | " ]
        if self._next is None :
            return super( Model, self ).shortName()

        return ( super( Model, self ).shortName( ) +
                 opname[self._next._operation] +
                 self._next.shortName( ) )


    #  *****GET/SET*************************************************************
    def getNumberOfParameters( self ):
        """ Returns the number of parameters of the ( compound ) model.  """
        return len( self._head.parameters )

    #  *****NUMERIC DERIVATIVE*************************************************
    def numDerivative( self, xdata, param ):
        """
        Returns numerical derivatives (df/dx) of the model function.

        Parameters
        ----------
        xdata : array_like
            input data
        param : array_like
            a parameters vector

        """
        return self.derivative( xdata, param, useNum=True )

    #  *****NUMERIC PARTIAL****************************************************
    def numPartial( self, xdata, param ):
        """
        Returns numerical partial derivatives of the model function.

        Parameters
        ----------
        xdata : array_like
            input data
        param : array_like
            a parameters vector

        """
        return self.partial( xdata, param, useNum=True )

    def isDynamic( self ) :
        """
        Return whether the model can change the number of parameters dynamically.
        """
        if super( Model, self ).isDynamic() :
            return True
        elif self._next is not None :
            return self._next.isDynamic()

        return False

    #  *****PRIOR***************************************************************
    def hasPriors( self, isBound=True ) :
        """
        Return True when the model has priors for all its parameters.

        Parameters
        ----------
        isBound : bool
            Also check if the prior is bound.
        """
        hasp = True
        mdl = self
        while mdl is not None and hasp :
            hasp = hasp and super( Model, self ).hasPriors( isBound=isBound )
            mdl = mdl._next
        return hasp

    def getPrior( self, k ):
        """
        Return the prior of the indicated parameter.

        Parameters
        ----------
        k : int
            parameter number.

        Raises
        ------
        IndexError when k is larger than the number of parameters.

        """
        np = self.npbase
        if k < np:
            return super( Model, self  ).getPrior( k )
        elif self._next != None:
            return self._next.getPrior( k - np )
        else:
            raise IndexError( "The (compound) model does not have %d parameters"%( k + 1 ) )

    def setPrior( self, k, prior=None, **kwargs ):
        """
        Set the prior for the indicated parameter.

        Parameters
        ----------
        k : int
            parameter number.
        prior : Prior
            prior for parameter k
        kwargs : keyword arguments
            attributes to be passed to the prior

        Raises
        ------
        IndexError when k is larger than the number of parameters.

        """
        np = self.npbase
        if k < np:
            super( Model, self  ).setPrior( k, prior=prior, **kwargs )
        elif self._next != None:
            self._next.setPrior( k - np, prior=prior, **kwargs )
        else:
            raise IndexError( "The (compound) model does not have %d parameters"%( k + 1 ) )


    #  ***PARAMETER NAME *******************************************************
    def getParameterName( self, k ):
        """
        Return the name of the indicated parameter.

        Parameters
        ----------
        k : int
            parameter number.

        Raises
        ------
        IndexError when k is larger than the number of parameters.

        """
        np = self.npbase
        if k < np:
            return super( Model, self ).getParameterName( k )
        elif self._next != None:
            return self._next.getParameterName( k - np )
        else:
            raise IndexError( "The (compound) model does not have " + str( k + 1 ) +
                              " parameters." )

    def getParameterUnit( self, k ):
        """
        Return the unit of the indicated parameter.

        Parameters
        ----------
        k : int
            parameter number.

        Raise
        -----
        IndexError when k is > number of parameters

        """
        np = self.npbase
        if k < np:
            return super( Model, self).getParameterUnit( k )
        elif self._next != None:
            return self._next.getParameterUnit( k - np )
        else:
            raise IndexError( "The (compound) model does not have " + str( k + 1 ) +
                              " parameters." )

    def getIntegralUnit( self ):
        """ Return the unit of the integral of the model over x. """

        unit = self.yUnit
        if isinstance( self.xUnit, list ) :
            for u in self.xUnit : unit *= u
        else : unit *= self.xUnit
        return unit


    #  *****LIMITS**************************************************************
    def setLimits( self, lowLimits=None, highLimits=None ):
        """
        Sets the limits for the parameters of the compound model.

        1. It is valid to insert for either parameter a None value
        indicating no lower/upper limits.
        2. When a lowerlimit >= upperlimit no limits are enforced.
        It only works in *Fitter classes which support it.

        Parameters
        ----------
        lowLimits : array_like
            lower limits on the parameters
        highLimits : array_like
            upper limits on the parameters

        """
        ml = 0
        if lowLimits is not None :
            lowLimits = Tools.toArray( lowLimits )
            nl = len( lowLimits )
            ml = max( ml, nl )
        if highLimits is not None :
            highLimits = Tools.toArray( highLimits )
            nh = len( highLimits )
            ml = max( ml, nh )
        if ml == 0 : return

        if ml > self.npchain :
            warnings.warn( "More limits given than parameters present: %d < %d" %
                            ( ml, self.npchain ) )

        for k in range( ml ) :
            lo = None if k >= nl else lowLimits[k]
            hi = None if k >= nh else highLimits[k]
            self.setPrior( k, limits=[lo,hi] )

    def getLimits( self ) :
        """
        Return the limits stored in the priors

        Returns
        -------
        limits : tuple of 2 array-like or of 2 None (if `self.priors` is None)
            (lowlimits, highlimits)

        """
        lolim = []
        hilim = []
        mdl = self
        while mdl is not None :
            if not super( Model, mdl ).hasLimits( ) :
                return [None,None]
            lolim += [p.lowLimit for p in mdl.priors]
            hilim += [p.highLimit for p in mdl.priors]

            mdl = mdl._next
        return (lolim, hilim)

    #  *************************************************************************
    def hasLimits( self, fitindex=None ):
        """
        Return true if limits has been set for this model.

        Parameters
        ----------
        fitindex    list of indices to be fitted.

        """
        haslim = True
        mdl = self
        while mdl is not None :
            haslim = haslim and super( Model, mdl ).hasLimits( fitindex=fitindex )
            if not haslim :
                return haslim
            if fitindex is not None :
                q = numpy.where( fitindex >= mdl.npbase )
                fitindex = fitindex[q] - mdl.npbase
            mdl = mdl._next
        return haslim

    def stayInLimits( self, oldpar, trypar, parlist ):
        """
        Return parameters guaranteed to be within limits.

        If out-of-limits the parameter is replaced by the midpoint between the
        old parameter and the limit.

        Parameters
        ----------
        oldpar : array_like
            old parameters ( must be inside limits )
        trypar : array_like
            parameters to be checked
        fitindex : array_like
            indices of the fitable parameters

        """
        if self.priors is None : return trypar

        i = 0
        for k in parlist :
            if self.getPrior( k ).isOutOfLimits( trypar[i] ) :
                 trypar[i] = 0.5 * ( self.getPrior( k ).stayInLimits( trypar[i] ) + oldpar[i] )
            i += 1
        return trypar

    #  *************************************************************************
    def checkLimits( self, param, parlist ):
        """
        Checks whether the parameters are within limits.

        It only works in *Fitter classes which can support it. Eg. Fitter
        itself cannot handle limits as it would turn the problem into
        a *non*-linear one.

        Parameters
        ----------
        param : list of float
            parameters to be checked
        parlist : list of int
            indices of the fitable parameters

        Raises
        ------
        ValueError when out of limits.

        """
        if self.priors is None : return

        i = 0
        for k in parlist :
            self.getPrior( k ).checkLimit( param[i] )
            i += 1

    #  ****** UNIT <--> DOMAIN ********************************************
    def unit2Domain( self, uvalue, kpar=None ):
        """
        Convert a value in [0,1] to one inside the limits of the parameter.

        Parameters
        ----------
        uvalue : (list of) float
            value in [0,1]
        kpar : int
            index of the parameter

        """
        if kpar is not None :
            return self.getPrior( kpar ).unit2Domain( uvalue )

        pgen = self.nextPrior(  )
        dval = numpy.fromiter( ( next( pgen ).unit2Domain( uv ) for uv in uvalue ), float )
        try :
            pgen.close()
        except :
            pass
        return dval


    def domain2Unit( self, dvalue, kpar=None ):
        """
        Convert a value within the domain of the parameter to one in [0,1].

        Parameters
        ----------
        dvalue : (list of) float
            value of parameter
        kpar : int
            index of the parameter

        """
        if kpar is not None :
            return self.getPrior( kpar ).domain2Unit( dvalue )

        pgen = self.nextPrior()
        uval = numpy.fromiter( ( next( pgen ).domain2Unit( dv ) for dv in dvalue ), float )
        try :
            pgen.close()
        except :
            pass
        return uval

    def partialDomain2Unit( self, dvalue ):
        """
        Return a the derivate of the domain2Unit function to dval.

        Parameters
        ----------
        dvalue : (list of) float
           parameter array

        """
        pgen = self.nextPrior(  )
        part = numpy.fromiter( ( next( pgen ).partialDomain2Unit( dv ) for dv in dvalue ), float )
        try :
            pgen.close()
        except :
            pass
        return part

    def nextPrior( self ) :
        mdl = self
        k = 0
        while True :
            try :
                yield mdl.priors[k]
            except :
                yield mdl.priors[-1]
            k += 1
            if k >= mdl.npbase :
                mdl = mdl._next
                k = 0


    #  *****SOME DUMMY METHODS TO ALLOW LINEAR MODELS IN NONLIN FITTERS********
    def isMixed( self ):
        """ Return false.  """
        return False

    def getLinearIndex( self ):
        """ Return null.  """
        return None

    #  ***** PYTHON INTERFACES ****************************************************
    def __getitem__( self, i ):
        """
        Return the i-th parameter.

        >>> p_i = model[i]

        Parameters
        ----------
        i : int
            index for the parameter.

        """
        return self.parameters[i]

    def __call__( self, x ):
        """
        Return the result of the model.

        >>> y = model( x )
        is equivalent to
        >>> y = model.result( x, model.parameters )

        Parameters
        ----------
        x : (list of) float
            apply the model to x (as input)

        """
        return self.result( x, self.parameters )

    def __iadd__( self, model ):
        """
        Method for making compound models using += operator.

        >>> model1 += model2

        Parameters
        ----------
        model : Model
            a model to add to self (the existing chain)

        """
        self.addModel( model )
        return self

    def __add__( self, model ):
        """
        Method for making compound models using + operator.

        >>> model1 = model2 + model3

        Parameters
        ----------
        model : Model
            a copy of model to add to a copy of self (the existing chain)

        """
        return self.copy().__iadd__( model.copy() )

    def __isub__( self, model ):
        """
        Method for making compound models using -= operator.

        >>> model1 -= model2

        Parameters
        ----------
        model : Model
            a model to subtract from self (the existing chain)

        """
        self.subtractModel( model )
        return self

    def __sub__( self, model ):
        """
        Method for making compound models using - operator.

        >>> model1 = model2 - model3

        Parameters
        ----------
        model : Model
            a copy of model to subtract from a copy of self (the existing chain)

        """
        return self.copy().__isub__( model.copy() )

    def __imul__( self, model ):
        """
        Method for making compound models using *= operator.

        >>> model1 *= model2

        Parameters
        ----------
        model : Model
            a model to multiply with self (the existing chain)

        """
        self.multiplyModel( model )
        return self

    def __mul__( self, model ):
        """
        Method for making compound models using * operator.

        >>> model1 = model2 * model3

        Parameters
        ----------
        model : Model
            a copy of model to multiply with a copy of self (the existing chain)

        """
        return self.copy().__imul__( model.copy() )

    def __itruediv__( self, model ):
        """
        Method for making compound models using /= operator.

        >>> model1 /= model2

        Parameters
        ----------
        model : Model
            a model to divide self with (the existing chain)

        """
        self.divideModel( model )
        return self

    def __truediv__( self, model ):
        """
        Method for making compound models using / operator.

        >>> model1 = model2 / model3

        Parameters
        ----------
        model : Model
            a copy of model to divide a copy of self with (the existing chain)

        """
        return self.copy().__itruediv__( model.copy() )

    def __ior__( self, model ):
        """
        Method for making compound models using |= operator.

        >>> model1 |= model2

        Use the results of model1 as input for model2.

        Parameters
        ----------
        model : Model
            a model to pipe the previous results through

        """
        self.pipeModel( model )
        return self

    def __or__( self, model ):
        """
        Method for making compound models using | (pipe) operator.

        >>> model1 = model2 | model3

        Parameters
        ----------
        model : Model
            a copy of model to pipe the results of a copy of self (the existing chain)

        """
        return self.copy().__ior__( model.copy() )

    #  *************************************************************************
    def testPartial( self, xdata, params ):
        """
        A test routine to check the calculation of the partial derivatives.

        It is compared to a numerical calculation.

        Parameters
        ----------
        xdata : (list of) float
            values of the independent variable
        params : list of floats
            parameters for the model

        Return
        ------
        The number of large deviations.

        """
        xdata  = Tools.toArray( xdata, ndim=self.ndim )
        params = Tools.toArray( params )
        sz = xdata.shape
        res = self.result( xdata, params )

        random.seed( 13010804 )
        try :
            df = self.derivative( xdata, params )
            numdf = self.numDerivative( xdata, params )
            hasderiv = True
        except :
            hasderiv = False

        partial = self.partial( xdata, params )
        numeric = self.numPartial( xdata, params )

        kerr = 0
        lrang = range( sz[0] )
        if sz[0] > 10 :
            lrang = random.sample( lrang, 10 )
        for k in lrang :
            print( "xdata[%2d]" % k, end='' )
            if self.ndim == 1 :
                print( " %8.3f"%(xdata[k]), end='' )
            else :
                for i in range( sz[1] ) :
                    print( " %8.3f"%(xdata[k,i]), end='' )
            print( "      result    %10.5f" % (res[k]) )

            print( "     par     value      partial   numpartial      numeric" )

            if hasderiv :
                snum = self.strictNumeric( xdata[k], params )
                print( "      df             %10.5f   %10.5f   %10.5f" %
                        ( df[k], numdf[k], snum ) )
                if ( abs( df[k] - snum ) > 0.001 or
                     abs( numdf[k] - snum ) > 0.001 ) :
                    kerr += 1

            for i in range( self.npchain ) :
                snum = self.strictNumeric( xdata[k], params, kpar=i )
                if ( abs( partial[k,i] - snum ) > 0.001 or
                     abs( numeric[k,i] - snum ) > 0.001 ) :
                    err = 1
                else :
                    err = 0
                kerr += err

                if err == 1 or random.random() < 5.0 / self.npchain :
                    print( "    %4d  %8.3f   %10.5f   %10.5f   %10.5f"%
                            (i, params[i], partial[k,i], numeric[k,i], snum) )

        return kerr

    def strictNumeric( self, x, param, kpar=None ) :
        """
        Strictly numeric calculation of derivative.

        For compound models it is different from numPartial and numDerivative.

        Parameters
        ----------
        x : float or array-like
            xdata
        param : array-like
            parameters
        kpar : None or int
            None : return derivative to x
            int  : return derivative to parameter kpar.
        """
        x = Tools.toArray( x, ndim=self.ndim )

        dx = self.deltaP[0]
        if kpar is None :
            r1 = self.result( x+dx, param )
            r2 = self.result( x-dx, param )
            return ( r1 - r2 ) / ( 2 * dx )

        pp = param.copy()
        pp[kpar] += dx
        r1 = self.result( x, pp )
        pp[kpar] -= 2 * dx
        r2 = self.result( x, pp )
        return ( r1 - r2 ) / ( 2 * dx )

##### End Model #########################################################


class Brackets( Model ):
    """
    Brackets is only for use in Model. Use BracketModel for independent uses.

    """

    #  *************************************************************************
    def __init__( self, model, copy=None, **kwargs ):

        super( Brackets, self ).__init__( model.npchain, ndim=model.ndim,
                                copy=copy, **kwargs )
        setatt( self, "model", model, type=Model )
        setatt( self, "parameters", model.parameters )

        next = model
        deep = 1
        while next is not None :
            if isinstance( next, Brackets ) :
                deep += 1
            next = next._next
        setatt( self, "deep", deep )
        setatt( self, "parNames", [] )

    def copy( self ):

        return Brackets( self.model.copy(), copy=self )


    #  *****Brackets RESULT**************************************************************
    def baseResult( self, xdata, param ):
        """
        Returns the result calculated at the xdatas.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        params : array_like
            values for the parameters.

        """
        return self.model.result( xdata, param )

    #  *****Brackets PARTIAL*************************************************************
    def basePartial( self, xdata, param, parlist=None ):
        """
        Returns the partial derivatives calculated at the xdatas.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the partials
        params : array_like
            values for the parameters.
        parlist : array_like
            Not in use

        """
        return self.model.partial( xdata, param )

    #  *****Brackets DERIVATIVE***********************************************************
    def baseDerivative( self, xdata, param ):
        """
        Returns the derivative (df/dx) calculated at the xdatas.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the partials
        params : array_like
            values for the parameters.

        """
        return self.model.derivative( xdata, param )


    def setLimits( self, lowLimits=None, highLimits=None ) :
        self.model.setLimits( lowLimits=lowLimits, highLimits=highLimits )

    def getLimits( self ) :
        return self.model.getLimits()

    def nextPrior( self ) :
        yield self.model.nextPrior()

    #  ******Brackets BASENAME*******************************************************************
    def baseName( self ):
        """ Returns a string representation of the model.  """
        indent = "  "
        return "{ " + self.model._toString( indent * self.deep ) + " }"

    def basePrior( self, k ) :
        """
        Return the prior of the indicated parameter.

        Parameters
        ----------
        k : int
            parameter number.
        """
        return self.model.getPrior( k )

    def hasPriors( self, isBound=True ) :
        """
        Return True when the model has priors for all its parameters.

        Parameters
        ----------
        isBound : bool
            Also check if the prior is bound.
        """
        return self.model.hasPriors( isBound=isBound )


    #  ******Brackets Parameter NAME*******************************************************************

    def baseParameterName( self, k ):
        """
        Return the name of the indicated parameter.

        Parameters
        ---------
        k : int
            parameter number.

        """
        return self.model.getParameterName( k )

    def baseParameterUnit( self, k ):
        """
        Return the unit of the indicated parameter.

        Parameters
        ---------
        k : int
            parameter number.

        """
        return self.model.getParameterUnit( k )


##### End Brackets #########################################################

