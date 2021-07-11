import numpy as numpy
from astropy import units
import re
import warnings
from . import Tools
from .Tools import setAttribute as setatt

from .Prior import Prior
from .UniformPrior import UniformPrior

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
#  *    2011 - 2014 Do Kester, SRON (Java code)
#  *    2016 - 2020 Do Kester

class BaseModel( object ):
    """
    BaseModel implements the common parts of simple Models.

    A simple model is a function in x with parameters p : f(x:p).
    The variable x is an array of points in a space of one or more
    dimensions; p can have any length, including 0 i.e. the model
    has no parameters.

    The result of the function for certain x and p is given by
    model.result( x, p )

    The partial derivatives of f to p (df/dp) is given by
    model.partial( x, p )
    Some fitters make use of the partials

    The derivative of f to x (df/dx) is given by
    model.derivative( x, p )

    BaseModel checks parameters for positivity and nonzero-ness, if such
    is indicated in the model itself.

    BaseModel also implements the numerical calculation of the (partial)
    derivatives to be used when they are not given in the model definition
    itself.

    Attributes
    ----------
    npbase : int
        number of params in the base model
    ndim : int
        number of dimensions (parallel streams) of input. (default : 1)
    priors : list of Prior
        pertaining to each of the parameters of the model.
        If the list is shorter than the number of parameters, the last one is repeated.
    posIndex : list of int
        list of indices indication positive-definite parameters.
    nonZero : list of int
        list of parameters that need a warning when they are equal to zero.
        Warnings will only be issued once. Values are replaced by self.tiny
    tiny : float
        very small value, replacing zero valued when found on NonZero.
        (default : 1e-20)
    deltaP : array_like
        (list of) width(s) for numerical partial calculation. (default : 0.00001)
    parNames : list of str
        list of parameter names. (default : "parameter_k")

    Author :         Do Kester

    """

    #  *************************************************************************
    def __init__( self, nparams=0, ndim=1, copy=None, posIndex=[], nonZero=[], **kwargs ):
        """
        BaseModel Constructor.
        <br>
        Parameters
        ----------
        nparams : int
            Number of parameters in the model (default: 0)
        ndim : int
            Number of dimensions of the input (default: 1)
        copy : BaseModel
            to be copied
        posIndex : list of int
            indices of parameters that need to be > 0
        nonZero : list of int
            indices of parameters that cannot be zero.
            they are replaced by self.tiny
        kwargs
            for internal use.
        """
        super( BaseModel, self ).__init__( )

        if copy is None :
            setatt( self, "npbase", nparams )
            setatt( self, "ndim", ndim )
            setatt( self, "priors", None )
            setatt( self, "posIndex", numpy.asarray( posIndex, dtype=int ) )
            setatt( self, "nonZero", numpy.asarray( nonZero, dtype=int ) )
            setatt( self, "deltaP", numpy.asarray( [0.00001] ) )
            setatt( self, "tiny", 1e-20 )
            parNames = ["parameter_%d"%k for k in range( self.npbase )]
            setatt( self, "parNames", parNames )
        else :
            setatt( self, "npbase", copy.npbase )
            setatt( self, "ndim", copy.ndim )
            setatt( self, "priors", None if copy.priors is None else copy.priors.copy() )
            setatt( self, "posIndex", copy.posIndex )
            setatt( self, "nonZero", copy.nonZero )
            setatt( self, "deltaP", copy.deltaP )
            setatt( self, "tiny", copy.tiny )
            setatt( self, "parNames", copy.parNames )

    def __setattr__( self, name, value ):
        """
        Set attributes.

        """
        if name == "priors" :
            setatt( self, name, value, type=Prior, isnone=True, islist=True )
            return
        if name == "parNames" :
            setatt( self, name, value, type=str, islist=True )
            return
        if name == "deltaP" :
            setatt( self, name, value, type=float, islist=True )
            return
        if name == "tiny" :
            setatt( self, name, value, type=float )
            return
        if name in ["posIndex", "nonZero"] :
            setatt( self, name, value, type=int, islist=True )
            return

        raise AttributeError(
            "Model has no attribute " + name + " of type " + str( value.__class__ ) )

    def __getattr__( self, name ) :
        """
        Return value belonging to attribute with name.

        Parameters
        ----------
        name : string
            name of the attribute
        """
        if name == "priors" :
            return None

        elif name == "parNames" :
#        if name == "parNames" :
            raise AttributeError( "Model has no defined parameter names." )

        for k,pn in enumerate( self.parNames ) :
            if name == pn :
                return self.parameters[k]

        ## Raise error for any remaining attributes, not found standardly
        raise AttributeError( "Model has no attribute %s." % name )


    #  *****RESULT**************************************************************
    def result( self, xdata, param ):
        """
        Returns the result calculated at the xdatas.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        param : array_like
            values for the parameters.

        """
        self.checkParameter( param )
        return self.baseResult( xdata, param )

    #  *****PARTIAL*************************************************************
    def partial( self, xdata, param, parlist=None ):
        """
        Returns the partial derivatives calculated at the inputs.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the result
        param : array_like
            values for the parameters.
        parlist : None or array_like
            indices of active parameters

        """
        self.checkParameter( param )
        return self.basePartial( xdata, param, parlist=parlist )

    def checkParameter( self, param ) :
        """
        Return parameters corrected for positivity and Non-zero.

        Parameters
        ----------
        param : array_like
            values for the parameters.

        """
        self.checkZeroParameter( param )
        self.checkPositive( param )
        return param

    def checkPositive( self, param ) :
        """
        Check parameters for positivity. Silently correct.

        Parameters
        ----------
        params : array_like
            values for the parameters

        """
        for k in self.posIndex :
            param[k] = abs( param[k] )

    def checkZeroParameter( self, param ):
        """
        Check parameters for Non-zero. Correct after one warning.

        Parameters
        ----------
        params : array_like
            values for the parameters

        """
        for k in self.nonZero :
            if param[k] == 0 :
                msg = ( ( self.shortName() + ": " + self.baseParameterName( k ) +
                            " ( =parameter[%d] ) equals zero."%k ) )
                warnings.warn( msg )
                param[k] = self.tiny

    def isDynamic( self ) :
        """
        Whether the model implements Dynamic
        """
        return False

    def isModifiable( self ) :
        """
        Whether the model implements Modifiable
        """
        return False

    #  *****TOSTRING***********************************************************
    def __str__( self ):
        """ Returns a string representation of the model.  """
        return self.baseName( )

    def shortName( self ):
        """
        Return a short version the string representation: upto first non-letter.

        """
        m = re.match( "^[a-zA-Z_]*", self.baseName() )
        return m.group(0)

    def derivative( self, xdata, param ) :
        """
        Returns the derivative of the model to xdata.

        It is a numeric derivative as the analytic derivative is not present
        in the model.

        Parameters
        ----------
        xdata : array_like
            values at which to calculate the derivative
        param : array_like
            values for the parameters. (default: model.parameters)

        """
        return self.baseDerivative( xdata, param )

    def setPrior( self, k, prior=None, **kwargs ) :
        """
        set the prior and/or limits for the indicated parameter.

        The prior (by default UniformPrior) is appended when k is equal to np, the
        length of the existing list of priors. It replaces the prior when k < np
        and it generates an error when k > np

        Parameters
        ---------
        k : int
            parameter number.
        prior : Prior
            prior for the parameter
        kwargs : dict
            attributes to be passed to the prior

        """
        np = Tools.length( self.priors )
        assert k <= np, "The par number %d is larger than the length of the priors %d"%(k,np)


        if prior is None :
            prior = UniformPrior( ) if np == 0 else self.basePrior( k ).copy()

        prior.setAttributes( **kwargs )

        if k == np :
            self.priors = [prior] if self.priors is None else self.priors + [prior]
            return

        self.priors[k] = prior

        return

    def hasPriors( self, isBound=True ) :
        """
        Return True when the model has priors for all its parameters.

        Parameters
        ----------
        isBound : bool
            Also check if the prior is bound.
        """
        return ( self.priors is not None and ( ( not isBound ) or
                 all( [p.isBound() for p in self.priors] ) ) )

    def getPrior( self, k ) :
        """
        Return the prior of the indicated parameter.

        Parameters
        ---------
        k : int
            parameter number.
        """
        return self.basePrior( k )

    def basePrior( self, k ) :
        """
        Return the prior of the indicated parameter.

        Parameters
        ----------
        k : int
            parameter number.
        """
        np = Tools.length( self.priors )
        if np == 0 :
            raise IndexError( "The model does not have priors." )

        if k < np:
            return self.priors[k]
        else :
            return self.priors[-1]

    def getParameterName( self, k ) :
        """
        Return the name of the indicated parameter.

        Parameters
        ---------
        k : int
            parameter number.
        """
        return self.baseParameterName( k )


    def baseParameterName( self, k ) :
        """
        Return the name of the indicated parameter.

        Parameters
        ---------
        k : int
            parameter number.
        """
        return self.parNames[k]

    def getParameterUnit( self, k ) :
        """
        Return the unit of the indicated parameter.

        Parameters
        ---------
        k : int
            parameter number.
        """
        return self.baseParameterUnit( k )

    def baseParameterUnit( self, k ) :
        """
        Return the name of the indicated parameter.

        Parameters
        ---------
        k : int
            parameter number.
        """
        return units.Unit( 1.0 )

    def hasLimits( self, fitindex=None ) :
        """
        Return True if the model has limits set.
        """
        if self.npbase == 0 :
            return True
        if self.priors is None :
            return False
        if fitindex is None :
            return all( p.hasLimits() for p in self.priors )
        else :
            haslim = True
            npr = len( self.priors )
            for k in fitindex :
                if k >= npr :
                    return haslim
                haslim = haslim and self.priors[k].hasLimits()
            return haslim

