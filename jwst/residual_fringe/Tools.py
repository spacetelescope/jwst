from __future__ import print_function

import numpy as numpy
import math as math
import numbers
import sys
import trace
import re

from astropy.table import Table

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
#  *    2016 - 2020 Do Kester


def getItem( ilist, k ) :
    """
    Return the k-th item of the ilist
        or the last when not enough
        or the ilist itself when it is not a list

    Parameters
    ----------
    ilist : an item or a list of items
        List to obtain an item from
    k : int
        item to be returned

    """
    return ( ilist if not isinstance( ilist, list ) else
             ilist[-1] if k >= len( ilist ) else ilist[k] )


def getColumnData( xdata, kcol ) :
    """
    Return the kcol-th column from xdata

    Parameters
    ----------
    xdata   2D array_like or Table
        the data array
    kcol    int
        column index
    """
    if isinstance( xdata, Table ) :
        return xdata.columns[kcol].data
    elif xdata.ndim == 2 :
        return xdata[:,kcol]
    else :
        return xdata

def isBetween( xs, x, xe ) :
    """
    Return True when x falls between xs and xe or on xs or xe.
        where the order of xs, xe is unknown.
    """
    return ( xs <= x <= xe ) or ( xe <= x <= xs )


def setAttribute( obj, name, value, type=None, islist=False, isnone=False ) :

    if type is None :                       ## no check on type
        object.__setattr__( obj, name, value )
        return True

    if isnone and value is None :           ## could be None
        object.__setattr__( obj, name, value )
        return True

    if islist :                             ## should be list of type
        isl = isList( value, type )
#        print( "TS  ", name, value, type, isl )
        if ( isl[0] ) :
            if type is not int and type is not float :
                array = value if isl[1] else [value]
            elif ( isl[1] ) : array = numpy.asarray( value, dtype=type )
            else : array = numpy.asarray( [value], dtype=type )
            object.__setattr__( obj, name, array )
            return True
        else :
            raise TypeError( name + ' has not the proper type: (list of) ' + str( type ) )
            return False

    if isInstance( value, type ) :          ## should be singular type
        object.__setattr__( obj, name, value )
        return True
    else :
        raise TypeError( name + ' has not the proper type: ' + str( type ) )
    return False


def setNoneAttributes( obj, name, value, listNone ) :
    """
    Set attribute contained in dictionary dictList into the attr-list.
    A list is a native list or a numpy.ndarray. It also checks the type.
    if values is a singular item of the proper type it will be inserted as [value].

    Parameters
    ----------
    obj : object
        to place the attribute in
    name : str
        of the attribute
    value : any
        of the attribute
    listNone : list of names
        that could have a None value

    Returns
    -------
        True on succesful insertion. False otherwise.

    Raises
    ------
        TypeError   if the type is not as in the dictionary
     """
    if name in listNone and value is None :
        object.__setattr__( obj, name, value )
        return True
    return False


def setListOfAttributes( obj, name, value, dictList ) :
    """
    Set attribute contained in dictionary dictList into the attr-list.
    A list is a native list or a numpy.ndarray. It also checks the type.
    if values is a singular item of the proper type it will be inserted as [value].

    Parameters
    ----------
    obj : object
        to place the attribute in
    name : str
        of the attribute
    value : any
        of the attribute
    dictList : dictionary
        of possible attributes {"name": type}

    Returns
    -------
        True on succesful insertion. False otherwise.

    Raises
    ------
        TypeError   if the type is not as in the dictionary

    """
    if name in dictList :
        _type = dictList[name]
        isl = isList( value, _type )
#        print( name, isl )
        if ( isl[0] ) :
            if _type is not int and _type is not float :
                _array = value if isl[1] else [value]
            elif ( isl[1] ) : _array = numpy.asarray( value, dtype=_type )
            else : _array = numpy.asarray( [value], dtype=_type )
            object.__setattr__( obj, name, _array )
            return True
        else :
            raise TypeError( name + ' has not the proper type: (list of) ' + str( dictList[name] ) )
    return False


def setSingleAttributes( obj, name, value, dictSingle ) :
    """
    Set a singular attribute contained in dictionary dictSingle into the attr-list.
    It also checks the type.

    Parameters
    ----------
    obj : object
        to place the attribute in
    name : str
        of the attribute
    value : any
        of the attribute
    dictSingle : dictionary
        of possible attributes {"name": type}

    Returns
    -------
        True on succesful insertion. False otherwise.

    Raises
    ------
        TypeError   if the type is not as in the dictionary
     """
    if name in dictSingle :
        if isInstance( value, dictSingle[name] ) :
            object.__setattr__( obj, name, value )
            return True
        else :
            raise TypeError( name + ' has not the proper type: ' + str( dictSingle[name] ) )
    return False

def makeNext( x, k ) :
    """
    Return next item of x, and last item if x is exhausted.
    Or x itself if x is singular.
    """
    try :
        _xnext = x[-1]
    except :
        _xnext = x
    while True :
        try :
            _xnext = x[k]
            yield _xnext
            k += 1
        except :
            yield _xnext


def length( x ) :
    """
    Return the length of any item. Singletons have length 1; None has length 0..
    """
    if x is None :
        return 0
    try :
        return len( x )
    except :
        return 1

def toArray( x, ndim=1, dtype=None ) :
    """
    Return a array of x when x is a number

    Parameters
    ----------
    x : any number, list/array of numbers or []
        to be converted to numpy.ndarray
    ndim : int
        minimum number of dimensions present
    dtype : type
        conversion to type (None : as is)

    """
    if isinstance( x, Table ) :
        return x
    return numpy.array( x, dtype=dtype, copy=False, ndmin=ndim )

def isList( item, cls ) :
    """
    Return (True,False) if item is a instance of cls
           (True,True)  if item is a (list|ndarray) of instances of cls
           (False,False) if not
    """
#    print( item, cls, item.__class__ )
    if isInstance( item, cls ) : return (True,False)
    islst = isinstance( item, list ) or isinstance( item, numpy.ndarray )
    if islst :
        for i in numpy.asarray( item ).flat :
            islst = islst and isInstance( i, cls )
    return (islst,islst)

def isInstance( item, cls ) :
    """
    Returns true when one of the following is true
    1. when cls is int   : item is an int or item is a numpy.integer.
    2. when cls is float : item is an float or item is an int.
    3. when cls is cls   : item is a cls.
    """
    if cls is int :
        return isinstance( item, int ) or isinstance( item, numpy.integer )
    elif cls is float :
        return isinstance( item, float ) or isInstance( item, int )
    else :
        return isinstance( item, cls )

def ndprint( x, form='{0:.3f}' ) :
    """
    Print a ndarray, formatted.
    """
    print( numpy.array2string(x, formatter={'float_kind':form.format}) )

def decorate( src, des, copy=True ) :
    """
    Transfer attributes from src to des.
    If copy is True try to copy the attributes, otherwise link it.

    Parameters
    ----------
    src : object
        source of the attributes
    des : object
        destiny for the attributes
    copy : bool
        if True: copy
    """
    atr = vars( src )
    ld = list( atr.keys() )
    for key in ld :
        value = atr[key]
        if copy :
            try :
                value = value.copy()
            except :
                pass
        object.__setattr__( des, key, value )


def printclass( cls, nitems=8 ) :
    """
    Print the attributes of a class.
    """
    numpy.set_printoptions( precision=3, threshold=10, edgeitems=4 )

    print( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
    print( cls, " at ", id( cls ) )
    print( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
    atr = vars( cls )
    ld = list( atr.keys() )
    ld.sort()
    for key in ld :
        print( "%-16.16s"%key, end="" )
        val = atr[key]
        if isinstance( val, (list,numpy.ndarray) ) :
            printlist( val )
        elif isinstance( val, str ) :
            print( shortName( val ) )
        elif key == "model" :
            print( val.shortName( ) )
        else :
            print( val )

def printlist( val, nitems=8 ) :
    nv = length( val )
    if nv == 0 or isinstance( val[0], numbers.Number ) :
        print( val, nv )
        return

    sep = "["
    for k in range( min( nv, nitems ) ) :
        try :
            print( "%s%s"%(sep, val[k].__name__ ), end="" )
        except :
            print( "%s%s"%(sep, shortName( str( val[k] ) ) ), end="" )
        sep = " "
    print( "%s"%("... ]" if nitems < nv else "]"), nv )

def shortName( val ):
    """
    Return a short version the string representation: upto first non-letter.
    """
    m = re.match( "^[0-9a-zA-Z_]*", val )
    return m[0]


def nicenumber( x ) :
    """
    Return a nice number close to x.
    """
    if x == 0 :
        return x
    n = 1.0
    while x > 10 :
        x /= 10
        n *= 10
    while x < 1 :
        x *= 10
        n /= 10
    return int( x ) * n

def fix2int( x ) :
    """
    Return integer array with values as in x
    Parameters
    ----------
    x : array_like
        array of integer floats
    """
    return numpy.fix( x + 0.000001 ).astype( int )

def track( statement ) :
    """

    Parameters
    ----------
    statement : str
        statement to be traced

    """
    # create a Trace object, telling it what to ignore, and whether to
    # do tracing or line-counting or both.
    tracer = trace.Trace( ignoredirs=[sys.prefix, sys.exec_prefix], trace=0, count=1 )

    # run the new command using the given tracer
    tracer.run( statement )

    # make a report, placing output in the current directory
    r = tracer.results()
    r.write_results(show_missing=True, coverdir=".")

def average( xx, weights=None, circular=None ) :
    """
    Return (weighted) average and standard deviation of input array.

    Parameters
    ----------
    xx : array_like
        input to be averaged
    weights : array_like
        if present these are the weights
    circular : list of 2 floats
        the input is a circular item between [low,high]
    """
    if circular is None :
        if weights is None :
            weights = numpy.ones_like( xx )
        sw = numpy.sum( weights )
        xw = xx * weights
        sx = numpy.sum( xw )
        s2 = numpy.sum( xw * xx )
        averx = sx / sw

        rr = xx - averx
        stdvx = math.sqrt( numpy.average( rr * rr, weights=weights ) )
#        stdvx = math.sqrt( s2 / sw - averx * averx )

    else :
        range = circular[1] - circular[0]
        d2r = 2 * math.pi / range
        rr = ( xx - circular[0] ) * d2r
        asx = numpy.average( numpy.sin( rr ), weights=weights )
        acx = numpy.average( numpy.cos( rr ), weights=weights )
        averx = math.atan2( asx, acx )

        ## make rr relative to average and put 0.0 in mid of range.
        rr = ( rr - averx + math.pi ) % ( 2 * math.pi ) - math.pi
        stdvx = math.sqrt( numpy.average( rr * rr, weights=weights ) ) / d2r

        ## put the average inside the circular domain
        averx = averx / d2r + circular[0]
        while averx < circular[0] : averx += range
        while averx > circular[1] : averx -= range

    return( averx, stdvx )
