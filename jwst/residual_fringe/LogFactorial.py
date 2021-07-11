import numpy as numpy
import math
from . import Tools

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
#  *    2006 - 2014 Do Kester, SRON (Java code)
#  *    2017        Do Kester



def logFactorial( k ):
    """
    logFactorial.  It provides the natural log of n!

    if k is float, it will be truncated to int

    Parameters
    ----------
    k : int or array_like of int
        the number(s) the factorial is wanted for.

    Return
    ------
    float : the ( natural ) log( k! ).


    Example
    -------
    >>> print( logFactorial( 0 ) )
    0
    >>> print( logFactorial( [3, 5, 10] ) )
    [1.7917594692280550, 4.7874917427820458, 15.1044125730755159]

    Author
    ------
    Do Kester, shamelessly copied from J.Skilling


    """
    logfact = [0.0000000000000000,   0.0000000000000000,   0.6931471805599453,
               1.7917594692280550,   3.1780538303479458,   4.7874917427820458,
               6.5792512120101012,   8.5251613610654147,  10.6046029027452509,
              12.8018274800814709,  15.1044125730755159,  17.5023078458738865,
              19.9872144956618882,  22.5521638531234245,  25.1912211827386834,
              27.8992713838408939,  30.6718601060806755,  33.5050734501368908,
              36.3954452080330526,  39.3398841871994946,  42.3356164607534851,
              45.3801388984769076,  48.4711813518352201,  51.6066755677643698,
              54.7847293981123187,  58.0036052229805179,  61.2617017610020014,
              64.5575386270063234,  67.8897431371815259,  71.2570389671680005,
              74.6582363488301581,  78.0922235533153071,  81.5579594561150287,
              85.0544670175815156,  88.5808275421976816,  92.1361756036870929,
              95.7196945421432019,  99.3306124547874276, 102.9681986145138097,
             106.6317602606434605, 110.3206397147573909, 114.0342117814616927,
             117.7718813997450553, 121.5330815154386244, 125.3172711493568841,
             129.1239336391272161, 132.9525750356163201, 136.8027226373263829,
             140.6739236482342790, 144.5657439463448952, 148.4777669517730487,
             152.4095925844973749, 156.3608363030787984, 160.3311282166309297,
             164.3201122631951989, 168.3274454484276816, 172.3527971391628171,
             176.3958484069973736, 180.4562914175437811, 184.5338288614495070,
             188.6281734236715977, 192.7390472878448975, 196.8661816728899794,
             201.0093163992815164, 205.1681994826411994, 209.3425867525368460,
             213.5322414945632659, 217.7369341139542200, 221.9564418191303332,
             226.1905483237275973, 230.4390435657769558, 234.7017234428182633,
             238.9783895618343195, 243.2688490029827051, 247.5729140961868779,
             251.8904022097231916, 256.2211355500095351, 260.5649409718632228,
             264.9216497985528349, 269.2910976510198680, 273.6731242856937456,
             278.0675734403661750, 282.4742926876304523, 286.8931332954270488,
             291.3239500942703444, 295.7666013507606522, 300.2209486470141542,
             304.6868567656687219, 309.1641935801469003, 313.6528299498790489,
             318.1526396202092997, 322.6634991267261512, 327.1852877037752023,
             331.7178871969284728, 336.2611819791984544, 340.8150588707990210,
             345.3794070622668642, 349.9541180407702541, 354.5390855194408459,
             359.1342053695754544]


    if Tools.length( k ) == 1 :
        return ( logfact[k] if k < 100 else
            0.9189385332046727 + ( k + 0.5 ) * math.log( k ) - k + 1.0 / ( 12.0 * k ) )

    logfactarray = numpy.asarray( logfact )

    k = numpy.asarray( k, dtype=int )

    lf = numpy.zeros( len( k ), dtype=float )
    q = numpy.where( k < 100 )[0]

    lf[q] = logfactarray[k[q]]

    b = numpy.where( k >= 100 )[0]
    x = numpy.asarray( k[b], dtype=float )
    lf[b] = 0.9189385332046727 + ( x + 0.5 ) * numpy.log( x ) - x + 1.0 / ( 12.0 * x )

    return lf


