import numpy as np

"""Definitions of basic 2D field polynomial of degree n by m.
n is the order of the 2D distortion, as a function of x,y
m is the degree of the relation for th t variable.
Functions can be accessed using POLY[(n,m)], INVPOLY[(n,m)], DPOLY[(n,m)].
The inverse functions INVPOLY are only defined for <2 degree functions.
The first derivative with respect of t of POLY[] are coded as DPOLY[]."""

def POLY10(e, x, y, t):
	return e[0, 0] + t*e[1, 0]

def DPOLY10(e, x, y, t):
	return e[1, 0]

def INVPOLY10(e, x, y, d):
	return (d - e[0, 0]) / e[1, 0]

def POLY11(e, x, y, t):
	return e[0,0] + x*e[0,1] + y*e[0,2] + t*(e[1,0] + x*e[1,1] + y*e[1,2])

def DPOLY11(e, x, y, t):
	return e[1,0] + x*e[1,1] + y*e[1,2]

def INVPOLY11(e, x, y, d):
	return (d - e[0,0] - x*e[0,1] - y*e[0,2])/(e[1,0] + x*e[1,1] + y*e[1,2])

def POLY12(e, x, y, t):
	return e[0,0] + x*e[0,1] + y*e[0,2] + x**2*e[0,3] + x*y*e[0,4] + y**2*e[0,5] \
               + t*(e[1,0] + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5])

def DPOLY12(e, x, y, t):
	return e[1,0] + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5] 

def INVPOLY12(e, x, y, d):
	return (d - e[0,0] - x*e[0,1] - y*e[0,2] - x**2*e[0,3] - x*y*e[0,4] \
                 - y**2*e[0,5])/(e[1,0] + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5])

def POLY13(e, x, y, t):
	return e[0,0] + x*e[0,1] + y*e[0,2] + x**2*e[0,3] + x*y*e[0,4] + y**2*e[0,5] \
               + x**3*e[0,6] + x**2*y*e[0,7] + x*y**2*e[0,8] + y**3*e[0,9] + t*(e[1,0] \
               + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5] + x**3*e[1,6] \
               + x**2*y*e[1,7] + x*y**2*e[1,8] + y**3*e[1,9])

def DPOLY13(e, x, y, t):
	return e[1,0] + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5] \
               + x**3*e[1,6] + x**2*y*e[1,7] + x*y**2*e[1,8] + y**3*e[1,9]

def INVPOLY13(e, x, y, d):
	return (d - e[0,0] - x*e[0,1] - y*e[0,2] - x**2*e[0,3] - x*y*e[0,4] - y**2*e[0,5] \
                - x**3*e[0,6] - x**2*y*e[0,7] - x*y**2*e[0,8] - y**3*e[0,9])/(e[1,0] + x*e[1,1] \
                + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5] + x**3*e[1,6] + x**2*y*e[1,7] + x*y**2*e[1,8] + y**3*e[1,9])

def POLY20(e, x, y, t):
	return e[0,0] + t*e[1,0] + t**2*e[2,0]

def DPOLY20(e, x, y, t):
	return e[1,0] + 2*t*e[2,0]

def POLY21(e, x, y, t):
	return e[0,0] + x*e[0,1] + y*e[0,2] + t*(e[1,0] + x*e[1,1] + y*e[1,2]) + t**2*(e[2,0] + x*e[2,1] + y*e[2,2])

def DPOLY21(e, x, y, t):
	return (e[1,0] + x*e[1,1] + y*e[1,2]) + 2*t*(e[2,0] + x*e[2,1] + y*e[2,2])

def POLY22(e, x, y, t):
	return e[0,0] + x*e[0,1] + y*e[0,2] + x**2*e[0,3] + x*y*e[0,4] + y**2*e[0,5] + t*(e[1,0] + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5]) + t**2*(e[2,0] + x*e[2,1] + y*e[2,2] + x**2*e[2,3] + x*y*e[2,4] + y**2*e[2,5])

def DPOLY22(e, x, y, t):
	return (e[1,0] + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5]) + 2*t*(e[2,0] + x*e[2,1] + y*e[2,2] + x**2*e[2,3] + x*y*e[2,4] + y**2*e[2,5])

def POLY23(e, x, y, t):
	return e[0,0] + x*e[0,1] + y*e[0,2] + x**2*e[0,3] + x*y*e[0,4] + y**2*e[0,5] + x**3*e[0,6] + x**2*y*e[0,7] + x*y**2*e[0,8] + y**3*e[0,9] + t*(e[1,0] + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5] + x**3*e[1,6] + x**2*y*e[1,7] + x*y**2*e[1,8] + y**3*e[1,9]) + t**2*(e[2,0] + x*e[2,1] + y*e[2,2] + x**2*e[2,3] + x*y*e[2,4] + y**2*e[2,5] + x**3*e[2,6] + x**2*y*e[2,7] + x*y**2*e[2,8] + y**3*e[2,9])

def DPOLY23(e, x, y, t):
	return (e[1,0] + x*e[1,1] + y*e[1,2] + x**2*e[1,3] + x*y*e[1,4] + y**2*e[1,5] + x**3*e[1,6] + x**2*y*e[1,7] + x*y**2*e[1,8] + y**3*e[1,9]) + 2*t*(e[2,0] + x*e[2,1] + y*e[2,2] + x**2*e[2,3] + x*y*e[2,4] + y**2*e[2,5] + x**3*e[2,6] + x**2*y*e[2,7] + x*y**2*e[2,8] + y**3*e[2,9])


POLY = {}
DPOLY = {}
INVPOLY = {}

# A Dictionary containing the functions we know, labeled by the size of the parameter array e
POLY[(2,1)] = POLY10
POLY[(2,3)] = POLY11
POLY[(2,6)] = POLY12
POLY[(2,10)] = POLY12

POLY[(3,1)] = POLY20
POLY[(3,3)] = POLY21
POLY[(3,6)] = POLY22
POLY[(3,10)] = POLY23

DPOLY[(2,1)] = DPOLY10
DPOLY[(2,3)] = DPOLY11
DPOLY[(2,6)] = DPOLY12
DPOLY[(2,10)] = DPOLY12

DPOLY[(3,1)] = DPOLY20
DPOLY[(3,3)] = DPOLY21
DPOLY[(3,6)] = DPOLY22
DPOLY[(3,10)] = DPOLY23

INVPOLY[(2,1)] = INVPOLY10
INVPOLY[(2,3)] = INVPOLY11
INVPOLY[(2,6)] = INVPOLY12
INVPOLY[(2,10)] = INVPOLY12


def npol(m):
	return int((m + 1)**2/2. + (m + 1)/2.)
