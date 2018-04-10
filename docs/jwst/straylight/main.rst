
Description
===========

Assumption
----------
The current stray-light correction is only valid for MIRI MRS Short
wavelength data. 

Overview
--------
This routine removes and/or otherwise corrects for stray-light that may
contaminate a MIRI MRS short-wavelength spectrum, due a bright source
in the MRS slice
gaps. The current routine determines the stray-light by using signal
in-between slices and linearly interpolates over the slice.

The chief source of the MIRI MRS stray-light appears to be  caused
by scattering in optical components within the SMO. The stray-light is
manifested as a signal that extends in the detector row direction. Its
magnitude is proportional to that of bright illuminated regions of the
spectral image, at a ratio that falls with increasing wavelength,
from about 1 % in Channel 1A to undetectable low levels long-ward of Channel 2B.
Two components of the stray-light have been observed, a smooth and a structured
distribution. 

Algorithm
---------
The basic idea of the stray-light removal algorithm is to deal with the 
smooth component of the stray-light only. Due to the extended nature of the
stray-light we use the detected signal in the slice gaps, where nominally no photons
should it the detectors, and assume that all detected light is the stray-light. 
Using this measurement, we can interpolate the gap flux within the slice to
estimate the amount of the stray-light in the slice. 

There are two possible algorithms in the stray-light step. The first algorithm is
a more simplistic approach by dealing with the stray-light estimation row-by-row
and interpolating the gap flux linearly. An intermediate stray-light map is 
generated row-by-row and then this map is further smoothed to remove row-by-row
variations. This algorithm uses a stray-light mask reference file that contains
1s for gap pixels and 0s for science pixels.  

Given the extended nature of the smooth component of the MRS stray-light, it
is obvious that a row-by-row handling of the stray-light could be replaced
by a two-dimensional approach such that not additional smoothing is required.
For the second algorithm we improve the technique by using the Modified Shepard's
Method to interpolate the gap fluxes two dimensionally. 

The second algorithm uses the MRS distortion reference file. The first extension
of this file contains the "Slice_Number" and consists of a detector pixel map 
indicating the slice numbers or in case of a non-science pixel, the value zero. 
Consequently, the gap pixels can be easily found by masking all pixels that are 
non-science.The assign_wcs step takes the "Slice_Number" information and stores it in
a regions file.  This second algorithm takes each slice pixel and determines the 
amount of stray-light :math:`s`  by interpolating the fluxes :math:`p_i` measured
by the gap pixels, taking the distance :math:`d_i` to the slice pixel into account. 
The Modified Shepard’s Method uses this distance to weight the different
 contributors according the equation:

.. math::

   s = \frac{ \sum_{i=1}^n p_i w_i}{\sum_{i=1}^n w_i}

where,
.. math::

   w_i =\frac{ max(0,R-d_i)} {R d_i}^ k

The radius of influence :math: 'R' and the exponent :math: 'k' are variables that 
can be adjusted to the actual problem. The default values for these parameters are
:math:`R = 50` pixels and :math:`k = 1`.