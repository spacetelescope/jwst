Description
===========

Assumption
----------
The ``straylight`` correction is only valid for MIRI MRS channel 1 and 2 data.
The step uses information in the REGIONS reference file about which pixels
belong to a slice and which pixels are located in the slice gaps.

Overview
--------
This routine removes stray light that may contaminate MIRI MRS channel 1/2
spectra by interpolating the measured signal in the inter-slice
regions of the detector. The inter-slice regions nominally should not receive
light from the sky and therefore should serve as a good measure of the
amount of stray light within the exposure.

The chief source of the MIRI MRS stray light appears to be caused
by scattering in the optical components within the SMO. The stray light is
manifested as a signal that extends in the detector row direction. Its
magnitude is proportional to that of bright illuminated regions of the
spectral image, at a ratio that falls with increasing wavelength,
from about 1% in Channel 1A to undetectably low levels longward of Channel 2B.
Two components of the stray light have been observed: a smooth and a structured
distribution. 

Algorithm
---------
The basic idea of the stray light removal algorithm is to only deal with the 
smooth component of the stray light. Due to the extended nature of the
stray light we use the detected signal in the slice gaps, where nominally no photons
should hit the detectors, and assume that all detected light is the stray light. 
Using this measurement, we can interpolate the gap flux within the slice to
estimate the amount of the stray light in the slice. 

There are two possible algorithms in the stray light step, both of which use the
REGIONS reference file (20% throughput threshold) to define slice and inter-slice pixels.

The first algorithm is a simplistic approach to deal with the stray light estimation row-by-row
and interpolate the gap flux linearly. An intermediate stray light map is 
generated row-by-row and then this map is further smoothed to remove row-by-row
variations. 

.. _msm_equations:

Given the extended nature of the smooth component of the MRS stray light, it
is obvious that a row-by-row handling of the stray light could be replaced
by a two-dimensional approach, so that no additional smoothing is required.
For the second algorithm we use the Modified Shepard's
Method to interpolate the gap fluxes two dimensionally. This approach has the benefit
of being more robust to outliers that can bias the first algorithm.  The stray light correction
for each science pixel is based on the flux of the gap pixels with a "region of influence"
from the science pixel. The algorithm takes each science pixel and determines the 
amount of stray light :math:`s` to remove from the pixel by interpolating the fluxes
:math:`p_i` measured by the gap pixels. The gap pixel flux is weighted by the distance
:math:`d_i` between the science pixel and gap pixel. 
The Modified Shepard’s Method uses this distance to weight the different contributors according
the equations

.. math::
 s = \frac{ \sum_{i=1}^n p_i w_i}{\sum_{i=1}^n w_i}

and

.. math::
 w_i =\frac{ max(0,R-d_i)} {R d_i}^ k

The radius of influence :math:`R` and the exponent :math:`k` are variables that 
can be adjusted to the actual problem. The default values for these parameters are
:math:`R = 50` pixels and :math:`k = 1`.
