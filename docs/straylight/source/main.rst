
Description
===========

Assumption
----------
The current straylight correction is only valid for MIRI MRS Short
wavelength data. In the future this  routine may also remove and/or 
otherwise correct for stray-light that may contaminate a MIRI LRS spectrum, 
due to a bright source in the imager part of the FOV. The current 
implementation for MIRI LRS data is a no-op  (nothing is performed), 
because no known algorithm exists yet for performing this correction.

Overview
--------
This routine removes and/or otherwise corrects for stray-light that may
contaminate a MIRI MRS short-wavelength spectrum, due a bright source 
in the MRS slice 
gaps. The current routine determines the stray-light by using signal
in-between slices and linearly interpolates over the slice. 
 
The source of the MIRI MRS stray-light has been identified as being caused 
by scattering in optical components within the SMO. The stray-light is 
manifested as a signal that extends in the detector row direction. Its 
magnitude is proportional to that of bright illuminated regions of the 
spectral image, at a ratio that falls with increasing wavelength, 
from about 2 % in Channel 1A to undetectably low levels longward of Channel 2B.

Algorithm
---------

The MIRI MRS algorithm uses a stray-light MASK reference file to determine
which pixels are science pixels and which pixels fall in-between the
slices. Each illuminated pixel on the array is has a signal that is the
sum of direct illumination and the scattering from neighboring areas.
Only the pixels located between the slices are areas of indirect illumination.
The  illumination on the inter-slice pixels are used  to determine a 
stray-light component to subtract from each science pixel. 
