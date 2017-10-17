Description
===========

Overview
--------

The ``barshadow`` step calculates the correction to be applied to
NIRSpec MSA data for uniform sources due to the bar that separates
adjacent microshutters.  This correction is applied to multislit
data after the pathloss correction has been applied in the calspec2
pipeline.

Input details
-------------
The input data should be from after the extract_2d step, so that it contains
cutouts around each slitlet.

Algorithm
---------
The reference file contains the correction as a function of Y and wavelength
for a single open shutter (the DATA1X1 extension), and for 2 adjacent open
shutters (DATA1X3).  This allows on-the-fly construction of a model for any combination
of open and closed shutters.  The shutter configuration of a slitlet is contained
in the attribute shutter_state, which shows whether the shutters of the slitlet are open,
closed or contain the source.  Once the correction as a function of Y and wavelength is
calculated, the WCS transformation from the detector to the slit frame is used
to calculate Y and wavelength for each pixel in the cutout.  The Y values are scaled from shutter
heights to shutter spacings, and then the Y and wavelength values are interpolated
into the model to determine the correction for each pixel.

Output product
--------------
The output product has the barshadow correction attached to each slit of the multislit
datamodel in the BARSHADOW extension.
