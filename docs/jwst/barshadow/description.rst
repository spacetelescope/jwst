Description
===========

:Class: `jwst.barshadow.BarShadowStep`
:Alias: barshadow

Overview
--------
The ``barshadow`` step calculates the correction to be applied to
NIRSpec MSA data for extended sources due to the bar that separates
adjacent microshutters.  This correction is applied to MultiSlit
data after the :ref:`pathloss <pathloss_step>` correction has been applied
in the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline.

Input details
-------------
The input data must have been processed through the
:ref:`extract_2d <extract_2d_step>` step, so that cutouts have been created
for each of the slitlets used in the exposure. Hence the input must be in the
form of a `~jwst.datamodels.MultiSlitModel`.

It is also assumed that the input data have been processed through the
:ref:`srctype <srctype_step>` step, which for NIRSpec MSA exposures sets the
SRCTYPE keyword value for each slit to "POINT", "EXTENDED", or "UNKNOWN." If the
source type is "EXTENDED" or "UNKNOWN", or the SRCTYPE keyword is not present,
the default action is to treat the source as extended and apply the ``barshadow``
correction. If SRCTYPE="POINT" for a given slit, the correction is not applied.

Algorithm
---------
The step loops over all slit instances contained in the input exposure, computing
and applying the barshadow correction to each slit for which the source type has
been determined to be extended.

The :ref:`barshadow_reffile` contains the correction as a function of Y
and wavelength for a single open shutter (the DATA1X1 extension), and for 2 adjacent open
shutters (DATA1X3).  This allows on-the-fly construction of a model for any combination
of open and closed shutters.  The shutter configuration of a slitlet is contained
in the attribute shutter_state, which shows whether the shutters of the slitlet are open,
closed, or contain the source.  Once the correction as a function of Y and wavelength is
calculated, the WCS transformation from the detector to the slit frame is used
to calculate Y and wavelength for each pixel in the cutout.  The Y values are scaled from shutter
heights to shutter spacings, and then the Y and wavelength values are interpolated
into the model to determine the correction for each pixel.

Once the 2-D correction array for a slit has been computed, it is applied to the
science (SCI), error (ERR), and variance (VAR_POISSON, VAR_RNOISE, and
VAR_FLAT) data arrays of the slit.
The correction values are divided into the SCI and ERR arrays, and the square of the
correction values are divided into the variance arrays.

Output product
--------------
The output is a new copy of the input `~jwst.datamodels.MultiSlitModel`, with the
corrections applied to the slit data arrays. The 2-D correction array for each slit
is also added to the datamodel in the "BARSHADOW" extension.
