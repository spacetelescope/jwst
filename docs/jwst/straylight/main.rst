Description
===========

:Class: `jwst.straylight.StraylightStep`
:Alias: straylight

Assumption
----------
The ``straylight`` correction is only valid for MIRI MRS data.

Overview
--------
This routine removes contamination of MIRI MRS spectral data
by the MIRI cross artifact feature produced by internal reflections
within the detector arrays.  As discussed in depth for the MIRI Imager
by A. Gáspár et al. 2021 (PASP, 133, 4504), the cross artifact manifests
as a signal extending up to hundreds of pixels along the detector column and row directions from
bright sources.  This signal has both smooth and structured components whose
profiles vary as a function of wavelength.
Although the peak intensity of the cross artifact is at
most 1% of the source intensity in Channel 1 (decreasing toward longer wavelengths),
the total integrated light in this feature can be of order 20% of the total light from a given source.


In the MIRI MRS, such a signal extending along detector rows is more disruptive
than for the MIRI imager.
Since the individual IFU slices are interleaved on the detector
and staggered in wavelength from each other, the cross artifact signal thus contaminates
non-local regions in reconstructed data cubes (both non-local on the sky and offset in wavelength
space from bright emission lines).
The purpose of this routine is thus to model the cross artifact feature in a given science exposure
and subtract it at the detector level prior to reformatting
the data into three-dimensional cubes.

At the same time, this step also ensures that the median count rate (in DN/s) in regions of the detector that
see no direct light from the sky is zero for consistency with the applied flux calibration vectors.

Algorithm
---------
The basic idea of the cross artifact correction is to convolve a given science detector image with a
kernel function that has been pre-calibrated based on observations
of isolated sources and subtract the corresponding convolved image.
As such, there are no free parameters in this step when applied to science data.

In Channel 1, the kernel function is based on engineering observations of isolated bright stars and
consists of a broad low-amplitude Lorentzian function plus two pairs
of double Gaussians.
The low-amplitude Lorentzian describes the broad wings of the kernel, and typically
has a FWHM of 100 pixels or more:

.. math::
 f_{Lor} = \frac{A_{Lor} \gamma^2}{\gamma^2 + (x - x_0)^2}

where :math:`\gamma = FWHM/2` and :math:`x_0` is the column coordinate of a given pixel.

The two double Gaussian functions describe the structured component of the profile,
in which two peaks are seen on each side of a bright spectral trace on the detector.  The relative offsets of
these Gaussians (:math:`dx`) are observed to be fixed with respect to each other, with the separation of
the secondary Gaussian from the bright trace being double the separation of the first Gaussian and both
increasing as a function of wavelength.  The widths of the Gaussians (:math:`\sigma`)
are also tied, with the secondary Gaussian
having double the width of the first.  The inner Gaussians are thus described by:

.. math::
 f_{G1} = A_{G1} exp^{\frac{- (x-x_0-dx)^2}{2 \sigma^2}}

.. math::
 f_{G3} = A_{G1} exp^{\frac{- (x-x_0+dx)^2}{2 \sigma^2}}

while the outer Gaussians are described by:

.. math::
 f_{G2} = A_{G2} exp^{\frac{- (x-x_0-2 dx)^2}{8 \sigma^2}}

.. math::
 f_{G4} = A_{G2} exp^{\frac{- (x-x_0+2 dx)^2}{8 \sigma^2}}


The best-fit parameters of these models derived from engineering data are recorded in the
:ref:`MRSXARTCORR <mrsxartcorr_reffile>` reference file and applied in a pixelwise
manner to the detector data.

The kernel functions for Channels 2 and 3 rely upon engineering observations of bright extended sources,
as the magnitude of the correction is typically too small to be visible from point sources.  These
channels use only a Lorentzian kernel with the Gaussian amplitudes set to zero as such structured components are less
obvious at these longer wavelengths.  In Channel 4 no correction appears to be necessary,
and the amplitudes of all model components are set equal to zero.
