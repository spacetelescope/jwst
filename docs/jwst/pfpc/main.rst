Description
===========

:Classes: `jwst.pfpc.pfpc_step.PFPCStep`
:Alias: pfpc


Overview
--------
The ``pfpc`` step applies point fixed pattern corrections (PFPC) to dithered spectral data.
These corrections are intended to account for residual instrumental signatures from flat
field or flux calibration uncertainties, sampling artifacts, and residual fringes.
These signatures are stable for a fixed instrument configuration, detector location, and
calibration context.  Corrections for them are derived from observations of
stars or asteroids at standard dither positions in each instrument configuration of interest.
The derived corrections can be applied to point source observations taken at the same dither
positions and configurations, as long as a target acquisition was performed to accurately
center the source.  Applying these corrections can significantly improve the signal-to-noise
ratio for extracted point source spectra.

Wavelength-dependent correction vectors for each supported mode and dither position
are stored in a :ref:`PFPC reference file <pfpc_reffile>`.

This step is currently available for MIRI MRS exposures only.
It is incorporated into the :ref:`calwebb_spec3 <calwebb_spec3>` pipeline, after
the :ref:`pixel_replace <pixel_replace_step>` step is performed.  If a correction is
successfully applied, the step produces a set of spectra in :ref:`pfpc <pfpc>` format,
one per band.  These spectra are saved as a side product of the pipeline, and are not
propagated through any further processing steps.  If no correction can be applied,
no :ref:`pfpc <pfpc>` products are produced.

Algorithm
---------

PFPC corrections must be applied to each dither position separately, so they require
an extracted spectrum for each exposure.  The input for the step is a set of cleaned
:ref:`cal <cal>` exposures for a single source. Each exposure is processed with a
set of standard pipeline steps, configured with default recommended parameters.

For each input exposure, the correction process is:

1. Check that the exposure can be corrected. If any of the following conditions
   are not met, no further processing is performed:

   - The :ref:`PFPC reference file <pfpc_reffile>` must contain a matching correction.
   - The exposure must be a point source (``SRCTYPE = "POINT"``).
   - The exposure must have an associated target acquisition.

   It is also recommended that the input exposure is reduced in the same calibration
   context used to create the reference file. If the calibration file names stored in the
   reference file table do not match the input, a warning is issued, but processing
   will proceed.

#. Run the :ref:`cube_build <cube_build_step>` step to create a resampled spectral cube
   for each input channel and band.

#. Run the :ref:`extract_1d <extract_1d_step>` step on each cube to extract a spectrum
   for each channel and band.

#. Run the :ref:`spectral_leak <spectral_leak_step>` step on the extracted spectra at
   each dither position, to correct channel 3A spectra for light leakage from channel 1B.

#. Match a correction vector from the reference file table to the input and interpolate
   it onto the extracted wavelengths.

#. Divide the spectral flux and surface brightness and their associated errors by the
   correction vector.

After all spectra are corrected, they are averaged across all dither positions to
create one final spectrum for each band. Residual fringes are corrected with the
:func:`~jwst.residual_fringe.utils.fit_residual_fringes_1d` function. The output spectra are
stored in a format identical to the :ref:`x1d <x1d>` product, with the flux, surface brightness,
residual fringe corrected flux, and residual fringe corrected surface brightness stored in
separate columns in a binary table.


References
----------
The PFPC algorithm and reference files are based on work by K. Gordon and D. Law,
"JWST MIRI Medium Resolution Spectrometer Point-fixed Pattern Corrections:
Cleaner and Higher Signal-to-noise Spectra of Point Sources"
(`2026, AJ, 172(4), 204 <https://ui.adsabs.harvard.edu/abs/2026AJ....172..204G/abstract>`__).

Step Arguments
--------------
The ``pfpc`` step has no step-specific arguments.
