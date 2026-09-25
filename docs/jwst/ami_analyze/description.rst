Description
===========

:Class: `jwst.ami.ami_analyze_step.AmiAnalyzeStep`
:Alias: ami_analyze

The ``ami_analyze`` step is one of the AMI-specific steps in the ``ami``
sub-package that is part of Stage 3 :ref:`calwebb_ami3 <calwebb_ami3>`
processing. It applies the Lacour-Greenbaum (LG) image plane
modeling algorithm to a NIRISS AMI image.
The routine computes a number of parameters, including a model fit (and
residuals) to the image, fringe amplitudes and phases, and closure phases
and amplitudes.

The JWST AMI observing template allows for exposures to be obtained using
either full-frame (SUBARRAY="FULL") or subarray (SUBARRAY="SUB80") readouts.
When processing a full-frame exposure, the ``ami_analyze`` step extracts
and processes a region from the image corresponding to the size and location of
the SUB80 subarray, in order to reduce execution time.

Inputs
------

3D calibrated image
^^^^^^^^^^^^^^^^^^^
:Data model: `~stdatamodels.jwst.datamodels.JwstDataModel`
:File suffix: _calints

The ``ami_analyze`` step takes a single calibrated image cube as input, which should be
the "_calints" product resulting from :ref:`calwebb_image2 <calwebb_image2>` processing.
Multiple exposures can be processed via use of an ASN file that is used as input
to the :ref:`calwebb_ami3 <calwebb_ami3>` pipeline.

.. note::

   The ``ami_analyze`` step will also
   accept a 2D "_cal" product but errors will not be computed in the output.
   The ``ami_analyze`` step itself does not accept an ASN as input.

Outputs
-------

The ``ami_analyze`` step produces three output files. The first two (``_ami-oi.fits`` and ``_amimulti-oi.fits``) contain the interferometric observables, and the third (``_amilg.fits``) contains the data, LG model, and residuals. All products contain a ``PRIMARY`` HDU containing header keywords but no science data, as well as an ``ASDF`` extension. The files are described in more detail below.

The output file name syntax is exposure-based, using the input file name as the root, with
the addition of the association candidate ID and the "_ami-oi", "_amimulti-oi", or "amilg" product type suffix, e.g.
"jw87600027001_02101_00002_nis_a3001_ami-oi.fits."

Interferometric observables
^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~stdatamodels.jwst.datamodels.AmiOIModel`
:File suffix: _ami-oi.fits, _amimulti-oi.fits

The inteferometric observables are saved as OIFITS files, a registered FITS format
for optical interferometry, containing the following list of extensions:

1. ``OI_ARRAY``: AMI subaperture information
2. ``OI_TARGET``: target properties
3. ``OI_T3``: extracted triple-product amplitudes, closure phases
4. ``OI_VIS``: extracted visibility (fringe) amplitudes, phases
5. ``OI_VIS2``: squared visibility (fringe) amplitudes
6. ``OI_Q4``: closure amplitudes, four-hole phases
7. ``OI_WAVELENGTH``: filter information

For more information on the format and contents of OIFITS files, see the `OIFITS2 standard <https://doi.org/10.1051/0004-6361/201526405>`_.

The _ami-oi.fits file contains tables of observables averaged over all integrations of the input file. The error is taken to be the standard error of the mean, where the variance is the covariance between amplitudes and phases (e.g. fringe amplitudes and fringe phases, closure phases and triple-product amplitudes, closure amplitudes and four-hole phases).
The _amimulti-oi.fits file contains observables for each integration, and does not contain error estimates. The
structure is the same as the _ami-oi.fits file, but the following data columns are 2D, with the second dimension being
the number of integrations: "PISTONS", "PIST_ERR", "VISAMP", "VISAMPERR", "VISPHI", "VISPHIERR", "VIS2DATA", "VIS2ERR", "T3AMP", "T3AMPERR", "T3PHI", "T3PHIERR", "Q4AMP","Q4AMPERR", "Q4PHI", "Q4PHI_ERR".

LG model parameters
^^^^^^^^^^^^^^^^^^^
:Data model: `~stdatamodels.jwst.datamodels.AmiLgFitModel`
:File suffix: _amilg.fits

The _amilg.fits output file contains the cropped and cleaned data, model, and residuals (data - model) as well as
the parameters of the best-fit LG model. It contains the following extensions:

1. ``CTRD``: a 3D image of the centered, cropped data
2. ``N_CTRD``: a 3D image CTRD normalized by data peak
3. ``FIT``: a 3D image of the best-fit model
4. ``N_FIT``: a 3D image of FIT normalized by data peak
5. ``RESID``: a 3D image of the fit residuals
6. ``N_RESID``: a 3D image of RESID normalized by data peak
7. ``SOLNS``: table of fringe coefficients
