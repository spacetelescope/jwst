Description
-----------

:Class: `jwst.ami.AmiAnalyzeStep`
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

Arguments
---------
The ``ami_analyze`` step has four optional arguments:

:--oversample: The oversampling factor to be used in the model fit (default=3).

:--rotation: Initial guess for the rotation of the PSF in the input image, in
             units of degrees (default=0.0).

:--psf_offset: List of PSF offset values to use when creating the model array
               (default='0.0 0.0').

:--rotation_search: List of start, stop, and step values that define the list of
                    rotation search values. The default setting of '-3 3 1'
                    results in search values of [-3, -2, -1, 0, 1, 2, 3].
                    

Inputs
------

2D calibrated image
^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _cal

The ``ami_analyze`` step takes a single calibrated image as input, which should be
the "_cal" product resulting from :ref:`calwebb_image2 <calwebb_image2>` processing.
Multiple exposures can be processed via use of an ASN file that is used as input
to the :ref:`calwebb_ami3 <calwebb_ami3>` pipeline. The ``ami_analyze`` step itself does
not accept an ASN as input.

Outputs
-------

LG model parameters
^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiLgModel`
:File suffix: _ami

The ``ami_analyze`` step produces a single output file, containing the
following list of extensions:

1) ``FIT``: a 2D image of the fitted model
2) ``RESID``: a 2D image of the fit residuals
3) ``CLOSURE_AMP``: table of closure amplitudes
4) ``CLOSURE_PHA``: table of closure phases
5) ``FRINGE_AMP``: table of fringe amplitudes
6) ``FRINGE_PHA``: table of fringe phases
7) ``PUPIL_PHA``: table of pupil phases
8) ``SOLNS``: table of fringe coefficients

The output file name syntax is exposure-based, using the input file name as the root, with
the addition of the association candidate ID and the "_ami" product type suffix, e.g.
"jw87600027001_02101_00002_nis_a3001_ami.fits."

Reference Files
---------------
The ``ami_analyze`` step uses a THROUGHPUT reference file.

.. include:: ../references_general/throughput_reffile.inc
