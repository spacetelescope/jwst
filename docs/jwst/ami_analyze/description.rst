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
The ``ami_analyze`` step has several optional arguments. In most cases the 
default arguments will be suitable but more advanced users may wish to test
other options:

:--oversample: The oversampling factor to be used in the model fit (default=3).

:--rotation: Initial guess for the rotation of the PSF in the input image, in
             units of degrees (default=0.0).

:--psf_offset: List of PSF offset values to use when creating the model array
               (default='0.0 0.0').

:--rotation_search: List of start, stop, and step values that define the list of
                    rotation search values. The default setting of '-3 3 1'
                    results in search values of [-3, -2, -1, 0, 1, 2, 3].

:--bandpass: ASDF file containing suitable array to override filter/source 
             (default=None)

:--usebp: If True, exclude pixels marked DO_NOT_USE from fringe fitting 
          (default=True)

:--firstfew: If not None, process only the first few integrations (default=None)

:--chooseholes: If not None, fit only certain fringes e.g. ['B4','B5','B6','C2']
                (default=None)

:--affine2d: ASDF file containing user-defined affine parameters (default='commissioning')

:--run_bpfix: Run Fourier bad pixel fix on cropped data (default=True)


Note that the `affine2d` default argument is a special case; 'commissioning' is currently the only string other than an ASDF filename that is accepted. If `None` is passed, it will perform a rotation search (least-squares fit to a PSF model) and use that for the affine transform.


Creating ASDF files
^^^^^^^^^^^^^^^^^^^
The optional arguments `bandpass` and `affine2d` must be written to `ASDF <https://asdf-standard.readthedocs.io/>`_ 
files to be used by the step. The step expects the contents to be stored with particular keys but the format is not currently
enforced by a schema; incorrect ASDF file contents will cause the step to revert back to the defaults for each argument.

Examples of how to create ASDF files containing the properly formatted information for each of the arguments follows.

.. code-block:: python

   # Create a F480M filter + Vega bandpass ASDF file

   import asdf
   from jwst.ami import utils
   from stdatamodels.jwst import datamodels
   from synphot import SourceSpectrum

   # F480M throughput reference file from JWST CRDS
   throughput_file = 'jwst_niriss_throughput_0012.fits'
   nspecbin=19
   throughput_model = datamodels.open(throughput_file)

   filt_spec = utils.get_filt_spec(throughput_model)
   src_spec = SourceSpectrum.from_vega()  
   bandpass = utils.combine_src_filt(filt_spec,
                                    src_spec,
                                    trim=0.01,
                                    nlambda=nspecbin)

   # This bandpass has shape (19, 2); each row is [throughput, wavelength]
   asdf_name = 'bandpass_f480m_vega.asdf'
   tree = {"bandpass": bandpass}
   with open(asdf_name, 'wb') as fh:
        af = asdf.AsdfFile(tree)
        af.write_to(fh)
   af.close()
   throughput_model.close()


.. code-block:: python

   # Create an affine transform ASDF file to use for the model

   import asdf
   tree = {
         'mx': 1., # dimensionless x-magnification
         'my': 1., # dimensionless y-magnification
         'sx': 0., # dimensionless x shear
         'sy': 0., # dimensionless y shear
         'xo': 0., # x-offset in pupil space
         'yo': 0., # y-offset in pupil space
         'rotradccw': None 
         }

   affineasdf = 'affine.asdf'

   with open(affineasdf, 'wb') as fh:
        af = asdf.AsdfFile(tree)
        af.write_to(fh)
   af.close()



Inputs
------

3D calibrated image
^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.DataModel`
:File suffix: _calints

The ``ami_analyze`` step takes a single calibrated image cube as input, which should be
the "_calints" product resulting from :ref:`calwebb_image2 <calwebb_image2>` processing.
Multiple exposures can be processed via use of an ASN file that is used as input
to the :ref:`calwebb_ami3 <calwebb_ami3>` pipeline. **Note:** The ``ami_analyze`` step will also
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
:Data model: `~jwst.datamodels.AmiOIModel`
:File suffix: _ami-oi.fits, _amimulti-oi.fits

The inteferometric observables are saved as OIFITS files, a registered FITS format
for optical interferometry, containing the following list of extensions:

1)  ``OI_ARRAY``: AMI subaperture information  
2)  ``OI_TARGET``: target properties  
3)  ``OI_T3``: extracted closure amplitudes, triple-product phases 
4)  ``OI_VIS``: extracted visibility (fringe) amplitudes, phases
5)  ``OI_VIS2``: squared visibility (fringe) amplitudes
6)  ``OI_WAVELENGTH``: filter information

For more information on the format and contents of OIFITS files, see the `OIFITS2 standard <https://doi.org/10.1051/0004-6361/201526405>`_.

The _ami-oi.fits file contains tables of observables averaged over all integrations of the input file. The error is taken to be the standard error of the mean, where the variance is the covariance between amplitudes and phases (e.g. fringe amplitudes and fringe phases, closure phases and triple-product amplitudes).
The _amimulti-oi.fits file contains observables for each integration, and does not contain error estimates. The
structure is the same as the _ami-oi.fits file, but the following data columns are 2D, with the second dimension being 
the number of integrations: "PISTONS", "PIST_ERR", "VISAMP", "VISAMPERR", "VISPHI", "VISPHIERR", "VIS2DATA", "VIS2ERR", "T3AMP", "T3AMPERR", "T3PHI", "T3PHIERR".

LG model parameters
^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiLgFitModel`
:File suffix: _amilg.fits

The _amilg.fits output file contains the cropped and cleaned data, model, and residuals (data - model) as well as 
the parameters of the best-fit LG model. It contains the following extensions:

1) ``CTRD``: a 3D image of the centered, cropped data
2) ``N_CTRD``: a 3D image CTRD normalized by data peak
3) ``FIT``: a 3D image of the best-fit model
4) ``N_FIT``: a 3D image of FIT normalized by data peak
5) ``RESID``: a 3D image of the fit residuals
6) ``N_RESID``: a 3D image of RESID normalized by data peak
7) ``SOLNS``: table of fringe coefficients

Reference Files
---------------
The ``ami_analyze`` step uses a THROUGHPUT reference file and NRM reference file.

.. include:: ../references_general/throughput_reffile.inc
.. include:: ../references_general/nrm_reffile.inc
