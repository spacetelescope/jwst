.. _calwebb_tso3:

calwebb_tso3: Stage 3 Time-Series Observation(TSO) Processing
=============================================================

:Class: `jwst.pipeline.Tso3Pipeline`
:Alias: calwebb_tso3

The stage 3 TSO pipeline is to be applied to associations of calibrated TSO exposures
(e.g. NIRCam TS imaging, NIRCam TS grism, NIRISS SOSS, NIRSpec BrightObj, MIRI LRS Slitless)
and is used to produce calibrated time-series photometry or spectra of the source object.

The steps applied by the ``calwebb_tso3`` pipeline for Imaging and Spectroscopic TSO
exposures are shown below:

.. |check| unicode:: U+2713

.. checkmark

+---------------------------------------------------+---------+--------------+
| calwebb_tso3                                      | Imaging | Spectroscopy |
+===================================================+=========+==============+
| :ref:`outlier_detection <outlier_detection_step>` | |check| | |check|      |
+---------------------------------------------------+---------+--------------+
| :ref:`tso_photometry <tso_photometry_step>`       | |check| |              |
+---------------------------------------------------+---------+--------------+
| :ref:`extract_1d <extract_1d_step>`               |         | |check|      |
+---------------------------------------------------+---------+--------------+
| :ref:`white_light <white_light_step>`             |         | |check|      |
+---------------------------------------------------+---------+--------------+

The logic that decides whether to apply the imaging or spectroscopy steps is based
on the EXP_TYPE and TSOVISIT keyword values of the input data. Imaging steps are
applied if either of the following is true:

 - EXP_TYPE = 'NRC_TSIMAGE'
 - EXP_TYPE = 'MIR_IMAGE' and TSOVISIT = True

The spectroscopy steps will be applied in all other cases.

Inputs
------

3D calibrated images
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _calints

The input to ``calwebb_tso3`` is in the form of an ASN file that lists multiple
exposures or exposure segments of a science target. The individual inputs should be in
the form of 3D calibrated ("_calints") products from either :ref:`calwebb_image2 <calwebb_image2>`
or :ref:`calwebb_spec2 <calwebb_spec2>` processing. These products contain 3D stacks of
per-integration images. Each pipeline step will loop over all of the integrations in each
input.

Many TSO exposures may contain a sufficiently large number of integrations (NINTS) so as to make
their individual exposure products too large (in terms of file size on disk) to be able to handle
conveniently. In these cases, the uncalibrated raw data for a given exposure are split into
multiple "segmented" products, each of which is identified with a segment number
(see :ref:`segmented products <segmented_files>`). The ``calwebb_tso3`` input ASN file includes
all "_calints" exposure segments. The :ref:`outlier_detection <outlier_detection_step>` step will
process a single segment at a time, creating one output "_crfints" product per segment. The
remaining ``calwebb_tso3`` steps, will process each segment and concatenate the results into a
single output product, containing the results for all exposures and segments listed in the ASN.

Outputs
-------

3D CR-flagged images
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _crfints

If the :ref:`outlier_detection <outlier_detection_step>` step is applied, a new version
of each input calibrated product is created, which contains a DQ array
that has been updated to flag pixels detected as outliers. This updated
product is known as a CR-flagged product and is saved as a "_crfints" product type.

Imaging photometry
^^^^^^^^^^^^^^^^^^
:Data model: N/A
:File suffix: _phot

For imaging TS observations, the :ref:`tso_photometry <tso_photometry_step>` step produces
a source catalog containing photometry results from all of the "_crfints" products, organized
as a function of integration time stamps.
This file is saved in ASCII "ecsv" format, with a product type of "_phot." The file naming is
source-based, using the output product name specified in the ASN file, e.g.
"jw93065-a3001_t1_nircam_f150w-wlp8_phot.ecsv."

1D extracted spectral data
^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.MultiSpecModel`
:File suffix: _x1dints

For spectroscopic TS observations, the :ref:`extract_1d <extract_1d_step>` step is applied to
all "_crfints" products, to create a single "_x1dints" product that contains 1D extracted
spectral data for all integrations contained in the input exposures. The file name is
source-based, using the output product name specified in the ASN file, e.g.
"jw87600-a3001_t001_niriss_clear-gr700xd_x1dints.fits."

Spectroscopic white-light photometry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: N/A
:File suffix: _whtlt

For spectroscopic TS observations, the :ref:`white_light <white_light_step>` step is applied
to all of the 1D extracted spectral data in the "_x1dints" product, to produce an ASCII catalog
in ``ecsv`` format containing the wavelength-integrated white-light photometry of the source.
The catalog lists the integrated white-light flux as a function of time, based on the
integration time stamps. The file name is source-based, using the output product name specified
in the ASN file, e.g.
"jw87600-a3001_t001_niriss_clear-gr700xd_whtlt.ecsv."
