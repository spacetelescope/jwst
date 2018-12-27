.. _calwebb_coron3:

calwebb_coron3: Stage 3 Coronagraphic Processing
================================================

:Config: calwebb_coron3.cfg
:Class: `~jwst.pipeline.Coron3Pipeline`

The stage 3 coronagraphic pipeline is to be applied to
associations of calibrated NIRCam coronagraphic and MIRI Lyot and 4QPM
exposures, and is used to produce PSF-subtracted, resampled, combined images
of the source object.

The steps applied by the ``calwebb_coron3`` pipeline are shown in the table
below.

+---------------------------------------------------+
| calwebb_coron3                                    |
+===================================================+
| :ref:`stack_refs <stack_refs_step>`               |
+---------------------------------------------------+
| :ref:`align_refs <align_refs_step>`               |
+---------------------------------------------------+
| :ref:`klip <klip_step>`                           |
+---------------------------------------------------+
| :ref:`outlier_detection <outlier_detection_step>` |
+---------------------------------------------------+
| :ref:`resample <resample_step>`                   |
+---------------------------------------------------+

Arguments
---------

The ``calwebb_coron3`` pipeline does not have any optional arguments.

Inputs
------

3D calibrated image stacks
^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _calints

The input to ``calwebb_coron3`` must be in the form of an ASN file that lists
one or more exposures of a science target and a reference PSF target.
The individual target and reference PSF exposures should be in the form of 3D
calibrated ("_calints") products from :ref:`calwebb_image2 <calwebb_image2>`
processing. Each pipeline step will loop over the 3D stack of per-integration images
contained in each exposure.

Outputs
-------

3D stacked PSF images
^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _psfstack

The data from each input PSF reference exposure are concatenated into a single
combined 3D stack by the :ref:`stack_refs <stack_refs_step>` step, for use by subsequent
steps. The stacked PSF data get written to disk in the form of a "_psfstack" product.
The output file name is source-based, using the product name specified in the
ASN file, e.g. "jw86073-a3001_t001_nircam_f140m-maskbar_psfstack.fits."

4D aligned PSF images
^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.QuadModel`
:File suffix: _psfalign

For each science target exposure, all of the reference PSF images in the
"_psfstack" product are aligned to each science target integration and saved to
a 4D "_psfalign" product by the :ref:`align_refs <align_refs_step>` step. The output file
name is exposure-based, with the addition of the associated candidate ID, e.g.
"jw8607342001_02102_00001_nrcb3_a3001_psfalign.fits."

3D PSF-subtracted images
^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _psfsub

For each science target exposure, the :ref:`klip <klip_step>` step applies PSF fitting and
subtraction for each integration, resulting in a 3D stack of PSF-subtracted
images. The data for each science target exposure are saved to a "_psfsub"
product, using exposure-based file names, e.g.
"jw8607342001_02102_00001_nrcb3_a3001_psfsub.fits."

CR-flagged images
^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _crfints

If the :ref:`outlier_detection <outlier_detection_step>` step is applied, a new version of
each "_psfsub" product is created, in which the DQ array is updated to flag pixels detected
as outliers. These files use the "_crfints" (CR-Flagged per integration)
product type suffix and include the association candidate ID, e.g.
"jw8607342001_02102_00001_nrcb3_a3001_crfints.fits."

2D resampled image
^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.DrizProductModel`
:File suffix: _i2d

The :ref:`resample <resample_step>` step is applied to the CR-flagged products to create a
single resampled and combined product for the science target. The file name is
source-based, using the product name specified in the ASN file, e.g.
"jw86073-a3001_t001_nircam_f140m-maskbar_i2d.fits."
