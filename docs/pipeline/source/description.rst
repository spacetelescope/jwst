Pipeline Classes and Configuration Files
========================================
The actual pipelines that call individual correction steps in various
orders are defined as python classes within python code modules. The pipelines
can be executed by referencing their class name or through the use of a
configuration (.cfg) file that in turn references the class. The table below
shows the pipeline classes that are currently available, the
corresponding pre-defined configurations that make use of those classes, and
the instrument modes to which they can be applied.

+----------------------+------------------------+------------------------------------------+
| Class Name           | Configuration File     | Used For                                 |
+======================+========================+==========================================+
| SloperPipeline       | calwebb_sloper.cfg     | Level-2a processing: all modes           |
+----------------------+------------------------+------------------------------------------+
| DarkPipeline         | calwebb_dark.cfg       | Level-2a processing: darks               |
+----------------------+------------------------+------------------------------------------+
| Image2Pipeline       | calwebb_image2.cfg     | Level-2b processing: imaging modes       |
+----------------------+------------------------+------------------------------------------+
| Spec2Pipeline        | calwebb_spec2.cfg      | Level-2b processing: spectroscopy modes  |
+----------------------+------------------------+------------------------------------------+
| Image3Pipeline       | calwebb_image3.cfg     | Level-3 processing: imaging modes        |
+----------------------+------------------------+------------------------------------------+
| Ami3Pipeline         | calwebb_ami3.cfg       | Level-3 processing: NIRISS AMI mode      |
+----------------------+------------------------+------------------------------------------+

Input Files, Output Files and Data Models
=========================================
An important concept used throughout the JWST pipeline is the :ref:`Data
Model`. Nearly all data used by any of the pipeline code is
encapsulated in a data model. Most input is read into a data model and
all output is produced by a data model. When possible, this document
will indicate the data model associated with a file type, usually as a
parenthetical link to the data model in question. For some steps, the
output file may represent different data models depending on the input
to those steps. As a result, the data models listed here will not be
an exhaustive list.

Level-2a Pipeline Step Flow (calwebb_sloper)
=============================================
Level-2a processing applies basic detector-level corrections to all exposure
types (imaging, spectroscopic, coronagraphic, etc.). It is applied to one
exposure at a time. The pipeline module for level-2a processing is
``calwebb_sloper`` (the equivalent pipeline class is ``SloperPipeline``). It is
often referred to as ``ramps-to-slopes`` processing, because the input raw data
are in the form of one or more ramps (integrations) containing accumulating
counts from the non-destructive detector readouts and the output is a corrected
countrate (slope) image. The list of steps applied by the Build 7 calwebb_sloper
pipeline is as follows.

==============  ==============
calwebb_sloper  calwebb_sloper
(All Near-IR)   (MIRI)
==============  ==============
dq_init         dq_init
saturation      saturation
ipc             ipc       
superbias       linearity 
refpix          rscd
linearity       lastframe    
dark_current    dark_current 
\               refpix
jump            jump
ramp_fit        ramp_fit
==============  ==============

Inputs
------

* Raw 4D product: The input to ``calwebb_sloper`` is a single raw exposure file,
  e.g. ``jw80600012001_02101_00003_mirimage_uncal.fits``, which contains the
  original raw data from all of the detector readouts in the exposure
  (ncols x nrows x ngroups x nintegrations).

Outputs
-------

* 2D Countrate product: All types of inputs result in a 2D countrate product,
  resulting from averaging over all of the integrations within the exposure.
  The output file will be of type ``_rate``, e.g.
  ``jw80600012001_02101_00003_mirimage_rate.fits``.

* 3D Countrate product: If the input exposure contains more than one integration
  (NINTS>1), a 3D countrate product is created that contains the individual
  results of each integration. The 2D countrate images for each integration are
  stacked along the 3rd axis of the data cubes (ncols x nrows x nints). This
  output file will be of type ``_rateints``.

Arguments
---------
The ``calwebb_sloper`` pipeline has one optional argument:

* ``save_calibrated_ramp``

which is a boolean argument with a default value of ``False``. If the user sets
it to ``True``, the pipeline will save intermediate data to a file as it
exists at the end of the ``jump`` step (just before ramp fitting). The data at
this stage of the pipeline are still in the form of the original 4D ramps
(ncols x nrows x ngroups x nints) and have had all of the detector-level
correction steps applied to it, including the detection and flagging of
Cosmic-Ray hits within each ramp (integration). If created, the name of the
intermediate file will be constructed from the root name of the input file, with
the new product type suffix ``_ramp`` appended
(e.g. ``jw80600012001_02101_00003_mirimage_ramp.fits``).

Dark Pipeline Step Flow (calwebb_dark)
======================================
The Level-2a dark (``calwebb_dark``) processing pipeline is intended for use
with dark exposures. It applies all of the same detector-level correction steps
as the ``calwebb_sloper`` pipeline, but stops just before the application of the
``dark_current`` step.

Inputs
------

* Raw 4D Dark product: The input to ``calwebb_dark`` is a single raw dark
  exposure.

Outputs
-------

* 4D Corrected product: The output is a 4D (ncols x nrows x ngroups x nints)
  product that has had all corrections up to, but not including, the
  ``dark_current`` step, with a product file type of ``_dark``.

Arguments
---------
The ``calwebb_dark`` pipeline does not have any optional arguments.

Level-2b Imaging Pipeline Step Flow (calwebb_image2)
====================================================
Level-2b imaging (``calwebb_image2``) processing applies additonal corrections
that result in a fully calibrated individual exposure. The list of correction
steps applied by the calwebb_image2 level-2b imaging pipeline is as follows.

+----------------+
| calwebb_image2 |
+================+
| assign_wcs     |
+----------------+
| flat_field     |
+----------------+
| photom         |
+----------------+

Inputs
------

* 2D or 3D Countrate product: The input to the ``calwebb_image2`` pipeline is
  a single level-2a exposure, in the form of either a ``_rate`` or ``_rateints``
  file. If the latter (data on a per-integration basis), the steps in the
  pipeline are applied individually to each integration, where appropriate.

Outputs
-------

* 2D or 3D Calibrated product: The output is a single calibrated exposure, using
  the product type suffix ``_cal`` or ``_calints``, depending on the type of
  input (e.g. ``jw80600012001_02101_00003_mirimage_cal.fits``).

Arguments
---------
The ``calwebb_image2`` pipeline does not have any optional arguments.

Level-2b Spectroscopic Pipeline Step Flow (calwebb_spec2)
==========================================================
Level-2b spectroscopic (``calwebb_spec2``) processing applies additional
corrections to level-2a products that result in fully calibrated individual
exposures.
The list of correction steps is shown below. Some steps are only applied to
certain instruments or instrument modes, as noted in the table.

Note that level-2b processing for NIRCam and NIRISS Wide-Field Slitless (grism)
Spectroscopy modes is not yet implemented.

+------------------+----+-----+-----+----+----+-----+--------+
| Instrument Mode  |     NIRSpec    |     MIRI      | NIRISS |
+------------------+----+-----+-----+----+----+-----+--------+
| Step             | FS | MOS | IFU | FS | SL | MRS |  SOSS  |
+==================+====+=====+=====+====+====+=====+========+
| assign_wcs       | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| bkg_subtract     | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| imprint_subtract |    |  X  |  X  |    |    |     |        |
+------------------+----+-----+-----+----+----+-----+--------+
| extract_2d       | X  |  X  |     |    |    |     |        |
+------------------+----+-----+-----+----+----+-----+--------+
| flat_field       | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| srctype          | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| straylight       |    |     |     |    |    |  X  |        |
+------------------+----+-----+-----+----+----+-----+--------+
| fringe           |    |     |     |    |    |  X  |        |
+------------------+----+-----+-----+----+----+-----+--------+
| photom           | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| resample_spec    | X  |  X  |     |    |    |     |        |
+------------------+----+-----+-----+----+----+-----+--------+
| cube_build       |    |     |  X  |    |    |  X  |        |
+------------------+----+-----+-----+----+----+-----+--------+
| extract_1d       | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+

The ``resamp_spec`` step produces a resampled/rectified product for non-IFU
modes of some kinds of spectroscopic exposures. If the resample_spec step is not
applied to a given exposure, the extract_1d operation will be performed on the
original (unresampled) data.
The ``cube_build`` step produces a resampled/rectified cube for IFU exposures.

Inputs
------
The input to the ``calwebb_spec2`` pipeline can be either a single level-2a
(``_rate`` or ``_rateints``) exposure or an Association (ASN) file
listing multiple exposures. The background subtraction (``bkg_subtract``) and
imprint subtraction (``imprint_subtract``) steps can only be executed when
the pipeline is supplied with an association of exposures, because they rely
on multiple exposures to perform their tasks. The ASN file must not only list
the input exposures, but must also contain information that indicates their
relationships to one another.

The background subtraction step can be applied to an assocation containing
nodded exposures, such as for MIRI LRS fixed-slit, NIRSpec fixed-slit, and
NIRSpec MSA observations, or an association that contains dedicated exposures
of a background source. The step will accomplish background subtraction by
doing direct subtraction of nodded exposures from one another or by direct
subtraction of dedicated background expsoures from the science target exposures.

The imprint subtraction step, which is only applied to NIRSpec MSA and IFU
exposures, also requires the use of an ASN file, in order to specify which of
the inputs is to be used as the imprint exposure. The imprint exposure will be
subtracted from all other exposures in the association.

If a single level-2a product is used as input, the background subtraction
and imprint subtraction steps will be skipped and only the remaining regular
calibration steps will be applied to the input exposure.

Outputs
-------
Two or three different types of outputs are created by ``calwebb_spec2``.

* Calibrated 2D product: All types of inputs result in a fully-calibrated 2D
  product at the end of the ``photom`` step, which use the ``_cal`` or
  ``_calints`` product type suffix, depending on whether the input was a
  ``_rate`` or ``_rateints`` product, respectively.

* Resampled 2D product: If the input is an exposure type that gets
  resampled/rectified by the ``resample_spec`` step, the rectified 2D spectral
  product created by the ``resample_spec`` step is saved as a ``_s2d`` file.

* Resampled 3D product: If the data are NIRSpec IFU or MIRI MRS, the
  results of the ``cube_build`` step will be saved as a ``_s3d`` file.

* 1D Extracted Spectrum product: All types of inputs result in a 1D extracted
  spectral data product, which is saved as a ``_x1d`` file.

If the input to ``calwebb_spec2`` is an ASN file, these products are created
for each input exposure.

Arguments
---------
The ``calwebb_spec2`` pipeline has one optional argument:

* ``save_bsub``

which is a Boolean argument with a default value of ``False``. If the user sets
it to ``True``, the results of the background subtraction step (if applied) are
saved to an intermediate file of type ``_bsub`` or ``_bsubints``, as appropriate.

Level-3 Imaging Pipeline Step Flow (calwebb_image3)
===================================================
Level-3 processing for imaging observations is intended for combining the data
from multiple exposures (e.g. a dither or mosaic pattern) into a single
rectified (distortion corrected) product.
Before being combined, the exposures receive additional corrections for the
purpose of astrometric alignment, background matching, and outlier rejection.
The steps applied by the ``calwebb_image3`` pipeline are shown below.

+-------------------+
| calwebb_image3    |
+===================+
| tweakreg_catalog  |
+-------------------+
| tweakreg          |
+-------------------+
| skymatch          |
+-------------------+
| outlier_detection |
+-------------------+
| resample          |
+-------------------+
| source_catalog    |
+-------------------+

Inputs
------

* Associated 2D Calibrated products: The inputs to ``calwebb_image3`` will
  usually be in the form of an ASN file that lists multiple exposures to be
  processed and combined into a single output product. The individual exposures
  should be in the form of Level-2b (``_cal``) products from ``calwebb_image2``
  processing.

* Single 2D Calibrated product: It is also possible use a single ``_cal`` file
  as input to ``calwebb_image3``, in which case only the ``resample`` and
  ``source_catalog`` steps will be applied.

Outputs
-------

* Resampled 2D Image product (:ref:`DrizProductModel`): A resampled/rectified 2D image product of type
  ``_i2d`` is created containing the rectified single exposure or the rectified
  and combined association of exposures, which is the direct output of the
  ``resample`` step. This is the 

* Source catalog: A source catalog produced from the ``_i2d`` product is saved
  as an ASCII file in ``ecsv`` format, with a product type of ``_cat``.

* Level-2c products: If the ``outlier_detection`` step is applied, a new version
  of each input Level-2b exposure product is created, which contains a DQ array
  that has been updated to flag pixels detected as outliers. This updated
  product is known as a Level-2c product and the file is identified by appending
  the association candidate ID to the original input ``_cal`` file name, e.g.
  ``jw96090001001_03101_00001_nrca2_cal-o001.fits``.

Level-3 Aperture Masking Interferometry (AMI) Pipeline Step Flow (calwebb_ami3)
===============================================================================
The Level-3 AMI pipeline (``calwebb_ami3``) is intended to be applied to
associations of calibrated NIRISS AMI exposures and is used to compute fringe
parameters and correct science target fringe parameters using observations of
reference targets.
The steps applied by the ``calwebb_ami3`` pipeline are shown below.

+---------------+
| calwebb_ami3  |
+===============+
| ami_analyze   |
+---------------+
| ami_average   |
+---------------+
| ami_normalize |
+---------------+

Inputs
------

* Associated 2D Calibrated products: The inputs to ``calwebb_ami3`` are assumed
  to be in the form of an ASN file that lists multiple science target exposures,
  and optionally reference target exposures as well. The individual exposures
  should be in the form of Level-2b (``_cal``) products from ``calwebb_image2``
  processing.

Outputs
-------

* LG product (:ref:`AmiLgModel`): For every input exposure, the fringe
  parameters and closure phases caculated by the ``ami_analyze`` step
  are saved to an ``_lg`` product type file.

* Averaged LG product (:ref:`AmiLgModel`): The LG results averaged over all science or reference
  exposures, calculated by the ``ami_average`` step, are saved to an ``_lgavgt``
  (for the science target) or ``_lgavgr`` (for the reference target) file. Note
  that these output products are only created if the pipeline argument
  ``save_averages`` (see below) is set to ``True``. 

* Normalized LG product (:ref:`AmiLgModel`): If reference target exposures are included in the input
  ASN, the LG results for the science target will be normalized by the LG
  results for the reference target, via the ``ami_normalize`` step, and will be
  saved to an ``_lgnorm`` product file.

Arguments
---------
The ``calwebb_ami3`` pipeline has one optional argument:

* ``save_averages``

which is a Boolean parameter set to a default value of ``False``. If the user
sets this agument to ``True``, the results of the ``ami_average`` step will be
saved, as described above.
