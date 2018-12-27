.. _calwebb_ami3:

calwebb_ami3: Stage 3 Aperture Masking Interferometry (AMI) Processing
======================================================================

:Config: calwebb_ami3.cfg
:Class: `~jwst.pipeline.Ami3Pipeline`

The stage 3 AMI pipeline is to be applied to
associations of calibrated NIRISS AMI exposures and is used to compute fringe
parameters and, optionally, correct science target fringe parameters using
observations of reference PSF targets.
The steps applied by the ``calwebb_ami3`` pipeline are shown below.

+------------------------------------------+
| calwebb_ami3                             |
+==========================================+
| :ref:`ami_analyze <ami_analyze_step>`    |
+------------------------------------------+
| :ref:`ami_average <ami_average_step>`    |
+------------------------------------------+
| :ref:`ami_normalize <ami_normalize_step>`|
+------------------------------------------+

If the input to the pipeline is an ASN containing multiple target and reference
PSF exposures:

 - the :ref:`ami_analyze <ami_analyze_step>` step is applied to each input exposure
   independently, computing fringe parameters for each
 - the :ref:`ami_average <ami_average_step>` step computes the average of the
   :ref:`ami_analyze <ami_analyze_step>` results
   for all of the target exposures, and an average for all of the reference
   PSF results (if present)
 - the :ref:`ami_normalize <ami_normalize_step>` step corrects the average target results using
   the average reference PSF results (if present)

Arguments
---------
The ``calwebb_ami3`` pipeline has one optional argument:
::
 --save_averages  boolean  default=False

If set to ``True``, the results of the :ref:`ami_average <ami_average_step>` step will be saved.

Inputs
------

2D calibrated images
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _cal

The inputs to ``calwebb_ami3`` need to be in the form of an ASN file that lists
multiple science target exposures, and optionally reference PSF exposures as well.
The individual exposures must be in the form of calibrated ("_cal") products from
:ref:`calwebb_image2 <calwebb_image2>` processing.

Outputs
-------

Fringe parameter tables
^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.AmiLgModel`
:File suffix: _ami

For every input exposure, the fringe parameters and closure phases caculated
by the :ref:`ami_analyze <ami_analyze_step>` step are saved to an "_ami" product file, which
is a FITS table containing the fringe parameters and closure phases. Product names
use the input "_cal" exposure-based file name, with the association candidate ID
included and the product type changed to "_ami", e.g.
"jw93210001001_03101_00001_nis_a0003_ami.fits."

Averaged fringe parameters table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.AmiLgModel`
:File suffix: _amiavg or _psf-amiavg

If multiple target or reference PSF exposures are used as input and the
"--save_averages" parameter is set to ``True``, the :ref:`ami_average <ami_average_step>` step
will save averaged results for the target in an "_amiavg" product and for the
reference PSF in a "_psf-amiavg" product. The file name root will use the
source-based output product name given in the ASN file. These files are the
same FITS table format as the "_ami" products.

Normalized fringe parameters table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.AmiLgModel`
:File suffix: _aminorm

If reference PSF exposures are included in the input ASN, the averaged AMI results
for the target will be normalized by the averaged AMI results for the reference PSF,
via the :ref:`ami_normalize <ami_normalize_step>` step, and will be saved to an "_aminorm"
product file. This file has the same FITS table format as the "_ami" products.
The file name root uses the source-based output product name given in the ASN file,
e.g. "jw93210-a0003_t001_niriss_f480m-nrm_aminorm.fits."
