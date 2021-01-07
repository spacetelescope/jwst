.. _calwebb_ami3:

calwebb_ami3: Stage 3 Aperture Masking Interferometry (AMI) Processing
======================================================================

:Config: calwebb_ami3.cfg
:Class: `~jwst.pipeline.Ami3Pipeline`

The stage 3 AMI pipeline is applied to associations of calibrated NIRISS AMI exposures.
It computes fringe parameters for individual exposures, averages the fringe results from
multiple exposures, and, optionally, corrects science target fringe parameters using the
fringe results from reference PSF targets.
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

When given an association file as input, which lists multiple science target and reference PSF
exposures, the pipeline will:

 - apply the :ref:`ami_analyze <ami_analyze_step>` step to each input exposure
   independently, computing fringe parameters for each
 - apply the :ref:`ami_average <ami_average_step>` step to compute the average of the
   :ref:`ami_analyze <ami_analyze_step>` results for all of the science target exposures,
   and the average for all of the reference PSF results (if present)
 - apply the :ref:`ami_normalize <ami_normalize_step>` step to correct the average science
   target results using the average reference PSF results (if present)

If no reference PSF target exposures are present in the input ASN file, the ``ami_normalize``
step is skipped.

Arguments
---------
The ``calwebb_ami3`` pipeline has one optional argument::

  --save_averages  boolean  default=False

If set to ``True``, the results of the :ref:`ami_average <ami_average_step>` step will be saved
to a file. If not, the results of the :ref:`ami_average <ami_average_step>` step are passed
along in memory to the :ref:`ami_normalize <ami_normalize_step>` step.

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

An example ASN file containing 2 science target and 2 reference PSF target exposures is
shown below. Only 1 product is defined, corresponding to the science target, with
members consisting of exposures for both the science target and the reference PSF target,
as indicated by the "exptype" values for each.
::

 {"asn_type": "ami3",
  "asn_rule": "discover_Asn_AMI",
  "program": "10005",
  "asn_id": "a3001",
  "target": "t001",
  "asn_pool": "jw10005_001_01_pool",
  "products": [
      {"name": "jw10005-a3001_t001_niriss_f277w-nrm",
       "members": [
           {"expname": "jw10005007001_02102_00001_nis_cal.fits",
            "exptype": "psf"
           },
           {"expname": "jw10005027001_02102_00001_nis_cal.fits",
            "exptype": "psf"
           },
           {"expname": "jw10005004001_02102_00001_nis_cal.fits",
            "exptype": "science"
           },
           {"expname": "jw10005001001_02102_00001_nis_cal.fits",
            "exptype": "science"
           }
       ]
      }
  ]
 }

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
