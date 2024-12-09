.. _calwebb_ami3:

calwebb_ami3: Stage 3 Aperture Masking Interferometry (AMI) Processing
======================================================================

:Class: `jwst.pipeline.Ami3Pipeline`
:Alias: calwebb_ami3

The stage 3 AMI pipeline is applied to associations of calibrated NIRISS AMI exposures.
It computes fringe parameters for individual exposures, and, optionally, corrects science target fringe parameters using the
fringe results from reference PSF targets.
The steps applied by the ``calwebb_ami3`` pipeline are shown below.

+------------------------------------------+
| calwebb_ami3                             |
+==========================================+
| :ref:`ami_analyze <ami_analyze_step>`    |
+------------------------------------------+
| :ref:`ami_normalize <ami_normalize_step>`|
+------------------------------------------+

When given an association file as input, which lists multiple science target and reference PSF
exposures, the pipeline will:

#. apply the :ref:`ami_analyze <ami_analyze_step>` step to each input exposure
   independently, computing fringe parameters for each
#. apply the :ref:`ami_normalize <ami_normalize_step>` step to correct the science
   target results using the reference PSF results (if present)

If no reference PSF target exposures are present in the input ASN file, the ``ami_normalize``
step is skipped.

Arguments
---------
The ``calwebb_ami3`` pipeline does not currently use any optional arguments.

Inputs
------

3D calibrated images
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.DataModel`
:File suffix: _calints

The inputs to ``calwebb_ami3`` need to be in the form of an ASN file that lists
multiple science target exposures, and optionally reference PSF exposures as well.
The individual exposures must be in the form of 3D calibrated ("_calints") products from
:ref:`calwebb_image2 <calwebb_image2>` processing.

An example ASN file containing one science target and one reference PSF target exposure is
shown below. Only 1 product is defined, corresponding to the science target, with
members consisting of exposures for both the science target and the reference PSF target,
as indicated by the "exptype" values for each.
::

 {"asn_type": "ami3",
  "asn_rule": DMS_Level3_Base",
  "program": "10005",
  "asn_id": "a3001",
  "target": "t001",
  "asn_pool": "jw10005_001_01_pool",
  "products": [
      {"name": "jw10005-a3001_t001_niriss_f277w-nrm",
       "members": [
           {"expname": "jw10005007001_02102_00001_nis_calints.fits",
            "exptype": "psf"
           },
           {"expname": "jw10005004001_02102_00001_nis_calints.fits",
            "exptype": "science"
           }
       ]
      }
  ]
 }

Outputs
-------

Interferometric observables
^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiOIModel`
:File suffix: _ami-oi.fits

For every input exposure, the interferometric observables calculated
by the :ref:`ami_analyze <ami_analyze_step>` step are saved to an "_ami-oi.fits" product file,
which is a FITS table of averaged observables over all integrations of the input file.
Product names use the input "_calints" exposure-based file name, with the association candidate ID
included and the product type changed to "_ami-oi.fits", e.g.
"jw93210001001_03101_00001_nis_a0003_ami-oi.fits."

Normalized interferometric observables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiOIModel`
:File suffix: _aminorm-oi.fits

If reference PSF exposures are included in the input ASN, the AMI results
for the target will be normalized by the AMI results for the reference PSF,
via the :ref:`ami_normalize <ami_normalize_step>` step, and will be saved to an "_aminorm-oi.fits"
product file. This file has the same FITS table format as the "_ami-oi.fits" products.
The file name root uses the source-based output product name given in the ASN file,
e.g. "jw93210-a0003_t001_niriss_f480m-nrm_aminorm-oi.fits."

.. note:: 
   
   Users may wish to run the :ref:`ami_analyze step <ami_analyze_step>` separately for access to flexible input arguments and to save additional diagnostic output products. See the step documentation for more details.

