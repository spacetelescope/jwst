Tasks in the Package
====================
The Aperture Masking Interferometry (AMI) package currently consists
of three tasks:

1) ``ami_analyze``: apply the LG algorithm to a NIRISS AMI exposure
2) ``ami_average``: average the results of LG processing for multiple exposures
3) ``ami_normalize``: normalize the LG results for a science target using LG
   results from a reference target

The three tasks can be applied to an association of AMI exposures using the
pipeline module ``calwebb_ami3``.

CALWEBB_AMI3 Pipeline
=====================

Overview
--------
The ``calwebb_ami3`` pipeline module can be used to apply all 3 steps of AMI
processing to an association (ASN) of AMI exposures. The processing flow through
the pipeline is as follows:

1) Apply the ``ami_analyze`` step to all products listed in the input
   association table. Output files will have a product type suffix of ``ami``.
   There will be one ``ami`` product per input exposure.

2) Apply the ``ami_average`` step to combine the above results for both
   science target and reference target exposures, if both types exist in the
   ASN table. If the optional parameter ``save_averages`` is set to true
   (see below), the results will be saved to output files with a product type
   suffix of ``amiavg``.
   There will be one ``amiavg`` product for the science target and one for
   the reference target.

3) If reference target results exist, apply the ``ami_normalize`` step to the
   averaged science target result, using the averaged reference target result
   to do the normalization. The output file will have a product type suffix of
   ``aminorm``.

Input
-----
The only input to the ``calwebb_ami3`` pipeline is the name of a json-formatted
association file. There is one optional parameter ``save_averages``. If set to
true, the results of the ``ami_average`` step will be saved to files.
It is assumed that the
ASN file will define a single output product for the science target result,
containing a list of input member file names, for both science target and
reference target exposures. An example ASN file is shown below.

::

 {"asn_rule": "NIRISS_AMI", "targname": "NGC-3603", "asn_pool": "jw00017_001_01_pool", "program": "00017",
 "products": [
     {"prodtype": "ami", "name": "jw87003-c1001_t001_niriss_f277w-nrm",
      "members": [
         {"exptype": "science", "expname": "targ14_cal.fits"},
         {"exptype": "science", "expname": "targ15_cal.fits"},
         {"exptype": "science", "expname": "targ16_cal.fits"},
         {"exptype": "psf", "expname": "ref1_cal.fits"},
         {"exptype": "psf", "expname": "ref2_cal.fits"},
         {"exptype": "psf", "expname": "ref3_cal.fits"}]}],
 "asn_type": "ami",
 "asn_id": "c1001"}

Note that the ``exptype`` attribute value for each input member is used to
indicate which files contain science target images and which contain reference
PSF images.

.. _ami_analyze_step:

AMI_Analyze
===========

Overview
--------
The ``ami_analyze`` step applies the Lacour-Greenbaum (LG) image plane
modeling algorithm to a NIRISS AMI image.
The routine computes a number of parameters, including a model fit (and
residuals) to the image, fringe amplitudes and phases, and closure phases
and amplitudes.

The JWST AMI observing template allows for exposures to be obtained using
either full-frame (SUBARRAY="FULL") or subarray (SUBARRAY="SUB80") readouts.
When processing a full-frame exposure, the ``ami_analyze`` step extracts
(on the fly) a region from the image corresponding to the size and location of
the SUB80 subarray, in order to keep the processing time to a reasonable level.

Inputs
------
The ``ami_analyze`` step takes a single input image, in the form of a simple 2D
ImageModel. There are two optional parameters:

1) ``oversample``: specifies the oversampling factor to be used in the model fit
   (default value = 3)
2) ``rotation``: specifies an initial guess, in degrees, for the rotation of the
   PSF in the input image (default value = 0.0)

Output
------
The ``ami_analyze`` step produces a single output file, which contains the
following list of extensions:

1) ``FIT``: a 2-D image of the fitted model
2) ``RESID``: a 2-D image of the fit residuals
3) ``CLOSURE_AMP``: table of closure amplitudes
4) ``CLOSURE_PHA``: table of closure phases
5) ``FRINGE_AMP``: table of fringe amplitudes
6) ``FRINGE_PHA``: table of fringe phases
7) ``PUPIL_PHA``: table of pupil phases
8) ``SOLNS``: table of fringe coefficients

.. _ami_average_step:

AMI_Average
===========

Overview
--------
The ``ami_average`` step averages the results of LG processing from the
``ami_analyze`` step for multiple exposures of a given target. It averages
all 8 components of the ``ami_analyze`` output files for all input exposures.

Inputs
------
The only input to the ``ami_average`` step is a list of input files to be
processed. These will presumably be output files from the ``ami_analyze`` step.
The step has no other required or optional parameters, nor does it use any
reference files.

Output
------
The step produces a single output file, having the same format as the input
files, where the data for the 8 file components
are the average of each component from the list of input files.

.. _ami_normalize_step:

AMI_Normalize
=============

Overview
--------
The ``ami_normalize`` step provides normalization of LG processing results for
a science target using LG results of a reference target. The algorithm
subtracts the reference target closure phases from the science target closure
phases and divides the science target fringe amplitudes by the reference target
fringe amplitudes.

Inputs
------
The ``ami_normalize`` step takes two input files: the first is the LG
processed results for a science target and the second is the LG processed
results for the reference target. There are no optional parameters and no
reference files are used.

Output
------
The output is a new LG product for the science target in which the closure
phases and fringe amplitudes have been normalized using the reference target
closure phases and fringe amplitudes. The remaining components of the science
target data model are left unchanged.
