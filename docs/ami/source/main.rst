Tasks in the Package
====================

The AMI pipeline package currently consists of three tasks:

1) ami_analyze: apply the LG algorithm to a NIRISS AMI exposure
2) ami_average: average the results of LG processing for multiple exposures
3) ami_normalize: normalize the LG results for a science target using LG
   results for a reference target

The three tasks can be applied to an association of AMI exposures using the
pipeline module 'calweb_ami3.'

AMI_Analyze
===========

Overview
--------
The 'ami_analyze' step applies the Lacour-Greenbaum (LG) image plane
modeling algorithm to a NIRISS Aperture Masking Interferometry (AMI) image.
The routine computes a number of parameters, including a model fit (and
residuals) to the image, fringe amplitudes and phases, and closure phases
and amplitudes.

Inputs
------
The ami_analyze step takes a single input image, in the form of a simple 2-D
ImageModel. Their are two optional parameters:

1) oversample: specifies the oversampling factor to be used in the model fit
   (default value = 3)
2) rotation: specifies an initial guess, in degrees, for the rotation of the
   PSF in the input image (default value = 0.0)

Output
------
The ami_analyze step produces a single output file, which contains the
following list of extensions:

1) ``FIT``: a 2-D image of the fitted model
2) ``RESID``: a 2-D image of the fit residuals
3) ``CLOSURE_AMP``: table of closure amplitudes
4) ``CLOSURE_PHA``: table of closure phases
5) ``FRINGE_AMP``: table of fringe amplitudes
6) ``FRINGE_PHA``: table of fringe phases
7) ``PUPIL_PHA``: table of pupil phases
8) ``SOLNS``: table of fringe coefficients

AMI_AVERAGE
===========

Overview
--------
The 'ami_average' step averages the results of LG processing (via the
ami_analyze step) for multiple exposures of a given target. It averages all
8 components of the ami_analyze output files for all input exposures.

Inputs
------
The only input to the 'ami_average' step is a list of input files to be
processed. These will presumably be output files from the 'ami_analyze' step.
The step has no other required or optional parameters, nor does it use any
reference files.

Output
------
The step produces a single output file, having the same format as the output
files from the 'ami_analyze' step, where the data for the 8 file components
are the average of each component from the list of input files.

AMI_NORMALIZE
=============

Overview
--------
The 'ami_normalize' step provides normalization of LG processing results for
a science target using AMI observations of a reference target that have also
had LG processing applied to them. The algorithm currently subtracts the
reference target closure phases from the science target closure phases and
divides the science target fringe amplitudes by the reference target fringe
amplitudes.

Inputs
------
The 'ami_normalize' step takes two input files: the first is the LG
processed results for a science target and the second is the LG processed
results for the reference target. There are no optional parameters and no
reference files are used.

Output
------
The output is a new LG product for the science target in which the closure
phases and fringe amplitudes have been normalized using the reference target
closure phases and fringe amplitudes. The remaining components of the science
target data model are left unchanged.

CALWEBB_AMI3 Pipeline
=======================

Overview
--------
The 'calwebb_ami3' pipeline module can be used to apply all 3 steps of AMI
processing to an association of AMI exposures. The processing flow through the
'calwebb_ami3' pipeline is as follows:

1) Apply the 'ami_analyze' step to all science target AMI exposures listed in
   the input association table. Output files will have a file name suffix of
   ``_lg.``

2) Apply the 'ami_average' step to combine the above results for the science
   target exposures. The output file will have a file name suffix of
   ``_lgavg.``

3) If there are reference target AMI exposures listed in the input association
   table, apply the 'ami_analyze' step to each of them. Output files will have
   a file name suffix of ``lg.``

4) If 'ami_analyze' was applied to reference target exposures, apply the
   'ami_average' step to combine the results for the reference target. The 
   output file will have a file name suffix of ``_lgavg.``

5) If reference target results exist, apply the 'ami_normalize' step to the
   averaged science target result, using the averaged reference target result.
   The output file will have a file name suffix of ``_lgnorm.``

Input
-----
The only input to the 'calwebb_ami3' pipeline is the name of a json-formatted
association table. There are no optional parameters.

