Description
-----------

:Class: `jwst.ami.AmiAverageStep`
:Alias: ami_average

The ``ami_average`` step is one of the AMI-specific steps in the ``ami``
sub-package and is part of Stage 3 :ref:`calwebb_ami3 <calwebb_ami3>` processing.
It averages the results of LG processing from the
:ref:`ami_analyze <ami_analyze_step>` step for multiple exposures of a given target.
It computes a simple average for all 8 components of the "ami" product files from all
input exposures.

For a given association of exposures, the "ami" products created by the `ami_analyze`
step may have `fit_image` and `resid_image` images that vary in size from one
exposure to another. If this is the case, the smallest image size of all the input
products is used for the averaged product and the averaged `fit_image` and
`resid_image` images are created by trimming extra rows/columns from the edges of
images that are larger.

Arguments
---------
The ``ami_average`` step does not have any step-specific arguments.

Inputs
------

LG model parameters
^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiLgModel`
:File suffix: _ami

The only input to the ``ami_average`` step is a list of one or more "ami" files to be
processed. These should be output files from the
:ref:`ami_analyze <ami_analyze_step>` step. The input to the step must be in the form
of a list of "ami" **file names**. Passing data models or ASN files is not supported
at this time. Use the :ref:`calwebb_ami3 <calwebb_ami3>` pipeline to conveniently
process multiple inputs.

Outputs
-------

Average LG model parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiLgModel`
:File suffix: _amiavg or _psf-amiavg

The ``ami_average`` step produces a single output file, having the same format as the input
files, where the data for the 8 file components are the averages from the list of input files.
If the inputs in the ASN file are designated as "science", the output product type will be
"_amiavg", whereas if the inputs are designated as "psf", the output product type will be
"_psf-amiavg." The output file name syntax is source-based, using the product name specified
in the input ASN file, e.g. "jw87600-a3001_t001_niriss_f480m-nrm_amiavg.fits."

Reference Files
---------------
The ``ami_average`` step does not use any reference files.
