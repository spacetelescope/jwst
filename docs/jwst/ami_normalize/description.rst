Description
-----------

:Class: `jwst.ami.AmiNormalizeStep`
:Alias: ami_normalize

The ``ami_normalize`` step is one of the AMI-specific steps in the ``ami``
sub-package and is used in Stage 3 :ref:`calwebb_ami3 <calwebb_ami3>`
processing. It provides normalization of LG processing results for
a science target using LG results of a reference PSF target. The algorithm
subtracts the PSF target closure phases from the science target closure
phases and divides the science target fringe amplitudes by the PSF target
fringe amplitudes.

Arguments
---------
The ``ami_normalize`` step does not have any step-specific arguments.

Inputs
------

Interferometric Observables
^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiOIModel`
:File suffix: _ami-oi.fits and/or _psf-ami-oi.fits

The ``ami_normalize`` step takes two inputs: the first is the 
interferometric observables for a science target and the second 
is the interferometric observables for the PSF target. These should
be the _ami-oi.fits and _psf-ami-oi.fits products resulting from the
:ref:`ami_analyze <ami_analyze_step>` step, or two _ami-oi.fits files if the steps 
are run independently. The inputs can be in the form of filenames or 
`~jwst.datamodels.AmiOIModel` data models.

Outputs
-------

Normalized Interferometric Observables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiOIModel`
:File suffix: _aminorm-oi.fits

The output is a new set of interferometric observables for the science target
in which the closure phases and fringe amplitudes have been normalized using the PSF target
closure phases and fringe amplitudes. The remaining components of the science
target data model are left unchanged. The output file name syntax is source-based,
using the product name specified in the input ASN file and having a product type
of "_aminorm-oi.fits", e.g. "jw87600-a3001_t001_niriss_f480m-nrm_aminorm-oi.fits."

Reference Files
---------------
The ``ami_normalize`` step does not use any reference files.
