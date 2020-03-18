Description
-----------
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

LG model parameters
^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiLgModel`
:File suffix: _amiavg and _psf-amiavg

The ``ami_normalize`` step takes two inputs: the first is the LG results for
a science target and the second is the LG results for the PSF target. These should
be the "_amiavg" and "_psf-amiavg" products resulting from the
:ref:`ami_average <ami_average_step>` step. The inputs can be in the form of file
names or `~jwst.datamodels.AmiLgModel` data models.

Outputs
-------

Normalized LG model parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.AmiLgModel`
:File suffix: _aminorm

The output is a new LG product for the science target in which the closure
phases and fringe amplitudes have been normalized using the PSF target
closure phases and fringe amplitudes. The remaining components of the science
target data model are left unchanged. The output file name syntax is source-based,
using the product name specified in the input ASN file and having a product type
of "_aminorm", e.g. "jw87600-a3001_t001_niriss_f480m-nrm_aminorm.fits."

Reference Files
---------------
The ``ami_normalize`` step does not use any reference files.
