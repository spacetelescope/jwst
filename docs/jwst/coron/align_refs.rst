.. _align_refs_step:

align_refs
==========

Description
-----------
The ``align_refs`` step is used to compute offsets between science target
images and reference PSF images, and shift the PSF images into
alignment. This is performed on a per-integration basis for both the science target
data and the reference PSF data. Each integration contained in the stacked PSF data
(the result of the :ref:`stack_refs <stack_refs_step>`) step is
aligned to each integration within a given science target exposure.
This results in a new product for each science target exposure that contains a stack
of individual PSF images that have been aligned to each integration in the science
target exposure.

Shifts between each PSF and target image are computed using the
``scipy.optimize.leastsq`` function. A 2D mask, supplied via a PSFMASK reference file, 
is used to indicate pixels to ignore when performing the minimization in the
``leastsq`` routine. The mask acts as a weighting function in performing the fit.
Alignment of a PSF image is performed using the ``scipy.ndimage.fourier_shift``
function and the computed sub-pixel offsets.

Arguments
---------
The ``align_refs`` step does not have any step-specific arguments.

Inputs
------

The ``align_refs`` step takes 2 inputs: a science target exposure containing a 3D
stack of calibrated per-integration images and a "_psfstack" product containing a 3D
stack of reference PSF images.

3D calibrated images
^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _calints

One of the science target exposures specified in the ASN file used as input to the
:ref:`calwebb_coron3 <calwebb_coron3>` pipeline. This should be a "_calints" product
from the :ref:`calwebb_image2 <calwebb_image2>` pipeline and contains a 3D stack of
per-integration images.

3D stacked PSF images
^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _psfstack

A "_psfstack" product created by the :ref:`stack_refs <stack_refs_step>` step, which
contains the collection of all PSF images to be used, in the form of a 3D image stack.

Outputs
-------

4D aligned PSF images
^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.QuadModel`
:File suffix: _psfalign

The output is a 4D data model, where the 3rd axis has length equal to the total number of
reference PSF images in the input PSF stack and the 4th axis has length equal to the number
of integrations in the input science target product (ncols x nrows x npsfs x nints).
Image[n,m] in the 4D data is the n :sup:`th` PSF image aligned to the m :sup:`th` science
target integration. The file name is exposure-based, using the input science target exposure
name as the root, with the addition of the association candidate ID and the "_psfalign"
product type suffix, e.g. "jw8607342001_02102_00001_nrcb3_a3001_psfalign.fits."

Reference Files
---------------
The ``align_refs`` step uses a PSFMASK reference file.

.. include:: ../references_general/psfmask_reffile.inc

.. automodapi:: jwst.coron.align_refs_step
