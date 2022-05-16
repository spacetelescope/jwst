Description
-----------

:Class: `jwst.coron.AlignRefsStep`
:Alias: align_refs

The ``align_refs`` step is one of the coronagraphic-specific steps in the ``coron``
sub-package that is part of Stage 3 :ref:`calwebb_coron3 <calwebb_coron3>` processing.
It computes offsets between science target
images and reference PSF images, and shifts the PSF images into
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
The ``align_refs`` step  has two optional arguments:

``--median_box_length`` (integer, default=3)
  The box size to use when replacing bad pixels with the median in a surrounding box.

``--bad_bits`` (string, default="DO_NOT_USE")
  The DQ bit values from the input image DQ arrays that should be considered bad
  and replaced with the median in a surrounding box. For example, setting to 
  ``"OUTLIER, SATURATED"`` (or equivalently ``"16, 2"`` or ``"18"``) will treat 
  all pixels flagged as OUTLIER or SATURATED as bad, while setting to ``""`` or  
  ``None`` will treat all pixels as good and omit any bad pixel replacement.

Inputs
------

The ``align_refs`` step takes 2 inputs: a science target exposure containing a 3D
stack of calibrated per-integration images and a "_psfstack" product containing a 3D
stack of reference PSF images. If the target or PSF images have any of the
data quality flags set to those specified by the ``bad_bits`` argument, these pixels 
are replaced with the median value of a region around the flagged data. The size of the
box region to use for the replacement can also be specified. These corrected images are 
used in the :ref:`align_refs <align_refs_step>` step and passed along for subsequent 
processing.

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
