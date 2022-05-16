Description
-----------

:Class: `jwst.coron.StackRefsStep`
:Alias: stack_refs

The ``stack_refs`` step is one of the coronagraphic-specific steps in the
``coron`` sub-package and is part of Stage 3 :ref:`calwebb_coron3 <calwebb_coron3>`
processing. It takes a list of reference PSF products and stacks all of the
per-integration images contained in each PSF product into a single 3D data cube.
This operation prepares the PSF images for use by subsequent steps in the
:ref:`calwebb_coron3 <calwebb_coron3>` pipeline. The image data are simply copied
and reformatted, without being modified in any way.

Arguments
---------
The ``stack_refs`` step does not have any step-specific arguments.

Inputs
------

3D calibrated images
^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _calints

The inputs to the ``stack_refs`` step are multiple calibrated products for the PSF
target, produced by the :ref:`calwebb_image2 <calwebb_image2>` pipeline. Each input
should be a 3D "_calints" product, containing a 3D stack of calibrated images for the
multiple integrations within each exposure.

It is assumed that the ``stack_refs`` step will be called from the
:ref:`calwebb_coron3 <calwebb_coron3>` pipeline, which is given an ASN file as input,
specifying one or more PSF target exposures.
The actual input passed to the ``stack_refs`` step will be a `~jwst.datamodels.ModelContainer`
created by the :ref:`calwebb_coron3 <calwebb_coron3>` pipeline, containing a
`~jwst.datamodels.CubeModel` data model for each PSF "_calints" exposure listed in the
ASN file. See :ref:`calwebb_coron3 <calwebb_coron3>` for more details on the contents of
the ASN file.

Outputs
-------

3D PSF image stack
^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _psfstack

The output of the ``stack_refs`` step will be a single 3D product containing a stack of
all the PSF images from the multiple input exposures. The size of the stack will be equal
to the sum of the number of integration (NINTS) in each input PSF exposure.
The output file name is source-based, using the product name specified in the ASN file,
e.g. "jw86073-a3001_t001_nircam_f140m-maskbar_psfstack.fits."

Reference Files
---------------
The ``stack_refs`` step does not use any reference files.
