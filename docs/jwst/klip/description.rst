Description
-----------
The ``klip`` step is one of the coronagraphic-specific steps in the ``coron``
sub-package and is used in Stage 3 :ref:`calwebb_coron3 <calwebb_coron3>` processing.
It applies the Karhunen-Loeve Image Plane (KLIP) algorithm to coronagraphic
images, using an accompanying set of reference PSF images, in order to fit and subtract an
optimal PSF from a source image. The KLIP algorithm uses a KL decomposition of the set of
reference PSF's, and generates a model PSF from the projection of the target on the KL vectors.
The model PSF is then subtracted from the target image (Soummer, Pueyo, and Larkin 2012).
KLIP is a Principle Component Analysis (PCA) method and is very similar to the Locally
Optimized Combination of Images (LOCI) method. The main advantages of KLIP over LOCI are
the possibility of direct forward modeling and a significant speed increase.

The KLIP algorithm consists of the following high-level steps:

1) Partition the target and reference PSF images in a set of search areas, and
   subtract their average values so that they have zero mean
2) Compute the KL transform of the set of reference PSF's
3) Choose the number of modes to keep in the estimated target PSF
4) Compute the best estimate of the target PSF from the projection of the
   target image on the KL eigenvectors
5) Calculate the PSF-subtracted target image

Arguments
---------
The ``klip`` step has one optional argument:

``--truncate``
  This is an integer parameter with a default value of 50 and is used to specify the number
  of KL transform rows to keep when computing the PSF fit to the target.

Inputs
------
The ``klip`` step takes two inputs: a science target exposure in the form of a 3D data
cube and a 4D aligned PSF image ("_psfalign") product.

3D calibrated images
^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _calints

A 3D calibrated science target product containing a stack of per-integration images.
This should be a "_calints" product created by the :ref:`calwebb_image2 <calwebb_image2>`
pipeline. Normally one of the science target exposures specified in the ASN file used
as input to the :ref:`calwebb_coron3 <calwebb_coron3>` pipeline.

4D aligned PSF images
^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.QuadModel`
:File suffix: _psfalign

A 4D collection of PSF images that have been aligned to each of the per-integration images
contained in the science target "_calints" product, created by the
:ref:`align_refs <align_refs_step>` step.

Outputs
-------

3D PSF-subtracted images
^^^^^^^^^^^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _psfsub

The output is a 3D stack of PSF-subtracted images of the science target, having the same
dimensions as the input science target ("_calints") product. The PSF fitting and subtraction
has been applied to each integration image independently. The file name syntax is
exposure-based, using the root of the input "_calints" product, with the addition of the
association candidate ID and the "_psfsub" product type suffix, e.g.
"jw8607342001_02102_00001_nrcb3_a3001_psfsub.fits."

Reference Files
---------------
The ``klip`` step does not use any reference files.
