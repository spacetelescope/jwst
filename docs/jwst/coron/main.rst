Tasks in the package
====================
The coronagraphic package currently consists of the following tasks:

* stack_refs
* align_refs
* klip
* hlsp

Briefly, the ``stack_refs`` step is used to load images of reference PSF
targets, as listed in an Association file, and stack the images into
a data cube in a single file to be used in subsequent processing steps.
The ``align_refs`` step is then used to align the stacked reference PSF
images with the images contained in a science target exposure.
The ``klip`` step applies the Karhunen-Loeve Image Plane (KLIP) algorithm to
the aligned reference PSF and science target images and produces
PSF-subtracted science target images. The ``hlsp`` task produces
high-level science products (HLSP's) from a KLIP-subtracted image.

CALWEBB_CORON3
==============
Currently the individual steps can only be run in a convenient way by
running the ``calwebb_coron3`` pipeline, which calls the individual steps
and takes care of all the necessary loading and passing of data
models for the input and output products of each step. The input to
the ``calwebb_coron3`` pipeline is expected to be an ASN file. The ASN file
should define a single output product, which will be the combined image
formed from the PSF-subtracted results of all the input science target
data. That output product should then define, as its members, the
various input reference PSF and science target files to be used in the
processing. An example ASN file is shown below.

::

 {"asn_rule": "CORON", "target": "NGC-3603", "asn_pool": "jw00017_001_01_pool", "program": "00017",
 "products": [
     {"prodtype": "coroncmb", "name": "jw89001-c1001_t001_nircam_f160w",
      "members": [
          {"exptype": "science", "expname": "targ1_calints.fits"},
          {"exptype": "science", "expname": "targ2_calints.fits"},
          {"exptype": "psf", "expname": "psf1_calints.fits"},
          {"exptype": "psf", "expname": "psf2_calints.fits"},
          {"exptype": "psf", "expname": "psf3_calints.fits"}]}],
 "asn_type": "coron",
 "asn_id": "c1001"}

In this example the output product "jw89001-c1001_t001_nircam_f160w"
is defined to consist of 2 science target inputs and 3 reference PSF
inputs. Note that the values of the ``exptype`` attribute for each
member are very important and used by the ``calwebb_coron3`` pipeline to
know which members are to be used as reference PSF data and which are
data for the science target. The output product name listed in the ASN
file is used as the root name for some of the products created by the
``calwebb_coron3`` pipeline. This includes:

- rootname_psfstack: the output of the ``stack_refs`` step
- rootname_i2d: the final combined target image

Other products will be created for each individual science target
member, in which case the root names of the original input science
target products will be used as a basis for the output products.
These products include:

- targetname_psfalign: the output of the ``align_refs`` step
- targetname_psfsub: the output of the ``klip`` step

Stack_refs
==========

Overview
--------

The ``stack_refs`` step takes a list of reference PSF products and stacks all
of the images in the PSF products into a single 3D data cube. It is
assumed that the reference PSF products are in the form of a data cube
(jwst CubeModel type data model) to begin with, in which images from
individual integrations are stacked along the 3rd axis of the data cube.
Each data cube from an input reference PSF file will be appended to a
new output 3D data cube (again a CubeModel), such that the dimension of
the 3rd axis of the output data cube will be equal to the total number
of integrations contained in all of the input files.

Inputs and Outputs
------------------
The ``stack_refs`` step is called from the ``calwebb_coron3`` pipeline module.
The ``calwebb_coron3`` pipeline will find all of the `psf` members listed
in the input ASN file, load each one into a CubeModel data model, and
construct a ModelContainer that is the list of all psf CubeModels. The
ModelContainer is passed as input to the ``stack_refs`` step. The output
of ``stack_refs`` will be a single CubeModel containing all of the
concatenated data cubes from the input psf files.

.. automodapi:: jwst.coron.stack_refs_step

Align_refs
==========

Overview
--------

The ``align_refs`` step is used to compute offsets between science target
images and the reference PSF images and shift the PSF images into
alignment. Each integration contained in the stacked PSF data is
aligned to each integration within a given science target product.
The ``calwebb_coron3`` pipeline applies the ``align_refs`` step to each input
science target product individually, resulting in a set of PSF images
that are aligned to the images in that science target product.

Inputs and Outputs
------------------
The ``align_refs`` step takes 2 inputs: a science target product, in the
form of a CubeModel data model, and the stacked PSF product, also in
the form of a CubeModel data model. The resulting output is a 4D data
model (QuadModel), where the 3rd axis has length equal to the total
number of reference PSF images in the input PSF stack and the 4th
axis has length equal to the number of integrations in the input
science target product.

.. automodapi:: jwst.coron.align_refs_step

Klip
====

Overview
--------
The ``klip`` task applies the KLIP algorithm to coronagraphic images, using an
accompanying set of reference PSF images, in order to
fit and subtract an optimal PSF from the source. The KLIP algorithm uses a KL
decomposition of the set of reference PSF's, and generates a model PSF from the
projection of the target on the KL vectors. The model PSF is then subtracted
from the target image (Soummer, Pueyo, and Larkin 2012). KLIP is a
Principle Component Analysis (PCA) method and is very similar to LOCI. The
main advantages of KLIP over LOCI is the possibility of direct forward
modeling and a significant speed increase.

The KLIP algorithm consists of the following steps:

1) Partition the target and reference images in a set of search areas, and
   subtract their average values so that they have zero mean.
2) Compute the KL transform of the set of reference PSF's
3) Choose the number of modes to keep in the estimated target PSF
4) Compute the best estimate of the target PSF from the projection of the
   target image on the KL eigenvectors
5) Calculate the PSF-subtracted target image

Inputs and Outputs
------------------
The ``klip`` task takes two inputs: a science target product, in the form of a 3D
CubeModel data model, and a set of aligned PSF images, in the form of a 4D
QuadModel data model. Each 'layer' in the 4th dimension of the PSF data
contains all of the aligned PSF images corresponding to a given integration
(3rd dimension) in the science target cube. The output from the klip step is
a 3D CubeModel data model, having the same dimensions as the input science
target product, and contains the PSF-subtracted images for every integration
of the science target product.

Arguments
---------

The task takes one optional argument, `truncate`, which is used to specify the
number of KL transform rows to keep when computing the PSF fit to the target.
The default value is 50.

.. automodapi:: jwst.coron.klip_step

HLSP
====

Overview
--------
The ``hlsp`` task produces high-level science products for KLIP-processed images.
The task currently produces two such products: a signal-to-noise ratio (SNR)
image and a table of contrast data. The SNR image is computed by simply taking
the ratio of the SCI and ERR arrays of the input target image. The contrast
data are in the form of azimuthally-averaged noise versus radius. The noise
is computed as the 1-sigma standard deviation within a set of concentric
annuli centered in the input image. The annuli regions are computed to the
nearest whole pixel; no sub-pixel calculations are performed.

Input Arguments
---------------
The ``hlsp`` task takes one input file name argument, which is the name of the
KLIP-processed target product to be analyzed. One optional argument is available,
`annuli_width`, which specifies the width (in pixels) of the annuli to use in
calculating the contrast data. The default value is 2 pixels.

Outputs
-------
The ``hslp`` task produces two output products. The first is the snr image (file
name suffix "_snr") and the second is the table of contrast data (file name
suffix "_contrast"). The contrast data are stored as a 2-column table giving
radius (in pixels) and noise (1-sigma).

.. automodapi:: jwst.coron.hlsp_step
