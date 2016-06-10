Tasks in the package
====================

The coronagraphic package currently consists of two prototype tasks: klip and
hlsp. The klip step applies the Karhunen-Loe`ve Image Plane (KLIP) algorithm to
a coronagraphic image to fit and subtract a PSF. The hlsp task produces
high-level science products (HLSP's) for a KLIP-subtracted images.

KLIP
====

Overview
--------

The klip task applies the KLIP algorithm to a coronagraphic image, using an
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

Input Arguments
---------------

The klip task takes two input file names, the first being the name of the
target image and the second being the file containing the reference PSF images.
The target image must be in the form of a regular 2-D ImageModel and the
reference PSF file must be in the form of a 3-D CubeModel, where each plane of
the data cube contains an individual 2-D PSF image.

The task takes one optional argument, "truncate", which is used to specify the
number of KL transform rows to keep when computing the PSF fit to the target.
The default value is 50.

Outputs
-------

The klip task produces two output images, the first being the best-fit target
PSF (file name suffix "_psf") and the second being the PSF-subtracted target
image (file name suffix "_klip"). Both are in the form of a simple 2-D 
ImageModel.

HLSP
====

Overview
--------

The hlsp task produces high-level science products for KLIP-processed images.
The task currently produces two such products: a signal-to-noise ratio (SNR)
image and a table of contrast data. The SNR image is computed by simply taking
the ratio of the SCI and ERR arrays of the input target image. The contrast
data are in the form of azimuthally-averaged noise versus radius. The noise
is computed as the 1-sigma standard deviation within a set of concentric
annuli centered in the input image. The annuli regions are computed to the
nearest whole pixel; no sub-pixel calculations are performed.

Input Arguments
---------------

The hlsp task takes one input file name argument, which is the name of the
KLIP-processed target image to be analyzed. One optional argument is available,
"annuli_width", which specifies the width (in pixels) of the annuli to use in
calculating the contrast data. The default value is 2 pixels.

Outputs
-------

The hslp task produces two output products. The first is the snr image (file
name suffix "_snr") and the second is the table of contrast data (file name
suffix "_contrast"). The contrast data are stored as a 2-column table giving
radius (in pixels) and noise (1-sigma).

