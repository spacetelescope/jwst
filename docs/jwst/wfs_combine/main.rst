Description
============

The ``wfs_combine`` step combines a pair of dithered Wavefront Sensing and Control (WFS&C) images.
The input images are aligned with one another and then combined using a pixel
replacement technique, described in detail below. The images are aligned to only the nearest
integer pixel in each direction. No sub-pixel resampling is done.

Algorithm
---------
Creation of the output combined image is a three-step process: first the offsets between the images
are computed, the offsets are used to shift image 2 to be in alignment with image 1, and finally
the aligned data from the two images are combined.

Computing Offsets
^^^^^^^^^^^^^^^^^
The WCS transforms of each image are used to compute the RA/Dec values for the center pixel
in image 1, and then the pixel indexes of those RA/Dec values are computed in image 2. The
difference in the pixel indexes, rounded to the nearest whole pixel, is used as the nominal
offsets in the x/y image axes.

If the optional argument "--do_refine" is set to ``True``, the nominal offsets are emperically
refined using a cross-correlation technique. The steps in the refinement are as follows:

1. Create a smoothed version of image 1 using a Gaussian kernel.
2. Find the approximate centroid of the source in image 1 by computing the mean pixel coordinates,
   separately in the x and y axes, of all pixel values that are above 50% of the peak signal
   in the smoothed image.
3. Create subarrays from image 1 and 2 centered on the computed source centroid.
4. For a range of +/- 2 pixels on either side of the nominal offsets computed from the WCS info,
   compute the cross-correlation between the two images.
5. Use the pixel offsets corresponding to the cross-correlation minimum as delta offsets to add
   to the nominal offsets computed from the WCS info.

Creating the Combined Image
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The image 2 data are shifted using the pixel offsets computed above, in order to align it with
image 1. For each pixel location in image 1, the output combined image is populated using the
following logic:

1. If the pixel values in both image 1 and 2 are good, i.e. DQ=0, the output SCI and ERR image
   values are the average of the input SCI and ERR values, respectively, and the output DQ is
   set to 0.

2. If the image 1 pixel is bad (DQ>0) and the image 2 pixel is good, the output SCI and ERR image
   values are copied from image 2, and the output DQ is set to 0.

3. If the image 1 pixel is good (DQ=0) and the image 2 pixel is bad, the output SCI and ERR image
   values are copied from image 1, and the output DQ is set to 0.

4. If both image 1 and 2 pixels are bad (DQ>0), the output SCI and ERR image values are set to
   0 and the output DQ contains the combination of input DQ values, as well as the "DO_NOT_USE"
   flag.

Upon successful completion of this step, the status keyword S_WFSCOM will be set to "COMPLETE"
in the output image header.

Inputs
------

2D calibrated images
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _cal

The input to ``wfs_combine`` is a pair of calibrated ("_cal") exposures, specified
via an ASN file. The ASN file may contain a list of several combined products to be created, in
which case the step will loop over each set of inputs, creating a combined output for each pair.

Outputs
-------

2D combined image
^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _wfscmb

The output is the combined image, using the product type suffix "_wfscmb."
