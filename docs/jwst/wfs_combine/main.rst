
Description
============

This step takes as input a set of 2 dithered wave front sensing images. The names of these images and the name of the output are given in an association table.  The table can contain a list of several combined products to be created (each from a separate pair of input files). Each pair of input images is 'combined' by:


1. WCS information is read from both images, from which the difference in pointings (in pixels) is calculated

2. Image #2 is aligned in the frame of image #1 using this WCS information in the input headers

3. For each pixel in the overlapped region, construct a 'combined' SCI image using:

   a) the pixel from image #1 if that pixel has a good DQ value, else
   b) the pixel from image #2 if that pixel has a good DQ value, else
   c) a default value (0).

4. For each pixel in the overlapped region, construct a 'combined' Data Quality image using:

   a) the DQ pixel from image #1 if that pixel has a good DQ value, else
   b) the DQ pixel from image #2 if that pixel has a good DQ value, else
   c) a default 'BAD_WFS' value added to the corresponding value in image #1.

5. For each pixel in the overlapped region, construct a 'combined' Error image using:

   a) the ERR pixel from image #1 if that pixel has a good DQ value, else
   b) the ERR pixel from image #2 if that pixel has a good DQ value, else
   c) a default value(0).


If the option to refine the estimate of the offsets is chosen (this is not the default) step #2 above becomes:

2.

   a) Interpolate over missing data (based on the corresponding DQ array) in both images
   b) Align these interpolated images to a common frame using the WCS information in the input headers
   c) Compare the 2 nominally aligned, interpolated images by varying the offsets to have values in the neighborhood of the nominal offsets to determine the best match.
   d) Align the original (pre-interpolated) image #2 in the frame of image #1 using this refined estimate of the offsets


Upon successful completion of this step, the status keyword S_WFSCOM will be set to COMPLETE.
