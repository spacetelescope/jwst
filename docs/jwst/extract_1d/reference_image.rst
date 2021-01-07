Reference Image Format
======================
An alternative EXTRACT1D reference format, an image, is also supported.
There are currently no files of this type in CRDS (there would be a conflict
with the current JSON-format reference files), but a user can create a file
in this format and specify that it be used as an override for the default
EXTRACT1D reference file.

This format is a `~jwst.datamodels.MultiExtract1dImageModel`, which is
loosely based on `~jwst.datamodels.MultiSlitModel`.  The file should
contain keyword DATAMODL, with value 'MultiExtract1dImageModel'; this is
not required, but it makes it possible to open the file simply with
`datamodels.open`.  The reference image file contains one or more images,
which are of type `~jwst.datamodels.Extract1dImageModel`, and one can
iterate over the list of these images to find one that matches the
observing configuration.  This iterable is the ``images`` attribute of
the model (``ref_model``, for purposes of discussion).  Each element of
``ref_model.images`` can contain a ``name`` attribute (FITS keyword
SLTNAME) and a ``spectral_order`` attribute (FITS keyword SPORDER), which
can be compared with the slit name and spectral order respectively in the
science data model in order to select the matching reference image.  The
wildcard for SLTNAME is "ANY", and any integer value for SPORDER greater
than or equal to 1000 is a wildcard for spectral order (SPORDER is an
integer, and an integer keyword may not be assigned a string value such as
"ANY").  For IFU data, the image to use is selected only on ``name``.

For non-IFU data, the shape of the reference image should match the shape
of the science data, although the step can either trim the reference image
or pad it with zeros to match the size of the science data, pinned at
pixel [0, 0].  For IFU data, the shape of the reference image can be 3-D,
exactly matching the shape of the IFU data, or it can be 2-D, matching
the shape of one plane of the IFU data.  If the reference image is 2-D,
it will be applied equally to each plane of the IFU data, i.e. it will be
broadcast over the dispersion direction.

The data type of each image is float32, but the data values may only
be +1, 0, or -1.  A value of +1 means that the matching pixel in the
science data will be included when accumulating data for the source
(target) region.  A value of 0 means the pixel will not be used for
anything.  A value of -1 means the pixel will be included for the
background; if there are no pixels with value -1, no background will be
subtracted.  A pixel will either be included or not; there is no option
to include only a fraction of a pixel.

For non-IFU data, values will be extracted column by column (if the
dispersion direction is horizontal, else row by row).  The gross count
rate will be the sum of the source pixels in a column (or row).  If
background region(s) were specified, the sum of those pixels will be
scaled by the ratio of the number of source pixels to the number of
background pixels (with possibly a different ratio for each column (row))
before being subtracted from the gross count rate.  The scaled background
is what will be saved in the output table.

For IFU data, the values will be summed over each plane in the dispersion
direction, giving one value of flux and optionally one value of background
per plane.  The background value will be scaled by the ratio of source
pixels to background pixels before being subtracted from the flux.
