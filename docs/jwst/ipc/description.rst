Description
===========

The ``ipc`` step corrects a JWST exposure for interpixel capacitance by
convolving with an IPC reference image.

The current implementation uses an IPC reference file that is normally
a small, rectangular image (e.g. 3 x 3 pixels), a deconvolution kernel.
The kernel may, however, be a 4-D array (e.g. 3 x 3 x 2048 x 2048),
to allow the IPC correction to vary across the detector.

For each integration in the input science data, the data are corrected
group-by-group by convolving with the kernel.  Reference pixels are not
included in the convolution; that is, their values will not be changed,
and when the kernel overlaps a region of reference pixels, those pixels
contribute a value of zero to the convolution.  The ERR and DQ arrays
will not be modified.

Subarrays
=========

Subarrays are treated the same as full-frame data, with the exception
that the reference pixels may be absent.
