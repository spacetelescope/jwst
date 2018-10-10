Reference File Types
====================
The IPC deconvolution current step uses an IPC reference file.

.. include:: ../includes/standard_keywords.rst

.. include:: ipc_selection.rst

IPC Reference File Format
--------------------------
IPC reference files are FITS files with one IMAGE extension, with EXTNAME
value of 'SCI'.  The FITS primary data array is assumed to be empty.
The SCI extension contains a floating-point data array.

Two formats are currently supported for the IPC kernel, a small 2-D array
or a 4-D array.  If the kernel is 2-D, its dimensions should be odd,
perhaps 3 x 3 or 5 x 5 pixels.  The value at the center pixel will be
larger than 1 (e.g. 1.02533), and the sum of all pixel values will be
equal to 1.

A 4-D kernel may be used to allow the IPC correction to vary from point
to point across the image.  In this case, the axes that are most rapidly
varying (the last two, in Python notation; the first two, in IRAF notation)
have dimensions equal to those of a full-frame image.  At each point in
that image, there will be a small, 2-D kernel as described in the previous
paragraph.
