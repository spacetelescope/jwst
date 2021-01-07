.. _ipc_reffile:

IPC Reference File
------------------

:REFTYPE: IPC
:Data model: `~jwst.datamodels.IPCModel`

The IPC reference file contains a deconvolution kernel.

.. include:: ../references_general/ipc_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for IPC
++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in IPC reference files,
because they are used as CRDS selectors
(see :ref:`ipc_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
=========  ==============================

Reference File Format
+++++++++++++++++++++
IPC reference files are FITS format, with 1 IMAGE extension.
The FITS primary HDU does not contain a data array.
The format and content of the file can be one of two forms,
as described below:

=======  ========  =====  =============================  =========
EXTNAME  XTENSION  NAXIS  Dimensions                     Data type
=======  ========  =====  =============================  =========
SCI      IMAGE       2    nkern x nkern                  float
or
SCI      IMAGE       4    ncols x nrows x nkern x nkern  float
=======  ========  =====  =============================  =========

Two formats are currently supported for the IPC kernel: a small 2-D array
or a 4-D array.  If the kernel is 2-D, its dimensions should be odd,
for example 3 x 3 or 5 x 5 pixels.  The value at the center pixel will be
larger than 1 (e.g. 1.02533) and the sum of all pixel values will be
equal to 1.

A 4-D kernel may be used to allow the IPC correction to vary from pixel
to pixel across the image.  In this case, the axes that are most rapidly
varying (the last two in Python notation; the first two in IRAF/FITS notation)
have dimensions equal to those of a full-frame image.  At each point in
that image, there will be a small, 2-D kernel as described in the previous
paragraph.
