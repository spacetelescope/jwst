.. _pictureframe_reffile:

PICTUREFRAME Reference File
---------------------------

:REFTYPE: PICTUREFRAME
:Data model: `~stdatamodels.jwst.datamodels.PictureFrameModel`

The PICTUREFRAME reference file contains data necessary for removing
thermal artifacts (the "picture frame" effect) from NIRSpec exposures.

.. include:: ../includes/standard_keywords.inc

PICTUREFRAME Reference File Format
++++++++++++++++++++++++++++++++++

NIRSpec PICTUREFRAME reference files are in FITS format with one
image extension. The primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    2048 x 2048     float
=======  ========  =====  ==============  =========

The values in the ``SCI`` array give rate values for the thermal background
artifact in full frame exposures.  The reference image is intended to be scaled
and offset before subtraction from science exposure images.
