Reference Files
===============

By default the ``straylight`` step does not use a reference file.
It instead uses meta information added by the ``assign_wcs`` step that
maps each pixel to an IFU slice number or the regions between slices
within the 2-D slice image.

If the step option ``method`` is set to "Nearest", a more simplistic
algorithm is used to compute the correction, which requires the use of
the ``STRAYMASK`` reference file. It contains a simple mask indicating
which pixels lie within the inter-slice regions of the image.

STRAYMASK Reference File
------------------------

:REFTYPE: STRAYMASK
:Data model: `~jwst.datamodels.StrayLightModel`

The STRAYMASK reference file contains a simple mask that indicates
which pixels lie within IFU slices and which ones lie between slices.

.. include:: straymask_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for STRAYMASK
++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in STRAYMASK reference files,
because they are used as CRDS selectors
(see :ref:`straymask_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
BAND       model.meta.instrument.band  
=========  ==============================

Reference File Format
+++++++++++++++++++++
STRAYMASK reference files are FITS format with 1 IMAGE extension.
The FITS primary HDU does not contain a data array.
The characteristics of the extension is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
MASK     IMAGE       2    ncols x nrows   integer
=======  ========  =====  ==============  =========

The ``MASK`` array provides a simple 0 or 1 mapping to indicate which pixels
fall in IFU slice gaps, with 1 indicating a gap pixel and 0 for non-gap pixels.
The simple form of the ``straylight`` correction step will use only the pixels
flagged with a 1.

