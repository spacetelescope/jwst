.. _psf_reffile:

PSF Reference File
^^^^^^^^^^^^^^^^^^

:REFTYPE: PSF
:Data model: `~jwst.datamodels.SpecPsfModel`

The PSF reference file contains a model of the 1-D point spread function
by wavelength, intended to support spectral modeling and extraction.


Reference Selection Keywords for PSF
++++++++++++++++++++++++++++++++++++
CRDS selects appropriate PSF references based on the following keywords.
PSF is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== =========================================================================
Instrument Keywords
========== =========================================================================
MIRI       INSTRUME, DETECTOR, FILTER, EXP_TYPE
========== =========================================================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
PSF reference files are in FITS format, with 2 IMAGE extensions.
The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
PSF      IMAGE       2    ncols x nrows   float
WAVE     IMAGE       1    ncols           float
=======  ========  =====  ==============  =========

The values in the ``PSF`` array give relative spectral flux values by cross-dispersion
position, at each dispersion element specified in the ``WAVE`` array.  Detector
pixels are subsampled by the amount specified in the SUBPIX keyword, and the PSF
is centered on the cross-dispersion element specified in either CENTCOL (for vertical
dispersion) or CENTROW (for horizontal dispersion).
