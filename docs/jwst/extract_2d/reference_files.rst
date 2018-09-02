Reference Files
===============

.. include:: ../includes/standard_keywords.rst

WAVECORR Reference
------------------

To apply the Nirspec wavelength zero-point correction, this step uses the
WAVECORR reference file. The zero-point correction is applied to observations
with EXP_TYPE of "NRS_FIXEDSLT", "NRS_BRIGHTOBJ" or "NRS_MSASPEC". This is an
optional correction (on by default). It can be turned off by specifying
``apply_wavecorr=False`` when running the step.

.. include:: wavecorr_selection.rst

WAVELENGTHRANGE Reference
-------------------------

NIRCAM WFSS and NIRISS WFSS observations use the wavelengthrange reference file
in order to construct the bounding boxes around each objects orders. If a list
of ``GrismObject`` is supplied, then no reference file is neccessary.

.. include:: wavelengthrange_selection.rst


