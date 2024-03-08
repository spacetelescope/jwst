.. _nrm_reffile:

NRM Reference File
^^^^^^^^^^^^^^^^^^

:REFTYPE: NRM
:Data model: `~jwst.datamodels.NRMModel`

The NRM reference file contains a 2-D model of the pupil mask.

.. include:: ../references_general/nrm_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for NRM
++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in NRM reference files,
because they are used as CRDS selectors
(see :ref:`nrm_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================


Reference File Format
+++++++++++++++++++++
NRM reference files are FITS format, with one IMAGE extension.
The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
NRM      IMAGE       2    1024 x 1024     float
=======  ========  =====  ==============  =========

The ``NRM`` array contains the 2-D image representing the geometry
of the non-redundant mask in the pupil wheel of NIRISS. This mask 
enables the Aperture Masking Interferometry (AMI) mode.
