Reference Files
===============
The ``ami_analyze`` step uses a THROUGHPUT reference file. The ``ami_average``
and ``ami_normalize`` steps do not use any reference files.

THROUGHPUT Reference File
-------------------------

:RETYPE: THROUGHPUT
:Data model: `~jwst.datamodels.ThroughputModel`

The THROUGHPUT reference file contains throughput data for the filter used
in the AMI image.

.. include:: throughput_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for THROUGHPUT
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in THROUGHPUT reference files,
because they are used as CRDS selectors
(see :ref:`throughput_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
FILTER     model.meta.instrument.filter
=========  ==============================

Reference File Format
+++++++++++++++++++++
THROUGHPUT reference files are FITS files with one BINTABLE
extension. The FITS primary data array is assumed to be empty.
The format of the file is as follows:

==========  ========  =====  ==============  =========
EXTNAME     XTENSION  NAXIS  Dimensions      Data type
==========  ========  =====  ==============  =========
THROUGHPUT  BINTABLE    2    TFIELDS = 2     N/A
==========  ========  =====  ==============  =========

The table extension contains two columns, giving wavelength and
throughput values for a particular filter:

===========  =========  ==========
Column name  Data type  Units
===========  =========  ==========
wavelength   float      Angstroms
throughput   float      (unitless)
===========  =========  ==========
