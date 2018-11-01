Reference Files
===============

The ``persistence`` step uses TRAPDENSITY_, PERSAT_, and
TRAPPARS_ reference files.

.. _TRAPDENSITY:

TRAPDENSITY Reference File
--------------------------

:REFTYPE: TRAPDENSITY
:Data model: `~jwst.datamodels.TrapDensityModel`

The TRAPDENSITY reference file contains a pixel-by-pixel map of the
trap density.

.. include:: trapdensity_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for TRAPDENSITY
++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in TRAPDENSITY reference files,
because they are used as CRDS selectors
(see :ref:`trapdensity_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
=========  ==============================

Reference File Format
+++++++++++++++++++++
TRAPDENSITY reference files are FITS format, with 2 IMAGE extensions
and 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   int
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

.. include:: ../includes/dq_def.rst

.. _PERSAT:

PERSAT Reference File
---------------------

:REFTYPE: PERSAT
:Data model: `~jwst.datamodels.PersistenceSatModel`

The PERSAT reference file contains a pixel-by-pixel map of the
persistence saturation (full well) threshold.

.. include:: persat_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for PERSAT
+++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in PERSAT reference files,
because they are used as CRDS selectors
(see :ref:`persat_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
=========  ==============================

Reference File Format
+++++++++++++++++++++
PERSAT reference files are FITS format, with 2 IMAGE extensions
and 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   int
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

.. include:: ../includes/dq_def.rst

.. _TRAPPARS:

TRAPPARS Reference File
-----------------------

:REFTYPE: TRAPPARS
:Data model: `~jwst.datamodels.TrapParsModel`

The TRAPPARS reference file contains default parameter values
used in the persistence correction.

.. include:: trappars_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for TRAPPARS
+++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in TRAPPARS reference files,
because they are used as CRDS selectors
(see :ref:`trappars_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
=========  ==============================

Reference File Format
+++++++++++++++++++++
TRAPPARS reference files are FITS format, with 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

========  ========  ===========
EXTNAME   XTENSION  Dimensions
========  ========  ===========
TRAPPARS  BINTABLE  TFIELDS = 4
========  ========  ===========

The format and contents of the table extension is as follows:

+--------------+-----------+------------------------------------------------+
| Column name  | Data type | Description                                    |
+==============+===========+================================================+
| capture0     | float     | Coefficient of exponential capture term        |
+--------------+-----------+------------------------------------------------+
| capture1     | float     | Minus the reciprocal of capture e-folding time |
+--------------+-----------+------------------------------------------------+
| capture2     | float     | The "instantaneous" capture coefficient        |
+--------------+-----------+------------------------------------------------+
| decay_param  | float     | Minus the reciprocal of decay e-folding time   |
+--------------+-----------+------------------------------------------------+

At the present time, there are no persistence reference files for MIRI
and NIRSpec. CRDS will return "N/A" for the names of the reference files
if the persistence step is run on MIRI or NIRSpec data, in which case
the input will be returned unchanged,
except that the primary header keyword S_PERSIS will will have been
set to 'SKIPPED'.

