.. _drizpars_reffile:

DRIZPARS Reference File
-----------------------

:REFTYPE: DRIZPARS
:Data model: `~jwst.datamodels.DrizParsModel`

The DRIZPARS reference file contains various drizzle parameter values that
control the characteristics of a drizzled image and how it is built.

.. include:: ../references_general/drizpars_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for DRIZPARS
+++++++++++++++++++++++++++++++++++++
No additional specific keywords are required in DRIZPARS reference
files, because CRDS selection is based only on the instrument name
(see :ref:`drizpars_selectors`).

Reference File Format
+++++++++++++++++++++
DRIZPARS reference files are FITS format, with 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

========  ========  =====  =======================  =========
EXTNAME   XTENSION  NAXIS  Dimensions               Data type
========  ========  =====  =======================  =========
DRIZPARS  BINTABLE    2    TFIELDS = 7              N/A
========  ========  =====  =======================  =========

The DRIZPARS extension contains various step parameter values to be
used when processing certain types of image collections. The first
two columns (numimages and filter) are used as row selectors within
the table. Image collections that match those selectors then use the
parameter values specified in the remainder of that table row.
The table contains the following 7 columns:

===========  =======  ===========================================
TTYPE        TFORM    Description
===========  =======  ===========================================
numimages    integer  The number of images to be combined
filter       integer  The filter used to obtain the images
pixfrac      float    The pixel "shrinkage" fraction
kernel       string   The kernel function used to distribute flux
fillval      float    Value assigned to pixels with no input flux
wht_type     sting    The input image weighting type
stepsize     integer  Output WCS grid interpolation step size
===========  =======  ===========================================
