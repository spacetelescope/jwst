Standard Keywords
+++++++++++++++++
The following table lists the keywords that are *required* to be present in all
reference files.  The first column gives the FITS keyword names (although these
reference files are ASDF).  The second column gives the model name, which is
needed when creating and populating a new reference file.

=========  ========================
Keyword    Model Name
=========  ========================
AUTHOR     meta.author
DATAMODL   meta.model_type
DATE       meta.date
DESCRIP    meta.description
FILENAME   meta.filename
INSTRUME   meta.instrument.name
PEDIGREE   meta.pedigree
REFTYPE    meta.reftype
TELESCOP   meta.telescope
USEAFTER   meta.useafter
=========  ========================

**NOTE:** More information on standard required keywords can be found here:
:ref:`Standard Required Keywords`
