Standard Keywords
+++++++++++++++++
The following table lists the keywords that are *required* to be present in all
reference files.  The first column gives the FITS keyword names. Note that some
reference files are in ASDF or JSON format, in which case the keyword names may
vary slightly (e.g. "INSTRUMENT" instead of "INSTRUME").
The second column gives the jwst data model
name for each keyword, which is useful when using data models in creating and
populating a new reference file.

=========  ========================
Keyword    Data Model Name
=========  ========================
AUTHOR     model.meta.author
DATAMODL   model.meta.model_type
DATE       model.meta.date
DESCRIP    model.meta.description
FILENAME   model.meta.filename
INSTRUME   model.meta.instrument.name
PEDIGREE   model.meta.pedigree
REFTYPE    model.meta.reftype
TELESCOP   model.meta.telescope
USEAFTER   model.meta.useafter
=========  ========================

**NOTE:** More information on standard required keywords can be found here:
:ref:`Standard Required Keywords`
