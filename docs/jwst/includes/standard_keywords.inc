Standard Keywords
+++++++++++++++++
The following table lists the keywords that are *required* to be present in all
reference files.  The first column gives the FITS keyword names.
The second column gives the jwst data model name for each keyword, which is
useful when using data models in creating and populating a new reference file.
The third column gives the equivalent meta tag in ASDF reference file headers,
which is the same as the name within the data model meta tree (second column).

============  ==========================  ========================
FITS Keyword  Data Model Name             ASDF meta tag
============  ==========================  ========================
AUTHOR        model.meta.author           author
DATAMODL      model.meta.model_type       model_type
DATE          model.meta.date             date
DESCRIP       model.meta.description      description
FILENAME      model.meta.filename         N/A     
INSTRUME      model.meta.instrument.name  instrument: {name}
PEDIGREE      model.meta.pedigree         pedigree
REFTYPE       model.meta.reftype          reftype
TELESCOP      model.meta.telescope        telescope
USEAFTER      model.meta.useafter         useafter
============  ==========================  ========================

**NOTE:** More information on standard required keywords can be found here:
:ref:`Standard Required Keywords`
