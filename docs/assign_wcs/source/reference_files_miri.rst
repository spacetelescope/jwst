Reference File Format for MIRI
------------------------------

MIRI Imaging Mode
:::::::::::::::::

There are two reference files required: distortion and filteroffset.

Distortion
~~~~~~~~~~

The distortion reference file must be an ASDF file that contains the distortion model.  The distortion model contains one field indexed by the key "model" and returns a distortion model. The distortion model contains the forward and reverse transforms to/ from pixel and sky frames (e.g. SIP coeffs).


Filter Offset
~~~~~~~~~~~~~

The filter offset reference file must be an ASDF file that contains a dictionary of row and column offsets for the MIRI imaging dataset. The filter offset reference file must contain a dictionary in the tree that is indexed by the instrument filter.  The dictionary must contain two fields needed from the filter offset reference file: row_offset and column_offset and must be in units of mm (or is it pixels ??).


MIRI LRS Mode
:::::::::::::

There are two reference files required: distortion and specwcs.

Distortion
~~~~~~~~~~

The distortion reference file must be an ASDF file that contains the distortion model.  The distortion model contains one field indexed by the key "model" and returns a distortion model. The distortion model contains the forward and reverse transforms to/ from pixel and sky frames.

SpecWCS
~~~~~~~

The reference file contains the zero point offset for the slit relative to the full field of view.  For the Fixed Slit exposure type the fields are stored in the header of the second HDU and are indexed by 'imx' and 'imy'.  For the Slitless exposure type the fields are stored in the header of the second HDU and are indexed by 'imxsltl' and 'imysltl'.  For both of the exposure types, the zero point offset is 1 based and the X (e.g., imx) refers to the column and Y refers to the row.


MIRI IFU
::::::::

There are 5 reference files required: disortion, specwcs, regions, wavelengthrange and v2v3.

Distortion
~~~~~~~~~~

alpha_model            beta_model            x_model            y_model            slice-model

specwcs
~~~~~~~

model  / lambda_model

regions
~~~~~~~

regions

wavelengthrange
~~~~~~~~~~~~~~~

wavelengthrange

channels

v2v3
~~~~

model / v2v3 model
