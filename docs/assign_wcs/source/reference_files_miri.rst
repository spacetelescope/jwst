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

Five models: alpha_model, beta_model, x_model, y_model, slice-model

The distortion reference file contains 5 models in order to transform the detector coordinates to sky coordinates.  The alpha and beta models convert the detector coordinates first to alpha / beta (MIRI-DD-10001-STSCI_v1.pdf) which is then the MRS-FOV.  The x_model and y_model convert the coordinates from alpha and beta to XAN and YAN (JWST-FOV, similar to V2/V3). The final slice model contains linear transformations that involve scaling, rotation and shift. 


specwcs
~~~~~~~

model  / lambda_model
The specwcs reference file contains a model that maps the detector coordinates to wavelength coordinates indexed by slice up to the alpha and beta process.

regions
~~~~~~~

The IFU takes a region reference file that defines the region over which the WCS is valid. The reference file should define a polygon and may consist of a set of X,Y coordinates that define the polygon.

wavelengthrange
~~~~~~~~~~~~~~~

Fields: wavelengthrange, channels

The wavelengthrange reference file consists of two models, one that defines the wavelength range and is indexed by 'wavelengthrange' and the second is a set of channels indexed in the file by 'channels'. The model defines, per channel, the wavelength mapping in going from alpha, beta to XAN, YAN. 

v2v3
~~~~

FIeld: model in v2v3

The model field in the tree contains N models, one per channel, that map the spatial coordinates from alpha, beta to V2, V3.
model / v2v3 model
