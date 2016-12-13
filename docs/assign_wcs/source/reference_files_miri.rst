Reference File Format for MIRI 
------------------------------  

MIRI Imaging Mode 
::::::::::::::::: 

There are two reference files required: distortion and filteroffset.

Distortion
~~~~~~~~~~

Required Fields:      model  


Filter Offset 
~~~~~~~~~~~~~  

The filter offset reference file must be an ASDF file that contains a dictionary of row and column offsets for the MIRI imaging dataset. The filter offset reference file must contain a dictionary in the tree that is indexed by the instrument filter.  The dictionary must contain two fields needed from the filter offset reference file: row_offset and column_offset and must be in units of mm.  

Required Fields:      row_offset     column_offset  

MIRI LRS Mode 
:::::::::::::  

There are two reference files required: distortion and specwcs.  

Distortion 
~~~~~~~~~~ 

 model  

SpecWCS 
~~~~~~~ 

"data"             imx  /  imxsltl             imy  /  imysltl  

MIRI IFU 
::::::::  

There are 5 reference files required: disortion, specwcs, regions, wavelengthrange and v2v3.

Distortion
~~~~~~~~~~

alpha_model             beta_model             x_model             y_model             slice-model  

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
