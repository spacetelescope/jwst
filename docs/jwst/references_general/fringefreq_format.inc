+-------------------+--------------+-----------+------------+
| Table             | Column       | Data type | Units      |
+===================+==============+===========+============+
| CUBEPAR           | CHANNEL      | shortint  | N/A        |
+                   +--------------+-----------+------------+
|                   | BAND         | ch*6      | N/A        |
+                   +--------------+-----------+------------+
|                   | WAVEMIN      | float     | microns    |
+                   +--------------+-----------+------------+
|                   | WAVEMAX      | float     | microns    |
+                   +--------------+-----------+------------+
|                   | SPAXELSIZE   | float     | arcseconds |
+                   +--------------+-----------+------------+
|                   | SPECTRALSTEP | double    | microns    |
+-------------------+--------------+-----------+------------+
| CUBEPAR_MSM       | CHANNEL      | shortint  | N/A        |
+                   +--------------+-----------+------------+
|                   | BAND         | ch*6      | N/A        |
+                   +--------------+-----------+------------+
|                   | ROISPATIAL   | float     | arcseconds |
+                   +--------------+-----------+------------+
|                   | ROISPECTRAL  | double    | microns    |
+                   +--------------+-----------+------------+
|                   | POWER        | double    | unitless   |
+                   +--------------+-----------+------------+
|                   | SOFTRAD      | double    | unitless   |
+-------------------+--------------+-----------+------------+
| MULTICHANNEL_MSM  | WAVELENGTH   | double    | microns    |
+                   +--------------+-----------+------------+
|                   | ROISPATIAL   | double    | arcseconds |
+                   +--------------+-----------+------------+
|                   | ROISPECTRAL  | double    | microns    |
+                   +--------------+-----------+------------+
|                   | POWER        | double    | unitless   |
+                   +--------------+-----------+------------+
|                   | SOFTRAD      | double    | unitless   |
+-------------------+--------------+-----------+------------+
| CUBEPAR_EMSM      | CHANNEL      | shortint  | N/A        |
+                   +--------------+-----------+------------+
|                   | BAND         | ch*6      | N/A        |
+                   +--------------+-----------+------------+
|                   | ROISPATIAL   | float     | arcseconds |
+                   +--------------+-----------+------------+
|                   | ROISPECTRAL  | double    | microns    |
+                   +--------------+-----------+------------+
|                   | SCALERAD     | double    | unitless   |
+-------------------+--------------+-----------+------------+
| MULTICHANNEL_EMSM | WAVELENGTH   | double    | microns    |
+                   +--------------+-----------+------------+
|                   | ROISPATIAL   | double    | arcseconds |
+                   +--------------+-----------+------------+
|                   | ROISPECTRAL  | double    | microns    |
+                   +--------------+-----------+------------+
|                   | SCALERAD     | double    | unitless   |
+-------------------+--------------+-----------+------------+


+----------------------+--------------+-----------+------------+
| Table                | Column       | Data type | Units      |
+======================+==============+===========+============+
| CUBEPAR              | DISPERSER    | ch*5      | N/A        |
+                      +--------------+-----------+------------+
|                      | FILTER       | ch*6      | N/A        |
+                      +--------------+-----------+------------+
|                      | WAVEMIN      | double    | microns    |
+                      +--------------+-----------+------------+
|                      | WAVEMAX      | double    | microns    |
+                      +--------------+-----------+------------+
|                      | SPAXELSIZE   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | SPECTRALSTEP | double    | microns    |
+----------------------+--------------+-----------+------------+
| CUBEPAR_MSM          | DISPERSER    | ch*5      | N/A        |
+                      +--------------+-----------+------------+
|                      | FILTER       | ch*6      | N/A        |
+                      +--------------+-----------+------------+
|                      | ROISPATIAL   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | ROISPECTRAL  | double    | microns    |
+                      +--------------+-----------+------------+
|                      | POWER        | double    | unitless   |
+                      +--------------+-----------+------------+
|                      | SOFTRAD      | double    | unitless   |
+----------------------+--------------+-----------+------------+
| MULTICHAN_PRISM_MSM  | WAVELENGTH   | double    | microns    |
+                      +--------------+-----------+------------+
|                      | ROISPATIAL   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | ROISPECTRAL  | double    | microns    |
+                      +--------------+-----------+------------+
|                      | POWER        | double    | unitless   |
+                      +--------------+-----------+------------+
|                      | SOFTRAD      | double    | unitless   |
+----------------------+--------------+-----------+------------+
| MULTICHAN_MED_MSM    | WAVELENGTH   | float     | microns    |
+                      +--------------+-----------+------------+
|                      | ROISPATIAL   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | ROISPECTRAL  | double    | microns    |
+                      +--------------+-----------+------------+
|                      | POWER        | double    | unitless   |
+                      +--------------+-----------+------------+
|                      | SOFTRAD      | double    | unitless   |
+----------------------+--------------+-----------+------------+
| MULTICHAN_HIGH_MSM   | WAVELENGTH   | float     | microns    |
+                      +--------------+-----------+------------+
|                      | ROISPATIAL   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | ROISPECTRAL  | double    | microns    |
+                      +--------------+-----------+------------+
|                      | POWER        | double    | unitless   |
+                      +--------------+-----------+------------+
|                      | SOFTRAD      | double    | unitless   |
+----------------------+--------------+-----------+------------+
| CUBEPAR_EMSM         | DISPERSER    | ch*5      | N/A        |
+                      +--------------+-----------+------------+
|                      | FILTER       | ch*6      | N/A        |
+                      +--------------+-----------+------------+
|                      | ROISPATIAL   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | ROISPECTRAL  | double    | microns    |
+                      +--------------+-----------+------------+
|                      | SCALERAD     | double    | unitless   |
+----------------------+--------------+-----------+------------+
| MULTICHAN_PRISM_EMSM | WAVELENGTH   | double    | microns    |
+                      +--------------+-----------+------------+
|                      | ROISPATIAL   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | ROISPECTRAL  | double    | microns    |
+                      +--------------+-----------+------------+
|                      | SCALERAD     | double    | unitless   |
+----------------------+--------------+-----------+------------+
| MULTICHAN_MED_EMSM   | WAVELENGTH   | float     | microns    |
+                      +--------------+-----------+------------+
|                      | ROISPATIAL   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | ROISPECTRAL  | double    | microns    |
+                      +--------------+-----------+------------+
|                      | SCALERAD     | double    | unitless   |
+----------------------+--------------+-----------+------------+
| MULTICHAN_HIGH_EMSM  | WAVELENGTH   | float     | microns    |
+                      +--------------+-----------+------------+
|                      | ROISPATIAL   | double    | arcseconds |
+                      +--------------+-----------+------------+
|                      | ROISPECTRAL  | double    | microns    |
+                      +--------------+-----------+------------+
|                      | SCALERAD     | double    | unitless   |
+----------------------+--------------+-----------+------------+
