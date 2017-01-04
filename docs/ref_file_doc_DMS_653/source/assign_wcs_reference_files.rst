Reference File Types
--------------------

WCS Reference files are in the Advanced Scientific Data Format (ASDF).
The best way to create the file is to programmatically create the model and then save it to a file.
A tutorial on creating reference files in ASDF format is available at:

https://github.com/spacetelescope/jwreftools/blob/master/docs/notebooks/referece_files_asdf.ipynb

Transforms are 0-based. The forward direction is from detector to sky.

There are 16 reference types used by assign_wcs:

===================    ==========================================================   ============================
reftype                                     description                              Instrument
===================    ==========================================================   ============================
**camera**             NIRSPEC Camera model                                          NIRSPEC
**collimator**         NIRSPEC Collimator Model                                      NIRSPEC
**disperser**          Disperser parameters                                          NIRSPEC
**distortion**         Spatial distortion model                                      MIRI, FGS, NIRCAM, NIRISS
**filteroffset**       MIRI Imager filter offsets                                    MIRI
**fore**               Transform through the NIRSPEC FORE optics                     NIRSPEC
**fpa**                Transform in the NIRSPEC FPA plane                            NIRSPEC
**ifufore**            Transform from the IFU slicer to the IFU entrance             NIRSPEC
**ifupost**            Transform from the IFU slicer to the back of the IFU          NIRSPEC
**ifuslicer**          IFU Slicer geometric description                              NIRSPEC
**msa**                Transform in the NIRSPEC MSA plane                            NIRSPEC
**ote**                Transform through the Optical Telescope Element               NIRSPEC
**specwcs**            Wavelength calibration models                                 MIRI, NIRCAM, NIRISS
**regions**            Stores location of the regions on the detector                MIRI
**v2v3**               Transform from MIRI instrument focal plane to V2V3 plane      MIRI
**wavelengthrange**    Typical wavelength ranges                                     MIRI, NIRSPEC
===================    ==========================================================   ============================


CRDS Selection Criteria For Each Reference File Type 
----------------------------------------------------


CAMERA
::::::
CAMERA reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 


COLLIMATOR
::::::::::
For NIRSPEC, COLLIMATOR reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 

DISPERSER
:::::::::
For NIRSPEC, DISPERSER reference files are currently selected based on the values of EXP_TYPE and GRATING in the input science data set. 


DISTORTION
::::::::::

For MIRI, DISTORTION reference files are currently selected based on the values of EXP_TYPE, DETECTOR, CHANNEL and BAND in the input science data set. 

For FGS, DISTORTION reference files are currently selected based on the values of EXP_TYPE and DETECTOR in the input science data set. 

For NIRCAM, DISTORTION reference files are currently selected based on the values of EXP_TYPE, DETECTOR, and CHANNEL in the input science data set. 

For NIRISS, DISTORTION reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 


FILTEROFFSET
::::::::::::
For MIRI, FILTEROFFSET reference files are currently selected based on the values of EXP_TYPE and DETECTOR in the input science data set. 


FORE
::::

For NIRSPEC, FORE reference files are currently selected based on the values of EXP_TYPE and FILTER in the input science data set. 

FPA
:::
For NIRSPEC, FPA reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 

IFUFORE
:::::::
For NIRSPEC, IFUFORE reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 


IFUPOST
:::::::
For NIRSPEC, IFUPOST reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 

IFUSLICER
:::::::::
For NIRSPEC, IFUSLICER reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 


MSA
:::
For NIRSPEC, MSA reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 

OTE
:::
For NIRSPEC, OTE reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 

SPECWCS
:::::::
For MIRI, SPECWCS reference files are currently selected based on the values of DETECTOR, CHANNEL, BAND, SUBARRAY, and EXP_TYPE in the input science data set. 

For NIRISS, SPECWCS reference files are currently selected based on the values of SUBARRAY and EXP_TYPE in the input science data set. 

REGIONS
:::::::
For MIRI, REGIONS reference files are currently selected based on the values of DETECTOR, CHANNEL, BAND, EXP_TYPE in the input science data set. 


V2V3
::::
For MIRI, V2V3 reference files are currently selected based on the values of DETECTOR, CHANNEL, BAND, EXP_TYPE in the input science data set. 


WAVELENGTHRANGE
:::::::::::::::
For NIRSPEC, WAVELENGTHRANGE reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 

For MIRI, WAVELENGTHRANGE reference files are currently selected based only on the value of EXP_TYPE in the input science data set. 



Reference File Formats For Each Reference File Type 
---------------------------------------------------

CAMERA
::::::

This reference file contains an astropy compound model made up of a polynomial models, rotation and translations. The forward direction is from the FPA to the GWA.
:model: Transform through the CAMERA.

COLLIMATOR
::::::::::

The collimator reference file contains an astropy compound model made up of a polynomial models, rotation and translations. The forward direction is from the GWA to the MSA.

:model: Transform through the COLLIMATOR.


DISPERSER
:::::::::
The disperser file contains reference data about the NIRSPEC dispersers (gratings or the prism). The reference data is described in the NIRSPEC Interface Control Document.

The following fields are common for all gratings and the prism:

:grating: Name of grating
:gwa_tiltx:
    :temperatures: Temperatures measured where the GWA sensor is
    :zeroreadings: Value of GWA sensor reading which corresponds to disperser model parameters
    :tilt_model: Model of the relation between THETA_Y vs GWA_X reading
:gwa_tilty:
    :temperatures: Temperatures measured where the GWA sensor is
    :zeroreadings: Value of GWA sensor reading which corresponds to disperser model parameters
    :tilt_model: Model of the relation between THETA_X vs GWA_Y reading
:tilt_x: Angle (in degrees) between the grating surface and the reference surface (the mirror)
:tilt_y: Angle (in degrees) between the grating surface and the reference surface (the mirror)
:theta_x: Element alignment angle in x-axis (in degrees)
:theta_y: Element alignment angle in y-axis (in degrees)
:theta_z: Element alignment angle in z-axis (in degrees)

The prism reference file has in addition the following fields:

:angle: Angle between the front and back surface of the prism (in degrees)
:kcoef: K coefficients of Selmeir equation, describing the material
:lcoef: L coefficients describing the material
:tcoef: Thermal coefficients describing the properties of the glass
:tref: Reference temperature (in K)
:pref: Reference pressure (in ATM)
:wbound: Min and Max wavelength (in meters) for which the model is valid

DISTORTION
::::::::::

The distortion reference file contains a combination of astropy models. For the MIRI Imager this file contains a polynomial and filter-dependent offsets.  For the MIRI MRS, NIRCAM, NIRISS, and FGS the model is a combination of polynomials. 
:model: Transform from detector to an intermediate frame (instrument dependent).


FILTEROFFSET
::::::::::::
The filter offset reference file is an ASDF file that contains a dictionary of row and column offsets for the MIRI imaging dataset. The filter offset reference file contains a dictionary in the tree that is indexed by the instrument filter. Each filter points to two fields - row_offset and column_offset. The format is

:miri_filter_name:
    :column_offset: Offset in x (in arcmin)
    :row_offset: Offset in y (in arcmin)

FORE
::::
The FORE reference file stores the transform through the Filter Wheel Assembly (FWA). It has two fields - “filter” and “model”. The transform through the FWA is chromatic. It is represented as a Polynomial of two variables whose coefficients are wavelength dependent. The compound model takes three inputs - x, y positions and wavelength.

:filter: Filter name.
:model: Transform through the Filter Wheel Assembly (FWA).


FPA
:::
The FPA reference file stores information on the metrology of the Focal Plane Array (FPA) which consists of two single chip arrays (SCA), named NRS1 and NRS2.

The reference file contains two fields : “NRS1” and “NRS2”. Each of them stores the transform (shift and rotation) to transform positions from the FPA to the respective SCA. The output units are in pixels.

:NRS1: Transform for the NRS1 detector.
:NRS2: Transform for the NRS2 detector.


IFUFORE
:::::::
The IFU reference file provides the parameters (Paraxial and distortions coefficients)
for the coordinate transforms from the MSA plane to the plane of the IFU slicer.
:model: Compound model, Polynomials



IFUPOST
:::::::
The IFUPOST reference file provides the parameters (Paraxial and distortions coefficients) for the coordinate transforms from the slicer plane to the MSA plane (out), that is the plane of the IFU virtual slits.

The reference file contains models made up based on an offset and a polynomial. There is a model for each of the slits and is indexed by the slit number. The models is used as part of the conversion from the GWA to slit.

:ifu_slice_number:
    :model: Polynomial and rotation models.



IFUSLICER
:::::::::
The IFUSLICER stores information about the metrology of the IFU slicer - relative positioning and size of the aperture of each individual slicer and the absolute reference with respect to the center of the field of view.
The reference file contains two fields - “data” and “model”.
The “data” field is an array with 30 rows pertaining to the 30 slices and the columns are

:data: Array with reference data for each slicer. It has 5 columns

          NO
            Slice number (0 - 29)
          x_center
            X coordinate of the center (in meters)
          y_center
            Y coordinate of the center (in meters)
          x_size
            X size of the aperture (in meters)
          y_size
            Y size of the aperture (in meters)
:model: Transform from relative positions within the IFU slicer to absolute positions within the field of view. It's a combination of shifts and rotation.



MSA
:::
The MSA reference file contains information on the metrology of the microshutter array and the associated fixed slits - relative positioning of each individual shutter (assumed to be rectangular)
And the absolute position of each quadrant within the MSA.

The MSA reference file has 5 fields, named

:1:
   :data: Array with reference data for each shutter in Quadrant 1.
          It has 5 columns

          NO
            Shutter number (1- 62415)
          x_center
            X coordinate of the center (in meters)
          y_center
            Y coordinate of the center (in meters)
          x_size
            X size of the aperture (in meters)
          y_size
            Y size of the aperture (in meters)
   :model: Transform from relative positions within Quadrant 1 to absolute positions within the MSA
:2:
   :data: Array with reference data for shutters in Quadrant 2, same as in 1 above
   :model: Transform from relative positions within Quadrant 2 to absolute positions within the MSA
:3:
   :data: Array with reference data for shutters in Quadrant 3, same as in 1 above
   :model: Transform from relative positions within Quadrant 3 to absolute positions within the MSA
:4:
   :data: Array with reference data for shutters in Quadrant 4, same as in 1 above
   :model: Transform from relative positions within Quadrant 4 to absolute positions within the MSA
:5:
   :data: Reference data for the fixed slits and the IFU, same as in 1, except NO is 6 rows (1-6)
          and the mapping is 1 - S200A1, 2 - S200A1, 3 - S400A1, 4 - S200B1, 5 - S1600A1, 6 - IFU
   :model: Transform from relative positions within each aperture to absolute positions within the MSA


OTE
:::
This reference file contains a combination of astropy models - polynomial, shift, rotation and scaling.


:model: Transform through the Telescope Optical Element (OTE), from the FWA to XAN, YAN telescope frame. The output units are in arcsec.


SPECWCS
:::::::
For the MIRI LRS mode the file is in FITS format.
The reference file contains the zero point offset for the slit relative to the full field of view.
For the Fixed Slit exposure type the zero points in X and Y are stored in the header of the second HDU in the
'IMX' and 'IMY' keywords. For the Slitless exposure type they are stored in the header of the second HDU in
FITS keywords 'IMXSLTl' and 'IMYSLTl'. For both of the exposure types, the zero point offset is 1 based and the
X (e.g., IMX) refers to the column and Y refers to the row.

For the MIRI MRS the file is in ASDF format with the following structure.

:channel: The MIRI channels in the observation, e.g. "12".
:band: The band for the observation (one of "LONG", "MEDIUM", "SHORT").
:model:
        :slice_number: The wavelength solution for each slice.
                       <slice_number> is the actual slice number (s), computed by s = channel * 100 + slice

For NIRISS SOSS mode the file is in ASDF format with the following structure.

:model: A tabular model with the wavelength solution.



REGIONS
:::::::

The IFU takes a region reference file that defines the region over which the WCS is valid. The reference file should define a polygon and may consist of a set of X,Y coordinates that define the polygon.

:channel: The MIRI channels in the observation, e.g. "12".
:band: The band for the observation (one of "LONG", "MEDIUM", "SHORT").
:regions: An array with the size of the MIRI MRS image where pixel values map to the MRS slice number. 0 indicates a pixel is not within any slice.



V2V3
::::
The model field in the tree contains N models, one per channel, that map the spatial coordinates from alpha, beta to XAN, YAN.

:channel: The MIRI channels in the observation, e.g. "12".
:band: The band for the observation (one of "LONG", "MEDIUM", "SHORT").
:model:
        :channel_band: Transform from alpha, beta to XAN, YAN for this channel.



WAVELENGTHRANGE
:::::::::::::::

For MIRI MRS the wavelengthrange file consists of two fields which define the wavelength range for each combination of a channel and band.

:channels: An ordered list of all possible channel and band combinations for MIRI MRS, e.g. "1SHORT".
:wavelengthrange: An ordered list of (lambda_min, lambda_max) for each item in the list above.

For NIRSPEC the file is a dictionary storing information about default wavelength range and spectral order for each combination of filter and grating.

:filter_grating:
                 :order: Default spectral order
                 :range: Default wavelength range


