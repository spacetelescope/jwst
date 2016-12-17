Reference File Descriptions - Build 7
=====================================

WCS Reference files are in the Advanced Scientific Data Format (ASDF).
The best way to create the file is to programmatically create the model and then save it to a file.
A tutorial on creating reference files in ASDF format is available at:

https://github.com/spacetelescope/jwreftools/blob/master/docs/notebooks/referece_files_asdf.ipynb



List of reference types used by assign_wcs
------------------------------------------



===================    ==========================================================   ============================
reftype                                     description                              Instrument
===================    ==========================================================   ============================
**camera**             NIRSPEC Camera model                                          NIRSPEC
**collimator**         NIRSPEC Collimator Model                                      NIRSPEC
**disperser**          Disperser parameters                                          NIRSPEC
**distortion**         Spatial distortion model                                      MIRI, FGS, NIRCAM, NIRISS
**filteroffset**       MIRI Imager fiter offsets                                     MIRI
**fore**               Transform through the NIRSPEC FORE optics                     NIRSPEC
**fpa**                Transform in the NIRSPEC FPA plane                            NIRSPEC
**ifufore**            Transform from the IFU slicer to the IFU entrance             NIRSPEC
**ifupost**            Transform from the IFU slicer to the back of the IFU          NIRSPEC
**ifuslicer**          FU Slicer geometric description                               NIRSPEC
**msa**                Transformin the NIRSPEC MSA plane                             NIRSPEC
**ote**                Transform through the Optical Telescope Element               NIRSPEC
**specwcs**            Wavelength calibration models                                 MIRI, NIRCAM, NIRISS
**regions**            Stores location of the regions on the detector                MIRI
**v2v3**               Transform from MIRI instrument focal plane to V2V3 plane      MIRI
**wavelengthrange**    Typical wavelength ranges                                     MIRI, NIRSPEC
===================    ==========================================================   ============================



CAMERA
------

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The camera reference file contains a field “model” which stores a compound model made up of a polynomial model and then rotations and translations.  The model represents the transform through the camera. The forward direction is from the FPA to the GWA.

COLLIMATOR
----------

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The collimator reference file contains a field “model” which stores a compound model made up of a polynomial model and then rotations and translations.  The model represents the transform through the collimator. The forward direction is from the GWA to the MSA.

DISPERSER
---------

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE, GRATING


Reference File Formats
::::::::::::::::::::::

The disperser reference file contains reference data about the NIRSPEC dispersers (gratings or the prism). The reference data is described in the NIRSPEC Interface Control Document.

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

:angle: Angle between the front and back surface of the prosm (in degrees)
:kcoef: K coefficients of Selmeir equation, describing the material
:lcoef: L coeffficients describing the material
:tcoef: Thermal coefficients describing the properties of the glass
:tref: Reference temperature (in K)
:pref: Reference pressure (in ATM)
:wbound: Min and Max wavelength (in meters) for which the model is valid

DISTORTION
----------

CRDS Selection Criteria
:::::::::::::::::::::::

:MIRI: DETECTOR, EXP_TYPE. CHANNEL, BAND
:FGS: DETECTOR, EXP_TYPE
:NIRCAM: DETECTOR, EXP_TYPE. CHANNEL.
:NIRISS: EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The distortion reference file is an ASDF file that contains the distortion model. The distortion model contains one field indexed by the key "model" and returns a distortion model. The distortion model contains the forward and reverse transforms to/ from pixel to intermediate frame (different for different instruments).

FILTEROFFSET
------------

CRDS Selection Criteria
:::::::::::::::::::::::

MIRI: DETECTOR, EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The filter offset reference file is an ASDF file that contains a dictionary of row and column offsets for the MIRI imaging dataset. The filter offset reference file must contain a dictionary in the tree that is indexed by the instrument filter. The dictionary must contain two fields needed from the filter offset reference file: row_offset and column_offset and must be in units of arc-minutes.

FORE
----

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE, FILTER

Reference File Formats
::::::::::::::::::::::

The FORE reference file stores the transform through the Filter Wheel Assembly (FWA). It has two fields - “filter” and “model”. The transform through the FWA is chromatic. It is represented as a Polynomial of two variables whose coefficients are wavelength dependent. The compound model takes three inputs - x, y positions and wavelength.

FPA
---

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The FPA reference file stores information on the metrology of the Focal Plane Array (FPA) which consists of two single chip arrays (SCA), named NRS1 and NRS2.

The reference file contains two fields : “NRS1” and “NRS2”. Each of them stores the transform (shift and rotation) to transform positions from the FPA to the respective SCA. The output units are in pixels.

IFUFORE
-------

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE

Reference File Formats
::::::::::::::::::::::

This file provides the parameters (Paraxial and distortions coefficients)
for the coordinate transforms from the MSA plane to the plane of the IFU slicer.


IFUPOST
-------

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The IFUPOST reference file provides the parameters (Paraxial and distortions coefficients) for the coordinate transforms from the slicer plane to the MSA plane (out), that is the plane of the IFU virtual slits.

The reference file contains models made up based on an offset and a polynomial. There is a model for each of the slits and is indexed by the slit number. The models is used as part of the conversion from the GWA to slit.


IFUSLICER
---------

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE


Reference File Formats
::::::::::::::::::::::

The IFUSLICER stores information about the metrology of the IFU slicer - relative positioning and size of the aperture of each individual slicer and the absolute reference with respect to the center of the field of view.
The reference file contains two fields - “data” and “model”.
The “data” field is an array with 30 rows pertaining to the 30 slices and the columns are

slice number - [0 - 29]
x center - in meters
y center - in meters
x size - in meters
y size - in meters

The “model” field stores the model transforming positions from relative frame within the IFU slicer to the absolute position in the field of view. It’s a combination of shifts and rotation.

MSA
---

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE

Reference File Formats
::::::::::::::::::::::

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
            X size of teh aperture (in meters)
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
   :data: Reference data for the fixed slits and the IFU, same as in 1, exceppt NO for which the maping is
   
           1
             S200A1
           2
             S200A2
           3
             S400A1
           4
             S200B1
           5
             S1600A1
           6
             IFU
   
   :model: Transform from relative positions within eac aperture to absolute positions within the MSA


OTE
---

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The OTE reference file contains the transform through the Optical telescope Element (OTE).
It has one field - “model” which stores the transform from the FWA to XAN, YAN telescope frame.
The output units are in arcsec.

SPECWCS
-------

CRDS Selection Criteria
:::::::::::::::::::::::

MIRI: DETECTOR, CHANNEL, BAND, SUBARRAY, EXP_TYPE
NIRISS: EXP_TYPE, SUBARRAY

Reference File Formats
::::::::::::::::::::::

The reference file contains the zero point offset for the slit relative to the full field of view. For the Fixed Slit exposure type the fields are stored in the header of the second HDU and are indexed by 'imx' and 'imy'. For the Slitless exposure type the fields are stored in the header of the second HDU and are indexed by 'imxsltl' and 'imysltl'. For both of the exposure types, the zero point offset is 1 based and the X (e.g., imx) refers to the column and Y refers to the row.

Regions
-------

CRDS Selection Criteria
:::::::::::::::::::::::

MIRI: DETECTOR, CHANNEL, BAND, EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The IFU takes a region reference file that defines the region over which the WCS is valid. The reference file should define a polygon and may consist of a set of X,Y coordinates that define the polygon.

V2V3
----

CRDS Selection Criteria
:::::::::::::::::::::::

MIRI: DETECTOR, CHANNEL, BAND, EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The model field in the tree contains N models, one per channel, that map the spatial coordinates from alpha, beta to V2, V3.

WAVELENGTHRANGE
---------------

CRDS Selection Criteria
:::::::::::::::::::::::

NIRSPEC: Match EXP_TYPE
MIRI: Match EXP_TYPE

Reference File Formats
::::::::::::::::::::::

The wavelengthrange reference file consists of two models, one that defines the wavelength range and is indexed by 'wavelengthrange' and the second is a set of channels indexed in the file by 'channels'. The model defines, per channel, the wavelength mapping in going from alpha, beta to XAN, YAN.



Observing modes supported in build 7
------------------------------------

:FGS_IMAGE:

  | reftypes: *distortion*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference file provided by NIRISS team

:MIR_IMAGE:

  | reftypes: *distortion*, *filteroffset*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: CDP6 reference data delivery, MIRI-TN-00070-ATC_Imager_distortion_CDP_Iss5.pdf


:MIR_LRS-FIXEDSLIT, MIR_LRS-SLITLESS:

  | reftypes: *specwcs*, *distortion*
  | CRDS rmap rules: SUBARRAY.name: GENERIC
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: CDP6 reference data delivery, MIRI-TR-10020-MPI-Calibration-Data-Description_LRSPSFDistWave_v4.0.pdf


:MIR_MRS:

  | reftypes: *distortion*, *specwcs*, *v2v3*, *wavelengthrange*, *regions*
  | CRDS rmap rules: EXP_TYPE, DETECTOR, CHANNEL, BAND
  | WCS pipeline coordinate frames: detector, miri_focal, xyan, v2v3, world
  | Implements: CDP4 reference data delivery, MIRI-TN-00001-ETH_Iss1-3_Calibrationproduct_MRS_d2c.pdf

:NRC_IMAGE:

  | reftypes: *distortion*
  | CRDS rmap rules: EXP_TYPE, DETECTOR, CHANNEL, BAND
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: Distortion file created from TEL team data.

:NIS_IMAGE:

  | reftypes: *distortion*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference file provided by NIRISS team

:NIS_SOSS:

  | reftypes: *distortion*, *specwcs*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference files provided by NIRISS team

:NRS_FIXEDSLIT:
:NRS_MSASPEC:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 2 delivery

:NRS_IFU:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*,
  | *ifufore*, *ifuslicer*, *ifupost*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 2 delivery

:NRS_IMAGING:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 2 delivery

