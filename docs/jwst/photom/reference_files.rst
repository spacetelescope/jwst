Reference Files
===============
The ``photom`` step uses PHOTOM_ and pixel AREA_ reference files.
The AREA reference file is only used when processing
imaging and NIRSpec IFU observations.

.. _PHOTOM:

PHOTOM Reference File
---------------------

:REFTYPE: PHOTOM

The PHOTOM reference file contains conversion factors for putting
pixel values into physical units.

.. include:: photom_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for PHOTOM
+++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in PHOTOM reference files,
because they are used as CRDS selectors
(see :ref:`photom_selectors`):

=========  ==============================  =========================
Keyword    Data Model Name                 Instruments
=========  ==============================  =========================
DETECTOR   model.meta.instrument.detector  FGS, MIRI, NIRCam, NIRISS
BAND       model.meta.instrument.band      MIRI
EXP_TYPE   model.meta.exposure.type        NIRSpec
=========  ==============================  =========================

Tabular PHOTOM Reference File Format
++++++++++++++++++++++++++++++++++++
PHOTOM reference files are FITS format.
For all modes except MIRI MRS, the PHOTOM file contains tabular data
in a BINTABLE extension with EXTNAME = 'PHOTOM'.
The FITS primary HDU does not contain a data array.
The contents of the table extension vary a bit for different
instrument modes, as shown in the tables below.

:Data model: `~jwst.datamodels.FgsPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| FGS        | Image | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | 5000       | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | 5000       | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.MiriImgPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| MIRI       | Image | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            | and   | subarray       | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            | LRS   | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | 500        | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | 500        | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NircamPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRCam     | All   | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | pupil          | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | order          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | 3000       | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | 3000       | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NirissPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRISS     | All   | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | pupil          | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | order          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | 5000       | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | 5000       | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NirspecPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRSpec    | IFU   | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            | and   | grating        | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            | MOS   | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | 150        | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | 150        | unitless               |
+            +       +----------------+-----------+------------+------------------------+
|            |       | reluncertainty | float     | 150        | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NirspecFSPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRSpec    | Fixed | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            | Slit  | grating        | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | slit           | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | 150        | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | 150        | unitless               |
+            +       +----------------+-----------+------------+------------------------+
|            |       | reluncertainty | float     | 150        | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

Row Selection
^^^^^^^^^^^^^
A row of data within the table is selected by the ``photom`` step based on
the optical elements in use for the exposure. The selection attributes are
always contained in the first few columns of the table. The remaining
columns contain the data needed for photometric conversion.
The row selection criteria for each instrument/mode are:

* FGS:
   - No selection criteria (table contains a single row)
* MIRI:
   - Imager and LRS: Filter and Subarray
   - MRS: Does not use table-based reference file (see below)
* NIRCam:
   - All: Filter and Pupil
* NIRISS:
   - Imaging: Filter and Pupil
   - Spectroscopic: Filter, Pupil, and Order number
* NIRSpec:
   - IFU and MOS: Filter and Grating
   - Fixed Slits: Filter, Grating, and Slit name

Note: If Nelem = 0 in a given row, then no data are read from the Wavelength
and Relresponse columns. If Nelem > 0, then Nelem entries are read from each
of the Wavelength and Relresponse arrays.

The primary header of the tabular PHOTOM reference files contains the keywords PIXAR_SR
and PIXAR_A2, which give the average pixel area in units of steradians and
square arcseconds, respectively.

MIRI MRS Photom Reference File Format
+++++++++++++++++++++++++++++++++++++

:Data model: `~jwst.datamodels.MiriMrsPhotomModel`

For MIRI MRS, the PHOTOM file contains 2-D arrays of conversion factors in
IMAGE extensions.
The FITS primary HDU does not contain a data array.
The format and content of the MIRI MRS PHOTOM reference file is as follows:

=======  ========  =====  ===========  =========
EXTNAME  XTENSION  NAXIS  Dimensions   Data type
=======  ========  =====  ===========  =========
SCI      IMAGE       2    1032 x 1024  float
ERR      IMAGE       2    1032 x 1024  float
DQ       IMAGE       2    1032 x 1024  integer
PIXSIZ   IMAGE       2    1032 x 1024  float
DQ_DEF   BINTABLE    2    TFIELDS = 4  N/A
=======  ========  =====  ===========  =========

The SCI extension contains a 2D array of spectral sensitivity factors
corresponding to each pixel in a 2D MRS slice image. The sensitivity factors
are in units of DN/sec/mJy/pixel. The ERR extension contains a 2D array of
uncertainties for the SCI values, in the same units. The DQ extension
contains a 2D array of bit-encoded data quality flags for the SCI values.
The DQ_DEF extension contains a table listing the definitions of the values
used in the DQ array. The PIXSIZ extension contains a 2D array of pixel
sizes, in units of square-arcsec.

The SCI and PIXSIZ array values are both divided into the science product
SCI and ERR arrays, yielding surface brightness in units of mJy/sq-arcsec.

Scalar PHOTMJSR and PHOTUJA2 values are stored in primary header keywords
in the MIRI MRS PHOTOM reference files and are copied into the science
product header by the photom step.

.. _AREA:

AREA Reference File
-------------------

:REFTYPE: AREA

The AREA reference file contains pixel area information for a given
instrument mode.

.. include:: area_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for AREA
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in AREA reference files,
because they are used as CRDS selectors
(see :ref:`area_selectors`):

=========  ==============================  =============================
Keyword    Data Model Name                 Instrument
=========  ==============================  =============================
DETECTOR   model.meta.instrument.detector  All
EXP_TYPE   model.meta.exposure.type        MIRI, NIRCam, NIRISS, NIRSpec
FILTER     model.meta.instrument.filter    MIRI, NIRCam, NIRISS, NIRSpec
PUPIL      model.meta.instrument.pupil     NIRCam, NIRISS
GRATING    model.meta.instrument.grating   NIRSpec
=========  ==============================  =============================

Reference File Format
+++++++++++++++++++++
AREA reference files are FITS format.
For imaging modes (FGS, MIRI, NIRCam, and NIRISS) the AREA reference files
contain 1 IMAGE extension, while reference files for NIRSpec spectroscopic
modes contain 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.

Imaging Modes
^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.PixelAreaModel`

The format of imaging mode AREA reference files is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
=======  ========  =====  ==============  =========

The SCI extension data array contains a 2-D pixel-by-pixel map of relative
pixel areas, normalized to a value of 1.0. The absolute value of the nominal
pixel area is given in the primary header keywords PIXAR_SR and PIXAR_A2, in
units of steradians and square arcseconds, respectively.
These keywords should have the same values as the
PIXAR_SR and PIXAR_A2 keywords in the header of the corresponding PHOTOM
reference file.

NIRSpec Fixed-Slit Mode
^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.NirspecSlitAreaModel`

The BINTABLE extension has EXTNAME='AREA' and has column characteristics
shown below. There is one row for each of the 5 fixed slits, with ``slit_id``
values of "S200A1", "S200A2", "S400A1", "S200B1", and "S1600A1". The pixel
area values are in units of square arcseconds.

===========  =========
Column name  Data type 
===========  =========
slit_id      char\*15
pixarea      float
===========  =========

NIRSpec MOS Mode
^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.NirspecMosAreaModel`

The BINTABLE extension has EXTNAME='AREA' and has column characteristics
shown below. There is one row for each shutter in each MSA quadrant. The
quadrant and shutter values are 1-indexed. The pixel area values are in
units of square arcseconds.

===========  =========
Column name  Data type 
===========  =========
quadrant     integer
shutter_x    integer
shutter_y    integer
pixarea      float
===========  =========

NIRSpec IFU Mode
^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.NirspecIfuAreaModel`

The BINTABLE extension has EXTNAME='AREA' and has column characteristics
shown below. There is one row for each of the 30 IFU slices, with the
``slice_id`` values being 0-indexed (i.e. range from 0 to 29). The pixel
area values are in units of square arcseconds.

===========  =========
Column name  Data type 
===========  =========
slice_id     integer
pixarea      float
===========  =========


Constructing a PHOTOM Reference File
------------------------------------
The most straight-forward way to construct a tabular PHOTOM reference file is to
populate a data model within python and then save the data model to a
FITS file. Each instrument has its own photom data model, as listed above,
which contains the columns of information unique to that instrument.

A NIRISS photom reference file, for example, could be constructed as follows
from within the python environment::

 >>> from jwst import models
 >>> import numpy as np
 >>> output=models.NirissPhotomModel()
 >>> filter=np.array(['F277W','F356W','CLEAR'])
 >>> pupil=np.array(['CLEARP','CLEARP','F090W'])
 >>> photf=np.array([1.e-15,2.e-15,3.e-15])
 >>> uncer=np.array([1.e-17,2.e-17,3.e-17])
 >>> nelem=np.zeros(3)
 >>> wave=np.zeros(3)
 >>> resp=np.zeros(3)
 >>> data=np.array(list(zip(filter,pupil,photf,uncer,nelem,wave,resp)),dtype=output.phot_table.dtype)
 >>> output.phot_table=data
 >>> output.save('niriss_photom_0001.fits')

