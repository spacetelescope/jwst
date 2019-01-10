.. _photom_reffile:

PHOTOM Reference File
---------------------

:REFTYPE: PHOTOM

The PHOTOM reference file contains conversion factors for putting
pixel values into physical units.

.. include:: ../references_general/photom_selection.rst

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
are in units of (DN/sec)/(mJy/pixel). The ERR extension contains a 2D array of
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

Constructing a PHOTOM Reference File
------------------------------------
The most straight-forward way to construct a tabular PHOTOM reference file is to
populate a data model within python and then save the data model to a
FITS file. Each instrument has its own photom data model, as listed above,
which contains the columns of information unique to that instrument.

A NIRISS photom reference file, for example, could be constructed as follows
from within the python environment::

 >>> from jwst import datamodels
 >>> import numpy as np
 >>> output = datamodels.NirissPhotomModel()
 >>> filter = np.array(['F277W','F356W','CLEAR'])
 >>> pupil = np.array(['CLEARP','CLEARP','F090W'])
 >>> photf = np.array([1.e-15,2.e-15,3.e-15])
 >>> uncer = np.array([1.e-17,2.e-17,3.e-17])
 >>> nelem = np.zeros(3, dtype=np.int16)
 >>> wave = np.zeros((3, 5000))
 >>> resp = np.zeros((3, 5000))
 >>> order = np.array([1, 2, 3], dtype=np.int16)
 >>> data = np.array(list(zip(filter, pupil, order, photf, uncer, nelem, wave, resp)), dtype=output.phot_table.dtype)
 >>> output.phot_table = data
 >>> output.save('niriss_photom_0001.fits')
 'niriss_photom_0001.fits'

