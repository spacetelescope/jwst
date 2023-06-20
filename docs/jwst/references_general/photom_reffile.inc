.. _photom_reffile:

PHOTOM Reference File
---------------------

:REFTYPE: PHOTOM

The PHOTOM reference file contains conversion factors for putting
pixel values into physical units.

.. include:: ../references_general/photom_selection.inc

.. include:: ../includes/standard_keywords.inc

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
EXP_TYPE   model.meta.exposure.type        All
BAND       model.meta.instrument.band      MIRI
=========  ==============================  =========================

Tabular PHOTOM Reference File Format
++++++++++++++++++++++++++++++++++++
PHOTOM reference files are FITS format.
For all modes except MIRI MRS, the PHOTOM file contains tabular data
in a BINTABLE extension with EXTNAME = 'PHOTOM'.
The FITS primary HDU does not contain a data array.
The contents of the table extension vary a bit for different
instrument modes, as shown in the tables below.

:Data model: `~jwst.datamodels.FgsImgPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| FGS        | Image | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.MirImgPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| MIRI       | Image | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | subarray       | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.MirLrsPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| MIRI       | LRS   | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | subarray       | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | array      | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | array      | unitless               |
+            +       +----------------+-----------+------------+------------------------+
|            |       | reluncertainty | float     | array      | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NrcImgPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRCam     | Image | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | pupil          | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NrcWfssPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRCam     | WFSS  | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | pupil          | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | order          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | array      | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | array      | unitless               |
+            +       +----------------+-----------+------------+------------------------+
|            |       | reluncertainty | float     | array      | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NisImgPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRISS     | Image | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | pupil          | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NisSossPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRISS     | SOSS  | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | pupil          | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | order          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmj         | float     | scalar     | MJy/(DN/sec)           |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/(DN/sec)           |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | array      | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | array      | unitless               |
+            +       +----------------+-----------+------------+------------------------+
|            |       | reluncertainty | float     | array      | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NisWfssPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRISS     | WFSS  | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | pupil          | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | order          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmjsr       | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/steradian/(DN/sec) |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | array      | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | array      | unitless               |
+            +       +----------------+-----------+------------+------------------------+
|            |       | reluncertainty | float     | array      | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NrsFsPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRSpec    | Fixed | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            | Slit  | grating        | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | slit           | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | photmj         | float     | scalar     | MJy/(DN/sec)           |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/(DN/sec)           |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | array      | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | array      | unitless               |
+            +       +----------------+-----------+------------+------------------------+
|            |       | reluncertainty | float     | array      | unitless               |
+------------+-------+----------------+-----------+------------+------------------------+

:Data model: `~jwst.datamodels.NrsMosPhotomModel`

+------------+-------+----------------+-----------+------------+------------------------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units                  |
+============+=======+================+===========+============+========================+
| NIRSpec    | MOS   | filter         | string    | 12         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            | and   | grating        | string    | 15         | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            | IFU   | photmj         | float     | scalar     | MJy/(DN/sec)           |
+            +       +----------------+-----------+------------+------------------------+
|            |       | uncertainty    | float     | scalar     | MJy/(DN/sec)           |
+            +       +----------------+-----------+------------+------------------------+
|            |       | nelem          | integer   | scalar     | N/A                    |
+            +       +----------------+-----------+------------+------------------------+
|            |       | wavelength     | float     | array      | microns                |
+            +       +----------------+-----------+------------+------------------------+
|            |       | relresponse    | float     | array      | unitless               |
+            +       +----------------+-----------+------------+------------------------+
|            |       | reluncertainty | float     | array      | unitless               |
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

Note: For spectroscopic data the Nelem column should be present. Its value must be
greater than 0, and Nelem entries are read from each of the Wavelength and Relresponse
arrays. Nelem is not used for imaging data because there are no columns containing arrays.

The primary header of the tabular PHOTOM reference files contains the keywords PIXAR_SR
and PIXAR_A2, which give the average pixel area in units of steradians and
square arcseconds, respectively.

MIRI MRS Photom Reference File Format
+++++++++++++++++++++++++++++++++++++

:Data model: `~jwst.datamodels.MirMrsPhotomModel`

For MIRI MRS, the PHOTOM file contains 2-D arrays of conversion factors in
IMAGE extensions.
The FITS primary HDU does not contain a data array.
The format and content of the MIRI MRS PHOTOM reference file is as follows:

==================  ========  =====  ===========  =========  
EXTNAME             XTENSION  NAXIS  Dimensions   Data type
==================  ========  =====  ===========  =========
SCI                 IMAGE       2    1032 x 1024  float
ERR                 IMAGE       2    1032 x 1024  float
DQ                  IMAGE       2    1032 x 1024  integer
PIXSIZ              IMAGE       2    1032 x 1024  float
DQ_DEF              BINTABLE    2    TFIELDS = 4  N/A
TIMECOEFF_CH1 [1]_  BINTABLE    2    TFIELDS = 5  N/A
TIMECOEFF_CH2 [1]_  BINTABLE    2    TFIELDS = 5  N/A
TIMECOEFF_CH3 [1]_  BINTABLE    2    TFIELDS = 5  N/A
TIMECOEFF_CH4 [1]_  BINTABLE    2    TFIELDS = 5  N/A
==================  ========  =====  ===========  =========

.. rubric:: Footnotes

.. [1] Optional extension. If present, the MRS time-dependent throughput correction can be applied.

The SCI extension contains a 2D array of inverse sensitivity factors
corresponding to each pixel in a 2D MRS slice image. The sensitivity factors
are in units of (MJy/pixel)/(DN/sec). The ERR extension contains a 2D array of
uncertainties for the SCI values, in the same units. The DQ extension
contains a 2D array of bit-encoded data quality flags for the SCI values.
The DQ_DEF extension contains a table listing the definitions of the values
used in the DQ array. The PIXSIZ extension contains a 2D array of pixel
sizes (i.e. solid angles), in units of square-arcsec.

The SCI and PIXSIZ array values are both divided into the science product
SCI and ERR arrays, yielding surface brightness in units of mJy/sq-arcsec.

Scalar PHOTMJSR and PHOTUJA2 values are stored in primary header keywords
in the MIRI MRS PHOTOM reference files and are copied into the science
product header by the photom step.

The TIMECOEFF_CH tables contain the parameters to correct the MRS time-dependent
throughput loss. If these tables do not exist in the reference file, then the MIRI
MRS time-dependent correction is skipped. 

Constructing a PHOTOM Reference File
------------------------------------
The most straight-forward way to construct a tabular PHOTOM reference file is to
populate a data model within python and then save the data model to a
FITS file. Each instrument mode has its own photom data model, as listed above,
which contains the columns of information unique to that instrument.

A NIRCam WFSS photom reference file, for example, could be constructed as follows
from within the python environment::

 >>> import numpy as np
 >>> from stdatamodels.jwst import datamodels
 >>> filter = np.array(['GR150C', 'GR150R'])
 >>> pupil = np.array(['F140M', 'F200W'])
 >>> order = np.array([1, 1], dtype=np.int16)
 >>> photf = np.array([1.e-15, 3.e-15], dtype=np.float32)
 >>> uncer = np.array([1.e-17, 3.e-17], dtype=np.float32)
 >>> nrows = len(filter)
 >>> nx = 437
 >>> nelem = np.zeros(nrows, dtype=np.int16) + nx
 >>> temp_wl = np.linspace(1.0, 5.0, nx, dtype=np.float32).reshape(1, nx)
 >>> wave = np.zeros((nrows, nx), np.float32)
 >>> wave[:] = temp_wl.copy()
 >>> resp = np.ones((nrows, nx), dtype=np.float32)
 >>> resp_unc = np.zeros((nrows, nx), dtype=np.float32)
 >>> data_list = [(filter[i], pupil[i], order[i], photf[i], uncer[i], nelem[i],
 ...               wave[i], resp[i], resp_unc[i]) for i in range(nrows)]
 >>> data = np.array(data_list,
 ...                 dtype=[('filter', 'S12'),
 ...                        ('pupil', 'S15'),
 ...                        ('order', '<i2'),
 ...                        ('photmjsr', '<f4'),
 ...                        ('uncertainty', '<f4'),
 ...                        ('nelem', '<i2'),
 ...                        ('wavelength', '<f4', (nx,)),
 ...                        ('relresponse', '<f4', (nx,)),
 ...                        ('reluncertainty', '<f4', (nx,))])
 >>> output = datamodels.NrcWfssPhotomModel(phot_table=data)
 >>> output.save('nircam_photom_0001.fits')
 'nircam_photom_0001.fits'
