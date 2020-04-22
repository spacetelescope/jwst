APCORR Reference File
---------------------

:REFTYPE: APCORR

The APCORR reference file contains data necessary for correcting
extracted imaging and spectroscopic photometry to the equivalent
of an infinite aperture. It is used within the
:ref:`source_catalog <source_catalog_step>` step for imaging and
within the :ref:`extract_1d <extract_1d_step>` step for
spectroscopic data.

.. include:: ../references_general/apcorr_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for APCORR
+++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in APCORR reference files,
because they are used as CRDS selectors
(see :ref:`apcorr_selectors`):

=========  ========================  ===========
Keyword    Data Model Name           Instruments
=========  ========================  ===========
EXP_TYPE   model.meta.exposure.type  All
=========  ========================  ===========

APCORR Reference File Format
++++++++++++++++++++++++++++
APCORR reference files are FITS format.
The APCORR reference file contains tabular data
in a BINTABLE extension with EXTNAME = 'APCORR'.
The FITS primary HDU does not contain a data array.
The contents of the table extension varies for different
instrument modes, as shown in the tables below.

:Data model: `~jwst.datamodels.FgsImgApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| FGS        | Image | eefraction  | float     | scalar     | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | radius      | float     | scalar     | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | scalar     | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | skyin       | float     | scalar     | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | skyout      | float     | scalar     | pixels   |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.MirImgApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| MIRI       | Image | filter      | string    | 12         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | subarray    | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | eefraction  | float     | scalar     | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | radius      | float     | scalar     | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | scalar     | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | skyin       | float     | scalar     | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | skyout      | float     | scalar     | pixels   |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.MirLrsApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| MIRI       | LRS   | subarray    | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | wavelength  | float     | 1D array   | micron   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_wl    | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | size        | integer   | 1D array   | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_size  | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | 2D array   | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr_err  | float     | 2D array   | unitless |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.MirMrsApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| MIRI       | MRS   | wavelength  | float     | 1D array   | micron   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_wl    | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | radius      | float     | 1D array   | arcsec   |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | 1D array   | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr_err  | float     | 1D array   | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | inner_bkg   | float     | 1D array   | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | outer_bkg   | float     | 1D array   | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | axis_ratio  | float     | 1D array   | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | axis_pa     | float     | 1D array   | N/A      |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.NrcImgApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| NIRCam     | Image | filter      | string    | 12         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | pupil       | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | eefraction  | float     | scalar     | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | radius      | float     | scalar     | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | scalar     | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | skyin       | float     | scalar     | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | skyout      | float     | scalar     | pixels   |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.NrcWfssApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| NIRCam     | WFSS  | filter      | string    | 12         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | pupil       | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | wavelength  | float     | 1D array   | micron   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_wl    | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | size        | integer   | 1D array   | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_size  | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | 2D array   | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr_err  | float     | 2D array   | unitless |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.NisImgApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| NIRISS     | Image | filter      | string    | 12         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | pupil       | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | eefraction  | float     | scalar     | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | radius      | float     | scalar     | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | scalar     | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | skyin       | float     | scalar     | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | skyout      | float     | scalar     | pixels   |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.NisWfssApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| NIRISS     | WFSS  | filter      | string    | 12         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | pupil       | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | wavelength  | float     | 1D array   | micron   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_wl    | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | size        | integer   | 1D array   | pixels   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_size  | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | 2D array   | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr_err  | float     | 2D array   | unitless |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.NrsFsApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| NIRSpec    |  FS   | filter      | string    | 12         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | grating     | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | slit        | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | wavelength  | float     | 1D array   | micron   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_wl    | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | size        | integer   | 2D array   | arcsec   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_size  | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | pixphase    | float     | 1D array   | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | 3D array   | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr_err  | float     | 3D array   | unitless |
+------------+-------+-------------+-----------+------------+----------+

:Data model: `~jwst.datamodels.NrsMosApcorrModel`

+------------+-------+-------------+-----------+------------+----------+
| Instrument | Mode  | Column name | Data type | Dimensions | Units    |
+============+=======+=============+===========+============+==========+
| NIRSpec    | MOS   | filter      | string    | 12         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            | and   | grating     | string    | 15         | N/A      |
+            +       +-------------+-----------+------------+----------+
|            | IFU   | wavelength  | float     | 1D array   | micron   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_wl    | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | size        | integer   | 2D array   | arcsec   |
+            +       +-------------+-----------+------------+----------+
|            |       | nelem_size  | integer   | scalar     | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | pixphase    | float     | 1D array   | N/A      |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr      | float     | 3D array   | unitless |
+            +       +-------------+-----------+------------+----------+
|            |       | apcorr_err  | float     | 3D array   | unitless |
+------------+-------+-------------+-----------+------------+----------+

Row Selection
^^^^^^^^^^^^^
A row of data within the reference table is selected by the pipeline step
based on the optical elements in use for the exposure. The selection
attributes are always contained in the first few columns of the table.
The remaining columns contain the data needed for the aperture correction.
The row selection criteria for each instrument+mode are:

* FGS Image:
   - None (table contains a single row)
* MIRI:
   - Image: Filter and Subarray
   - LRS: Subarray
   - MRS: None (table contains a single row)
* NIRCam:
   - Image: Filter and Pupil
   - WFSS: Filter and Pupil
* NIRISS:
   - Image: Filter and Pupil
   - WFSS: Filter and Pupil
* NIRSpec:
   - IFU and MOS: Filter and Grating
   - Fixed Slits: Filter, Grating, and Slit name

Note: Spectroscopic mode reference files contain the "nelem_wl" and "nelem_size"
entries, which indicate to the pipeline step how many valid elements are contained
in the "wavelength" and "size" arrays, respectively. Only the first "nelem_wl"
and "nelem_size" entries are read from each array.
