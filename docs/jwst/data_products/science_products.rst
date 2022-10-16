Science products
----------------
The following sections describe the format and contents of each of the JWST FITS science
products. Things to note in the descriptions include:

 - Not all FITS extensions occur in every data product of a given type. Many are either
   optional or dependent on the instrument or observing mode. Such optional extensions are
   noted with an asterisk in the tables below.

 - Because some extensions are optional, as well as the fact that the exact ordering of the
   extensions is not guaranteed, the FITS HDU index numbers of a given extension type can
   vary from one product to another. The only guarantee is that the ``SCI`` extension,
   containing the actual pixel values, will always be the first FITS extension (HDU=1).
   Other common extensions, like ``DQ`` and ``ERR``, usually immediately follow the ``SCI``,
   but the order is not guaranteed. Hence HDU index numbers are not listed for many
   extension types, because they can vary.

.. _uncal:

Uncalibrated raw data: ``uncal``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Exposure raw data products are designated by a file name
suffix of "uncal." These files usually contain only the raw detector pixel values
from an exposure, with the addition of some table extensions containing various types of
meta data associated with the exposure.
Additional extensions can be included for certain instruments and readout types, as noted
below.
The FITS file structure is as follows.

+-----+------------+----------+-----------+-----------------------------------+
| HDU | EXTNAME    | HDU Type | Data Type | Dimensions                        |
+=====+============+==========+===========+===================================+
|  0  | N/A        | primary  | N/A       | N/A                               |
+-----+------------+----------+-----------+-----------------------------------+
|  1  | SCI        | IMAGE    | uint16    | ncols x nrows x ngroups x nints   |
+-----+------------+----------+-----------+-----------------------------------+
|  2  | GROUP      | BINTABLE | N/A       | variable                          |
+-----+------------+----------+-----------+-----------------------------------+
|  3  | INT_TIMES  | BINTABLE | N/A       | nints (rows) x 7 cols             |
+-----+------------+----------+-----------+-----------------------------------+
|     | ZEROFRAME* | IMAGE    | uint16    | ncols x nrows x nints             |
+-----+------------+----------+-----------+-----------------------------------+
|     | REFOUT*    | IMAGE    | uint16    | ncols/4 x nrows x ngroups x nints |
+-----+------------+----------+-----------+-----------------------------------+
|     | ASDF       | BINTABLE | N/A       | variable                          |
+-----+------------+----------+-----------+-----------------------------------+

 - SCI: 4-D data array containing the raw pixel values. The first two dimensions are equal to
   the size of the detector readout, with the data from multiple groups (NGROUPS) within each
   integration stored along the 3rd axis, and the multiple integrations (NINTS) stored along
   the 4th axis.
 - GROUP: A table of meta data for some (or all) of the data groups.
 - INT_TIMES: A table of beginning, middle, and end time stamps for each integration in the
   exposure.
 - ZEROFRAME: 3-D data array containing the pixel values of the zero-frame for each
   integration in the exposure, where each plane of the cube corresponds to a given integration.
   Only appears if the zero-frame data were requested to be downlinked separately.
 - REFOUT: The MIRI detector reference output values. Only appears in MIRI exposures.
 - ADSF: The data model meta data.

This FITS file structure is the result of serializing a `~jwst.datamodels.Level1bModel`, but
can also be read into a `~jwst.datamodels.RampModel`, in which case zero-filled
ERR, GROUPDQ, and PIXELDQ data arrays will be created and stored in the model, having array
dimensions based on the shape of the SCI array (see `~jwst.datamodels.RampModel`).

.. _ramp:

Ramp data: ``ramp``
^^^^^^^^^^^^^^^^^^^
As raw data progress through the :ref:`calwebb_detector1 <calwebb_detector1>` pipeline
they are stored internally in a `~jwst.datamodels.RampModel`.
This type of data model is serialized to a ``ramp`` type FITS
file on disk. The original detector pixel values (in the SCI extension) are converted
from integer to floating-point data type. The same is true for the ZEROFRAME and REFOUT
data extensions, if they are present. An ERR array and two types of data quality arrays are
also added to the product. The FITS file layout is as follows:

+-----+------------+----------+-----------+-----------------------------------+
| HDU | EXTNAME    | HDU Type | Data Type | Dimensions                        |
+=====+============+==========+===========+===================================+
|  0  | N/A        | primary  | N/A       | N/A                               |
+-----+------------+----------+-----------+-----------------------------------+
|  1  | SCI        | IMAGE    | float32   | ncols x nrows x ngroups x nints   |
+-----+------------+----------+-----------+-----------------------------------+
|  2  | PIXELDQ    | IMAGE    | uint32    | ncols x nrows                     |
+-----+------------+----------+-----------+-----------------------------------+
|  3  | GROUPDQ    | IMAGE    | uint8     | ncols x nrows x ngroups x nints   |
+-----+------------+----------+-----------+-----------------------------------+
|  4  | ERR        | IMAGE    | float32   | ncols x nrows x ngroups x nints   |
+-----+------------+----------+-----------+-----------------------------------+
|     | ZEROFRAME* | IMAGE    | float32   | ncols x nrows x nints             |
+-----+------------+----------+-----------+-----------------------------------+
|     | GROUP      | BINTABLE | N/A       | variable                          |
+-----+------------+----------+-----------+-----------------------------------+
|     | INT_TIMES  | BINTABLE | N/A       | nints (rows) x 7 cols             |
+-----+------------+----------+-----------+-----------------------------------+
|     | REFOUT*    | IMAGE    | uint16    | ncols/4 x nrows x ngroups x nints |
+-----+------------+----------+-----------+-----------------------------------+
|     | ASDF       | BINTABLE | N/A       | variable                          |
+-----+------------+----------+-----------+-----------------------------------+

 - SCI: 4-D data array containing the pixel values. The first two dimensions are equal to
   the size of the detector readout, with the data from multiple groups (NGROUPS) within each
   integration stored along the 3rd axis, and the multiple integrations (NINTS) stored along
   the 4th axis.
 - PIXELDQ: 2-D data array containing DQ flags that apply to all groups and all integrations
   for a given pixel (e.g. a hot pixel is hot in all groups and integrations).
 - GROUPDQ: 4-D data array containing DQ flags that pertain to individual groups within individual
   integrations, such as the point at which a pixel becomes saturated within a given integration.
 - ERR: 4-D data array containing uncertainty estimates on a per-group and per-integration basis.
 - ZEROFRAME: 3-D data array containing the pixel values of the zero-frame for each
   integration in the exposure, where each plane of the cube corresponds to a given integration.
   Only appears if the zero-frame data were requested to be downlinked separately.
 - GROUP: A table of meta data for some (or all) of the data groups.
 - INT_TIMES: A table of beginning, middle, and end time stamps for each integration in the
   exposure.
 - REFOUT: The MIRI detector reference output values. Only appears in MIRI exposures.
 - ADSF: The data model meta data.
 
.. _rate:
.. _rateints:

Countrate data: ``rate`` and ``rateints``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Countrate products are produced by applying the :ref:`ramp_fitting <ramp_fitting_step>` step to
the integrations within an exposure, in order to compute count rates from the original
accumulating signal ramps. For exposures that contain multiple integrations (NINTS > 1) this
is done in two ways, which results in two separate products. First, countrates are computed
for each integration within the exposure, the results of which are stored in a ``rateints`` product.
These products contain 3-D data arrays, where each plane of the data cube contains the
countrate image for a given integration.

The results for each integration are also averaged together to form a single 2-D countrate
image for the entire exposure. These results are stored in a ``rate`` product.

The FITS file structure for a ``rateints`` product is as follows:

+-----+-------------+----------+-----------+-----------------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions            |
+=====+=============+==========+===========+=======================+
|  0  | N/A         | primary  | N/A       | N/A                   |
+-----+-------------+----------+-----------+-----------------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  2  | ERR         | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  3  | DQ          | IMAGE    | uint32    | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  4  | INT_TIMES   | BINTABLE | N/A       | nints (rows) x 7 cols |
+-----+-------------+----------+-----------+-----------------------+
|  5  | VAR_POISSON | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  6  | VAR_RNOISE  | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  7  | ASDF        | BINTABLE | N/A       | variable              |
+-----+-------------+----------+-----------+-----------------------+

 - SCI: 3-D data array containing the pixel values, in units of DN/s. The first two dimensions are equal to
   the size of the detector readout, with the data from multiple integrations stored along the 3rd axis.
 - ERR: 3-D data array containing uncertainty estimates on a per-integration basis. These values
   are based on the combined VAR_POISSON and VAR_RNOISE data (see below), given as
   standard deviation.
 - DQ: 3-D data array containing DQ flags. Each plane of the cube corresponds to a given integration.
 - INT_TIMES: A table of beginning, middle, and end time stamps for each integration in the
   exposure.
 - VAR_POISSON: 3-D data array containing the per-integration variance estimates for each pixel,
   based on Poisson noise only.
 - VAR_RNOISE: 3-D data array containing the per-integration variance estimates for each pixel,
   based on read noise only.
 - ADSF: The data model meta data.

These FITS files are compatible with the `~jwst.datamodels.CubeModel` data model.

The FITS file structure for a ``rate`` product is as follows:

+-----+-------------+----------+-----------+-----------------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions            |
+=====+=============+==========+===========+=======================+
|  0  | N/A         | primary  | N/A       | N/A                   |
+-----+-------------+----------+-----------+-----------------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows         |
+-----+-------------+----------+-----------+-----------------------+
|  2  | ERR         | IMAGE    | float32   | ncols x nrows         |
+-----+-------------+----------+-----------+-----------------------+
|  3  | DQ          | IMAGE    | uint32    | ncols x nrows         |
+-----+-------------+----------+-----------+-----------------------+
|  4  | VAR_POISSON | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  5  | VAR_RNOISE  | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  6  | ASDF        | BINTABLE | N/A       | variable              |
+-----+-------------+----------+-----------+-----------------------+

 - SCI: 2-D data array containing the pixel values, in units of DN/s.
 - ERR: 2-D data array containing uncertainty estimates for each pixel. These values
   are based on the combined VAR_POISSON and VAR_RNOISE data (see below), given as
   standard deviation.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - VAR_POISSON: 2-D data array containing the variance estimate for each pixel,
   based on Poisson noise only.
 - VAR_RNOISE: 2-D data array containing the variance estimate for each pixel,
   based on read noise only.
 - ADSF: The data model meta data.

These FITS files are compatible with the `~jwst.datamodels.ImageModel` data model.

Note that the ``INT_TIMES`` table does not appear in ``rate`` products, because the
data have been averaged over all integrations and hence the per-integration time stamps
are no longer relevant.

.. _bsub:
.. _bsubints:

Background-subtracted data: ``bsub`` and ``bsubints``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`calwebb_image2 <calwebb_image2>` and :ref:`calwebb_spec2 <calwebb_spec2>`
pipelines have the capability to perform background subtraction on countrate data.
In its simplest form, this consists of subtracting background exposures or a
CRDS background reference image from science images. This operation is performed by
the :ref:`background <background_step>` step in the stage 2 pipelines. If the pipeline
parameter ``save_bsub`` is set to ``True``, the result of the background subtraction
step will be saved to a file. Because this is a direct image-from-image operation, the
form of the result is identical to input. If the input is a ``rate`` product, the
background-subtracted result will be a ``bsub`` product, which has the exact same
structure as the rate_ product described above. Similarly, if the input is a ``rateints``
product, the background-subtracted result will be saved to a ``bsubints`` product, with
the exact same structure as the rateints_ product described above.

.. _cal:
.. _calints:

Calibrated data: ``cal`` and ``calints``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Single exposure calibrated products duplicate a lot of the format and content of
countrate products. There are two different high-level forms of calibrated products:
one containing results for all integrations in an exposure (``calints``) and one for
results averaged over all integrations (``cal``). These products are the main result of
Stage 2 pipelines like :ref:`calwebb_image2 <calwebb_image2>` and
:ref:`calwebb_spec2 <calwebb_spec2>`. There are many additional types of extensions
that only appear for certain observing modes or instruments, especially for spectroscopic
exposures.

The FITS file structure for a ``calints`` product is as follows:

+-----+-------------+----------+-----------+-----------------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions            |
+=====+=============+==========+===========+=======================+
|  0  | N/A         | primary  | N/A       | N/A                   |
+-----+-------------+----------+-----------+-----------------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  2  | ERR         | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  3  | DQ          | IMAGE    | uint32    | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|     | INT_TIMES   | BINTABLE | N/A       | nints (rows) x 7 cols |
+-----+-------------+----------+-----------+-----------------------+
|     | VAR_POISSON | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|     | VAR_RNOISE  | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|     | VAR_FLAT    | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|     | AREA*       | IMAGE    |           | ncols x nrows         |
+-----+-------------+----------+-----------+-----------------------+
|     | WAVELENGTH* | IMAGE    | float32   | ncols x nrows         |
+-----+-------------+----------+-----------+-----------------------+
|     | ASDF        | BINTABLE | N/A       | variable              |
+-----+-------------+----------+-----------+-----------------------+

 - SCI: 3-D data array containing the pixel values, in units of surface brightness, for
   each integration.
 - ERR: 3-D data array containing uncertainty estimates for each pixel, for each integration.
   These values are based on the combined VAR_POISSON and VAR_RNOISE data (see below),
   given as standard deviation.
 - DQ: 3-D data array containing DQ flags for each pixel, for each integration.
 - INT_TIMES: A table of beginning, middle, and end time stamps for each integration in the
   exposure.
 - VAR_POISSON: 3-D data array containing the variance estimate for each pixel,
   based on Poisson noise only, for each integration.
 - VAR_RNOISE: 3-D data array containing the variance estimate for each pixel,
   based on read noise only, for each integration.
 - VAR_FLAT: 2-D data array containing the variance estimate for each pixel,
   based on uncertainty in the flat-field.
 - AREA: 2-D data array containing pixel area values, added by the :ref:`photom <photom_step>`
   step, for imaging modes.
 - WAVELENGTH: 2-D data array of wavelength values for each pixel, for some spectroscopic modes.
 - ADSF: The data model meta data.

The FITS file structure for a ``cal`` product is as follows:

+-----+---------------------------+----------+-----------+---------------+
| HDU | EXTNAME                   | HDU Type | Data Type | Dimensions    |
+=====+===========================+==========+===========+===============+
|  0  | N/A                       | primary  | N/A       | N/A           |
+-----+---------------------------+----------+-----------+---------------+
|  1  | SCI                       | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|  2  | ERR                       | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|  3  | DQ                        | IMAGE    | uint32    | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|  4  | VAR_POISSON               | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|  5  | VAR_RNOISE                | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|  6  | VAR_FLAT                  | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | AREA*                     | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | WAVELENGTH*               | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | PATHLOSS_PS*              | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | PATHLOSS_UN*              | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | BARSHADOW*                | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | ASDF                      | BINTABLE | N/A       | variable      |
+-----+---------------------------+----------+-----------+---------------+

 - SCI: 2-D data array containing the pixel values, in units of surface brightness.
 - ERR: 2-D data array containing uncertainty estimates for each pixel.
   These values are based on the combined VAR_POISSON and VAR_RNOISE data (see below),
   given as standard deviation.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - VAR_POISSON: 2-D data array containing the variance estimate for each pixel,
   based on Poisson noise only.
 - VAR_RNOISE: 2-D data array containing the variance estimate for each pixel,
   based on read noise only.
 - VAR_FLAT: 2-D data array containing the variance estimate for each pixel,
   based on uncertainty in the flat-field.
 - AREA: 2-D data array containing pixel area values, added by the :ref:`photom <photom_step>`
   step, for imaging modes.
 - WAVELENGTH: 2-D data array of wavelength values for each pixel, for some spectroscopic modes.
 - PATHLOSS_PS: 2-D data array of point-source pathloss correction factors, added by
   the :ref:`pathloss <pathloss_step>` step, for some spectroscopic modes.
 - PATHLOSS_UN: 1-D data array of uniform-source pathloss correction factors, added by
   the :ref:`pathloss <pathloss_step>` step, for some spectroscopic modes.
 - BARSHADOW: 2-D data array of NIRSpec MSA bar shadow correction factors, added by the
   :ref:`barshadow <barshadow_step>` step, for NIRSpec MOS exposures only.
 - ADSF: The data model meta data.

For spectroscopic modes that contain data for multiple sources, such as NIRSpec MOS,
NIRCam WFSS, and NIRISS WFSS, there will be multiple tuples of the SCI, ERR, DQ, VAR_POISSON,
VAR_RNOISE, etc. extensions, where each tuple contains the data for a given source or
slit, as created by the :ref:`extract_2d <extract_2d_step>` step. FITS "EXTVER" keywords are
used in each extension header to segregate the multiple instances of each extension type by
source.

.. _crf:
.. _crfints:

Cosmic-Ray flagged data: ``crf`` and ``crfints``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Several of the stage 3 pipelines, such as :ref:`calwebb_image3 <calwebb_image3>` and
:ref:`calwebb_spec3 <calwebb_spec3>`, include the :ref:`outlier detection <outlier_detection_step>`
step, which finds and flags outlier pixel values within calibrated images. The results of this
process have the identical format and content as the input ``cal`` and ``calints`` products.
The only difference is that the DQ arrays have been updated to contain CR flags. If the inputs
are in the form of ``cal`` products, the CR-flagged data will be saved to a ``crf`` product, which
has the exact same structure and content as the cal_ product described above. Similarly, if the
inputs are ``calints`` products, the CR-flagged results will be saved to a ``crfints`` product,
which has the same structure and content as the calints_ product described above.

.. _i2d:
.. _s2d:

Resampled 2-D data: ``i2d`` and ``s2d``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Images and spectra that have been resampled by the :ref:`resample <resample_step>` step use a
different set of data arrays than other science products. Resampled 2-D images are stored in
``i2d`` products and resampled 2-D spectra are stored in ``s2d`` products.
The FITS file structure for ``i2d`` and ``s2d`` products is as follows:

+-----+-------------+----------+-----------+-------------------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions              |
+=====+=============+==========+===========+=========================+
|  0  | N/A         | primary  | N/A       | N/A                     |
+-----+-------------+----------+-----------+-------------------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows           |
+-----+-------------+----------+-----------+-------------------------+
|  2  | ERR         | IMAGE    | float32   | ncols x nrows           |
+-----+-------------+----------+-----------+-------------------------+
|  3  | CON         | IMAGE    | int32     | ncols x nrows x nplanes |
+-----+-------------+----------+-----------+-------------------------+
|  4  | WHT         | IMAGE    | float32   | ncols x nrows           |
+-----+-------------+----------+-----------+-------------------------+
|  5  | VAR_POISSON | IMAGE    | float32   | ncols x nrows           |
+-----+-------------+----------+-----------+-------------------------+
|  6  | VAR_RNOISE  | IMAGE    | float32   | ncols x nrows           |
+-----+-------------+----------+-----------+-------------------------+
|  7  | VAR_FLAT    | IMAGE    | float32   | ncols x nrows           |
+-----+-------------+----------+-----------+-------------------------+
|     | HDRTAB*     | BINTABLE | N/A       | variable                |
+-----+-------------+----------+-----------+-------------------------+
|     | ASDF        | BINTABLE | N/A       | variable                |
+-----+-------------+----------+-----------+-------------------------+

 - SCI: 2-D data array containing the pixel values, in units of surface brightness
 - ERR: 2-D data array containing resampled uncertainty estimates, given as standard deviation
 - CON: 3-D context image, which encodes information about which input images contribute
   to a specific output pixel
 - WHT: 2-D weight image giving the relative weight of the output pixels
 - VAR_POISSON: 2-D resampled Poisson variance estimates for each pixel
 - VAR_RNOISE: 2-D resampled read noise variance estimates for each pixel
 - VAR_FLAT: 2-D resampled flat-field variance estimates for each pixel
 - HDRTAB: A table containing meta data (FITS keyword values) for all of the input images
   that were combined to produce the output image. Only appears when multiple inputs are used.
 - ADSF: The data model meta data.

For spectroscopic exposure-based products that contain spectra for more than one source or slit
(e.g. NIRSpec MOS) there will be multiple tuples of the SCI, ERR, CON, WHT, and variance
extensions, one set for each source or slit. FITS "EXTVER" keywords are used in each
extension header to segregate the multiple instances of each extension type by
source.

For the context array, CON, though the schema represents it as an ``int32``,
users should interpret and recast the array as ``uint32`` post-processing. This
inconsistency will be dealt with in a later release.

.. _s3d:

Resampled 3-D (IFU) data: ``s3d``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
3-D IFU cubes created by the :ref:`cube_build <cube_build_step>` step are stored in FITS
files with the following structure:

+-----+-------------+----------+-----------+------------------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions             |
+=====+=============+==========+===========+========================+
|  0  | N/A         | primary  | N/A       | N/A                    |
+-----+-------------+----------+-----------+------------------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows x nwaves |
+-----+-------------+----------+-----------+------------------------+
|  2  | ERR         | IMAGE    | float32   | ncols x nrows x nwaves |
+-----+-------------+----------+-----------+------------------------+
|  3  | DQ          | IMAGE    | uint32    | ncols x nrows x nwaves |
+-----+-------------+----------+-----------+------------------------+
|  4  | WMAP        | IMAGE    | float32   | ncols x nrows x nwaves |
+-----+-------------+----------+-----------+------------------------+
|     | WCS-TABLE   | BINTABLE | N/A       | 2 cols x 1 row         |
+-----+-------------+----------+-----------+------------------------+
|     | HDRTAB*     | BINTABLE | N/A       | variable               |
+-----+-------------+----------+-----------+------------------------+
|     | ASDF        | BINTABLE | N/A       | variable               |
+-----+-------------+----------+-----------+------------------------+

 - SCI: 3-D data array containing the spaxel values, in units of surface brightness.
 - ERR: 3-D data array containing uncertainty estimates for each spaxel.
 - DQ: 3-D data array containing DQ flags for each spaxel.
 - WMAP: 3-D weight image giving the relative weights of the output spaxels.
 - WCS-TABLE: A table listing the wavelength to be associated with each plane of the
   third axis in the SCI, DQ, ERR, and WMAP arrays, in a format that conforms to the
   FITS spectroscopic WCS standards. Column 1 of the table ("nelem") gives the number of
   wavelength elements listed in the table and column 2 ("wavelength") is a 1-D array
   giving the wavelength values.
 - HDRTAB: A table containing meta data (FITS keyword values) for all of the input images
   that were combined to produce the output image. Only appears when multiple inputs are used.
 - ADSF: The data model meta data.

``s3d`` products contain several unique meta data elements intended to aid in the use
of these products in data analysis tools. This includes the following keywords located in
the header of the FITS primary HDU:

 - FLUXEXT: A string value containing the EXTNAME of the extension containing the IFU flux
   data. Normally set to "SCI" for JWST IFU cube products.
 - ERREXT: A string value containing the EXTNAME of the extension containing error estimates
   for the IFU cube. Normally set to "ERR" for JWST IFU cube products.
 - ERRTYPE: A string value giving the type of error estimates contained in ERREXT, with
   possible values of "ERR" (error = standard deviation), "IERR" (inverse error),
   "VAR" (variance), and "IVAR" (inverse variance). Normally set to "ERR" for JWST IFU
   cube products.
 - MASKEXT: A string value containing the EXTNAME of the extension containing the Data Quality
   mask for the IFU cube. Normally set to "DQ" for JWST IFU cube products.

In addition, the following WCS-related keywords are included in the header of the "SCI"
extension to support the use of the wavelength table contained in the "WCS-TABLE" extension.
These keywords allow data analysis tools that are compliant with the FITS spectroscopic WCS
standards to automatically recognize and load the wavelength information in the "WCS-TABLE"
and assign wavelengths to the IFU cube data.

 - PS3_0 = 'WCS-TABLE': The name of the extension containing coordinate data for axis 3.
 - PS3_1 = 'wavelength': The name of the table column containing the coordinate data.

The coordinate data (wavelength values in this case) contained in the "WCS-TABLE" override
any coordinate information normally computed from FITS WCS keywords like CRPIX3, CRVAL3,
and CDELT3 for coordinate axis 3.

.. _x1d:
.. _x1dints:

Extracted 1-D spectroscopic data: ``x1d`` and ``x1dints``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Extracted spectral data produced by the :ref:`extract_1d <extract_1d_step>` step are stored
in binary table extensions of FITS files. The overall layout of the FITS file is as follows:

+-----+-------------+----------+-----------+---------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions    |
+=====+=============+==========+===========+===============+
|  0  | N/A         | primary  | N/A       | N/A           |
+-----+-------------+----------+-----------+---------------+
|  1  | EXTRACT1D   | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+
|  2  | ASDF        | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+

 - EXTRACT1D: A 2-D table containing the extracted spectral data.
 - ADSF: The data model meta data.

Multiple "EXTRACT1D" extensions can be present if there is data for more than one source or
if the file is an ``x1dints`` product. For ``x1dints`` products, there is one "EXTRACT1D"
extension for each integration in the exposure.

The structure of the "EXTRACT1D" table extension is as follows:

+-------------------+-----------+--------------------+---------------+
| Column Name       | Data Type | Contents           | Units         |
+===================+===========+===================+================+
| WAVELENGTH        | float64   | Wavelength values  | :math:`\mu` m |
+-------------------+-----------+--------------------+---------------+
| FLUX              | float64   | Flux values        | Jy            |
+-------------------+-----------+--------------------+---------------+
| FLUX_ERROR        | float64   | Error values       | Jy            |
+-------------------+-----------+--------------------+---------------+
| FLUX_VAR_POISSON  | float64   | Error values       | Jy^2          |
+-------------------+-----------+--------------------+---------------+
| FLUX_VAR_RNOISE   | float64   | Error values       | Jy^2          |
+-------------------+-----------+--------------------+---------------+
| FLUX_VAR_FLAT     | float64   | Error values       | Jy^2          |
+-------------------+-----------+--------------------+---------------+
| SURF_BRIGHT       | float64   | Surface Brightness | MJy/sr        |
+-------------------+-----------+--------------------+---------------+
| SB_ERROR          | float64   | Surf. Brt. errors  | MJy/sr        |
+-------------------+-----------+--------------------+---------------+
| SB_VAR_POISSON    | float64   | Surf. Brt. errors  | (MJy/sr)^2    |
+-------------------+-----------+--------------------+---------------+
| SB_VAR_RNOISE     | float64   | Surf. Brt. errors  | (MJy/sr)^2    |
+-------------------+-----------+--------------------+---------------+
| SB_VAR_FLAT       | float64   | Surf. Brt. errors  | (MJy/sr)^2    |
+-------------------+-----------+--------------------+---------------+
| DQ                | uint32    | DQ flags           | N/A           |
+-------------------+-----------+--------------------+---------------+
| BACKGROUND        | float64   | Background signal  | MJy/sr        |
+-------------------+-----------+--------------------+---------------+
| BKGD_ERROR        | float64   | Background error   | MJy/sr        |
+-------------------+-----------+--------------------+---------------+
| BKGD_VAR_POISSON  | float64   | Background error   | (MJy/sr)^2    |
+-------------------+-----------+--------------------+---------------+
| BKGD_VAR_RNOISE   | float64   | Background error   | (MJy/sr)^2    |
+-------------------+-----------+--------------------+---------------+
| BKGD_VAR_FLAT     | float64   | Background error   | (MJy/sr)^2    |
+-------------------+-----------+--------------------+---------------+
| NPIXELS           | float64   | Number of pixels   | N/A           |
+-------------------+-----------+--------------------+---------------+

The table is constructed using a simple 2-D layout, using one row per extracted spectral
element in the dispersion direction of the data (i.e. one row per wavelength bin).
Note that for point sources observed with NIRSpec or NIRISS SOSS mode, it is not
possible to express the extracted spectrum as surface brightness and hence the
SURF_BRIGHT and SB_ERROR columns will be set to zero. NPIXELS gives the (fractional)
number of pixels included in the source extraction region at each wavelength bin.

.. _c1d:

Combined 1-D spectroscopic data: ``c1d``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Combined spectral data produced by the :ref:`combine_1d <combine_1d_step>` step are stored
in binary table extensions of FITS files. The overall layout of the FITS file is as follows:

+-----+-------------+----------+-----------+---------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions    |
+=====+=============+==========+===========+===============+
|  0  | N/A         | primary  | N/A       | N/A           |
+-----+-------------+----------+-----------+---------------+
|  1  | COMBINE1D   | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+
|  2  | ASDF        | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+

 - COMBINE1D: A 2-D table containing the combined spectral data.
 - ADSF: The data model meta data.

The structure of the "COMBINE1D" table extension is as follows:

+-------------+-----------+--------------------+----------------+
| Column Name | Data Type | Contents           | Units          |
+=============+===========+====================+================+
| WAVELENGTH  | float64   | Wavelength values  | :math:`\mu` m  |
+-------------+-----------+--------------------+----------------+
| FLUX        | float64   | Flux values        | Jy             |
+-------------+-----------+--------------------+----------------+
| ERROR       | float64   | Error values       | Jy             |
+-------------+-----------+--------------------+----------------+
| SURF_BRIGHT | float64   | Surface Brightness | MJy/sr         |
+-------------+-----------+--------------------+----------------+
| SB_ERROR    | float64   | Surf. Brt. errors  | MJy/sr         |
+-------------+-----------+--------------------+----------------+
| DQ          | uint32    | DQ flags           | N/A            |
+-------------+-----------+--------------------+----------------+
| WEIGHT      | float64   | Sum of weights     | N/A            |
+-------------+-----------+--------------------+----------------+
| N_INPUT     | float64   | Number of inputs   | N/A            |
+-------------+-----------+--------------------+----------------+

The table is constructed using a simple 2-D layout, using one row per extracted spectral
element in the dispersion direction of the data (i.e. one row per wavelength bin).

.. _cat:

Source catalog: ``cat``
^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`source_catalog <source_catalog_step>` step contained in the
:ref:`calwebb_image3 <calwebb_image3>` pipeline detects and quantifies sources within imaging
products. The derived data for the sources is stored in a ``cat`` product, which is in the form
of an ASCII table in `ECSV <http://docs.astropy.org/en/stable/_modules/astropy/io/ascii/ecsv.html>`_
(Enhanced Character Separated Values) format. It is a flat text file, containing meta data
header entries and the source data in a 2-D table layout, with one row per source.

.. _segm:

Segmentation map: ``segm``
^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`source_catalog <source_catalog_step>` step contained in the
:ref:`calwebb_image3 <calwebb_image3>` pipeline uses an image segmentation procedure
to detect sources, which is a process of assigning a label to every image pixel that
contains signal from a source, such that pixels belonging to the same source have the
same label. The result of this procedure is saved in a ``segm`` product. The product
is in FITS format, with a single image extension containing a 2-D image. The image
has the same dimensions as the science image from which the sources were detected,
and each pixel belonging to a source has an integer value corresponding to the
label listed in the source catalog (``cat`` product).
Pixels not belonging to a source have a value of zero.

.. _phot:

Photometry catalog: ``phot``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`tso_photometry <tso_photometry_step>` step in the :ref:`calwebb_tso3 <calwebb_tso3>`
pipeline produces light curve from TSO imaging observations by computing aperture photometry as a
function of integration time stamp within one or more exposures. The resulting photometric data
are stored in a ``phot`` product, which is in the form of an ASCII table in
`ECSV <http://docs.astropy.org/en/stable/_modules/astropy/io/ascii/ecsv.html>`_
(Enhanced Character Separated Values) format. It is a flat text file, containing meta data
header entries and the photometric data in a 2-D table layout, with one row per exposure
integration.

.. _whtlt:

White-light photometry catalog: ``whtlt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`white_light <white_light_step>` step in the :ref:`calwebb_tso3 <calwebb_tso3>`
pipeline produces a light curve from TSO spectroscopic observations by computing the
wavelength-integrated spectral flux as a function of integration time stamp within one or more
exposures. The resulting data
are stored in a ``whtlt`` product, which is in the form of an ASCII table in
`ECSV <http://docs.astropy.org/en/stable/_modules/astropy/io/ascii/ecsv.html>`_
(Enhanced Character Separated Values) format. It is a flat text file, containing meta data
header entries and the white-light flux data in a 2-D table layout, with one row per exposure
integration.

.. _psfstack:

Stacked PSF data: ``psfstack``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`stack_refs <stack_refs_step>` step in the :ref:`calwebb_coron3 <calwebb_coron3>`
pipeline takes a collection of PSF reference image and assembles them into a 3-D stack of
PSF images, which results in a ``psfstack`` product. The ``psfstack`` product uses the
`~jwst.datamodels.CubeModel` data model, which when serialized to a FITS file has the
structure shown below.

+-----+-------------+----------+-----------+-----------------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions            |
+=====+=============+==========+===========+=======================+
|  0  | N/A         | primary  | N/A       | N/A                   |
+-----+-------------+----------+-----------+-----------------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows x npsfs |
+-----+-------------+----------+-----------+-----------------------+
|  2  | DQ          | IMAGE    | uint32    | ncols x nrows x npsfs |
+-----+-------------+----------+-----------+-----------------------+
|  3  | ERR         | IMAGE    | float32   | ncols x nrows x npsfs |
+-----+-------------+----------+-----------+-----------------------+
|  4  | ASDF        | BINTABLE | N/A       | variable              |
+-----+-------------+----------+-----------+-----------------------+

 - SCI: 3-D data array containing a stack of 2-D PSF images.
 - DQ: 3-D data array containing DQ flags for each PSF image.
 - ERR: 3-D data array containing a stack of 2-D uncertainty estimates for each PSF image.
 - ADSF: The data model meta data.

.. _psfalign:

Aligned PSF data: ``psfalign``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`align_refs <align_refs_step>` step in the :ref:`calwebb_coron3 <calwebb_coron3>`
pipeline creates a 3-D stack of PSF images that are aligned to corresponding science target
images. The resulting ``psfalign`` product uses the `~jwst.datamodels.QuadModel` data model,
which when serialized to a FITS file has the structure and content shown below.

+-----+-------------+----------+-----------+-------------------------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions                    |
+=====+=============+==========+===========+===============================+
|  0  | N/A         | primary  | N/A       | N/A                           |
+-----+-------------+----------+-----------+-------------------------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows x npsfs x nints |
+-----+-------------+----------+-----------+-------------------------------+
|  2  | DQ          | IMAGE    | uint32    | ncols x nrows x npsfs x nints |
+-----+-------------+----------+-----------+-------------------------------+
|  3  | ERR         | IMAGE    | float32   | ncols x nrows x npsfs x nints |
+-----+-------------+----------+-----------+-------------------------------+
|  4  | ASDF        | BINTABLE | N/A       | variable                      |
+-----+-------------+----------+-----------+-------------------------------+

 - SCI: 4-D data array containing a stack of 2-D PSF images aligned to each integration
   within a corresponding science target exposure.
   each integration.
 - DQ: 4-D data array containing DQ flags for each PSF image.
 - ERR: 4-D data array containing a stack of 2-D uncertainty estimates for each PSF image,
   per science target integration.
 - ADSF: The data model meta data.

.. _psfsub:

PSF-subtracted data: ``psfsub``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`klip <klip_step>` step in the :ref:`calwebb_coron3 <calwebb_coron3>`
pipeline subtracts an optimized combination of PSF images from each integration in a
science target exposure. The resulting PSF-subtracted science exposure data uses the
`~jwst.datamodels.CubeModel` data model, which when serialized to a FITS file has the
structure shown below.

+-----+-------------+----------+-----------+-----------------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions            |
+=====+=============+==========+===========+=======================+
|  0  | N/A         | primary  | N/A       | N/A                   |
+-----+-------------+----------+-----------+-----------------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  2  | ERR         | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  3  | DQ          | IMAGE    | uint32    | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  4  | INT_TIMES   | BINTABLE | N/A       | nints (rows) x 7 cols |
+-----+-------------+----------+-----------+-----------------------+
|  5  | VAR_POISSON | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  6  | VAR_RNOISE  | IMAGE    | float32   | ncols x nrows x nints |
+-----+-------------+----------+-----------+-----------------------+
|  7  | ASDF        | BINTABLE | N/A       | variable              |
+-----+-------------+----------+-----------+-----------------------+

 - SCI: 3-D data array containing a stack of 2-D PSF-subtracted science target images, one per
   integration.
 - ERR: 3-D data array containing a stack of 2-D uncertainty estimates for each science target
   integration.
 - DQ: 3-D data array containing DQ flags for each science target integration.
 - INT_TIMES: A table of beginning, middle, and end time stamps for each integration in the
   exposure.
 - VAR_POISSON: 3-D data array containing the per-integration variance estimates for each pixel,
   based on Poisson noise only.
 - VAR_RNOISE: 3-D data array containing the per-integration variance estimates for each pixel,
   based on read noise only.
 - ADSF: The data model meta data.

.. _ami:
.. _amiavg:
.. _aminorm:

AMI data: ``ami``, ``amiavg``, and ``aminorm``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AMI derived data created by the :ref:`ami_analyze <ami_analyze_step>`,
:ref:`ami_average <ami_average_step>`, and :ref:`ami_normalize <ami_normalize_step>` steps,
as part of the :ref:`calwebb_ami3 <calwebb_ami3>` pipeline, are stored in FITS files that
contain a mixture of images and binary table extensions. The output format of all three
pipeline steps is the same, encapsulated within a `~jwst.datamodels.AmiLgModel` data model.
The overall layout of the corresponding FITS files is as follows:

+-----+-------------+----------+-----------+-----------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions      |
+=====+=============+==========+===========+=================+
|  0  | N/A         | primary  | N/A       | N/A             |
+-----+-------------+----------+-----------+-----------------+
|  1  | FIT         | IMAGE    | float32   | ncols x nrows   |
+-----+-------------+----------+-----------+-----------------+
|  2  | RESID       | IMAGE    | float32   | ncols x nrows   |
+-----+-------------+----------+-----------+-----------------+
|  3  | CLOSURE_AMP | BINTABLE | float64   | 1 col x 35 rows |
+-----+-------------+----------+-----------+-----------------+
|  4  | CLOSURE_PHA | BINTABLE | float64   | 1 col x 35 rows |
+-----+-------------+----------+-----------+-----------------+
|  5  | FRINGE_AMP  | BINTABLE | float64   | 1 col x 21 rows |
+-----+-------------+----------+-----------+-----------------+
|  6  | FRINGE_PHA  | BINTABLE | float64   | 1 col x 21 rows |
+-----+-------------+----------+-----------+-----------------+
|  7  | PUPIL_PHA   | BINTABLE | float64   | 1 col x 7 rows  |
+-----+-------------+----------+-----------+-----------------+
|  8  | SOLNS       | BINTABLE | float64   | 1 col x 44 rows |
+-----+-------------+----------+-----------+-----------------+
|  9  | ASDF        | BINTABLE | N/A       | variable        |
+-----+-------------+----------+-----------+-----------------+

 - FIT: A 2-D image of the fitted model.
 - RESID: A 2-D image of the fit residuals.
 - CLOSURE_AMP: A table of closure amplitudes.
 - CLOSURE_PHA: A table of closure phases.
 - FRINGE_AMP: A table of fringe amplitudes.
 - FRINGE_PHA: A table of fringe phases.
 - PUPIL_PHA: A table of pupil phases.
 - SOLNS: A table of fringe coefficients.
 - ADSF: The data model meta data.

