JWST FITS data product formats
==============================

Here we describe the structure and content of the most frequently used
forms of FITS files for JWST science data products. Each type of FITS
file is the result of serialization of a corresponding data model. All
JWST pipeline input and output products, with the exception of a few
reference files and catalogs, are serialized as FITS files. The ASDF
representation of the data model is serialized as a FITS BINTABLE extension
within the FITS file, with EXTNAME="ASDF" (essentialy a text character
serialization in yaml format.  The ASDF representation
is read from the extension when a FITS file is loaded into a data model,
to provide the initial instance of the data model. Values in the other
FITS extensions then either override this initial model or are added to it.

Common Features
---------------

All JWST FITS science products have a few common features in their structure
and organization:

1. The FITS primary Header-Data Unit (HDU) only contains header information,
   in the form of keyword records, with an empty data array, which is
   indicated by the occurence of NAXIS=0 in the primary header. Meta
   data that pertains to the entire product is stored in keywords in the
   primary header. Meta data related to specific extensions (see below)
   is stored in keywords in the headers of each extension.

2. All data related to the product are contained in one or more FITS
   IMAGE or BINTABLE extensions. The header of each extension may contain
   keywords that pertain uniquely to that extension.

Stage 1 and 2 exposure-based products, which contain the data
from an individual exposure on an individual detector, use the
following FITS file naming scheme::

 jw<ppppp><ooo><vvv>_<gg><s><aa>_<eeeee>_<detector>_<prodType>.fits

where

 - ppppp: program ID number
 - ooo: observation number
 - vvv: visit number
 - gg: visit group
 - s: parallel sequence ID (1=prime, 2-5=parallel)
 - aa: activity number (base 36)
 - eeeee: exposure number
 - detector: detector name (e.g. 'nrca1', 'nrcblong', 'mirimage')
 - prodType: product type identifier (e.g. 'uncal', 'rate', 'cal')

An example Stage 1 product FITS file name is::

 jw93065002001_02101_00001_nrca1_rate.fits

Stage 3 products, which consist of data from multiple exposures and/or
detectors, use a source-based file naming scheme::

 jw<ppppp>-<AC_ID>_[<"t"TargID | "s"SourceID">](-<"epoch"X>)_<instr>_<optElements>(-<subarray>)_<prodType>(-<ACT_ID>).fits

where

 - ppppp: program ID number
 - AC_ID: association candidate ID
 - TargID: a 3-digit target ID (either TargID or SourceID must be present)
 - SourceID: a 5-digit source ID
 - epochX: the text "epoch" followed by a single digit epoch number
 - instr: science instrument name (e.g. 'nircam', 'miri')
 - optElements: a single or hyphen-separated list of optical elements (e.g. filter, grating)
 - subarray: optional indicator of subarray name
 - prodType: product type identifier (e.g. 'i2d', 's3d', 'x1d')
 - ACT_ID: a 2-digit activity ID

An example Stage 3 product FITS file name is::

 jw87600-a3001_t001_niriss_f480m-nrm_amiavg.fits

Specific products
-----------------

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

Uncalibrated raw data: ``uncal``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Exposure raw data products are designated by a file name
suffix of "uncal." These files usually contain only the raw detector pixel values
from an exposure, with the addition of some table extensions containing various types of
meta data associated with the exposure.
Additional extensions can be included for certain instruments and readout types, as noted
below.
The FITS file structure is as follows.

+-----+------------+----------+-----------+---------------------------------+
| HDU | EXTNAME    | HDU Type | Data Type | Dimensions                      |
+=====+============+==========+===========+=================================+
|  0  | N/A        | primary  | N/A       | N/A                             |
+-----+------------+----------+-----------+---------------------------------+
|  1  | SCI        | IMAGE    | uint16    | ncols x nrows x ngroups x nints |
+-----+------------+----------+-----------+---------------------------------+
|  2  | GROUP      | BINTABLE | N/A       | variable                        |
+-----+------------+----------+-----------+---------------------------------+
|  3  | INT_TIMES  | BINTABLE | N/A       | nints (rows) x 7 cols           |
+-----+------------+----------+-----------+---------------------------------+
|     | ZEROFRAME* | IMAGE    | uint16    | ncols x nrows x nints           |
+-----+------------+----------+-----------+---------------------------------+
|     | REFOUT*    | IMAGE    | uint16    | ncols x 256 x ngroups x nints   |
+-----+------------+----------+-----------+---------------------------------+
|     | ASDF       | BINTABLE | N/A       | variable                        |
+-----+------------+----------+-----------+---------------------------------+

 - SCI: 4-D data array containing the raw pixel values. The first two dimensions are equal to
   the size of the detector readout, with the data from multiple groups (NGROUPS) within each
   integration stored along the 3rd axis, and the multiple integrations (NINTS) stored along
   the 4th axis.
 - GROUP: A table of meta data for some (or all) of the data groups.
 - INT_TIMES: A table of begining, middle, and end time stamps for each integration in the
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

Ramp data: ``ramp``
^^^^^^^^^^^^^^^^^^^
As raw data progress through the :ref:`calwebb_detector1 <calwebb_detector1>` pipeline
they are stored internally in a `~jwst.datamodels.RampModel` (or `~jwst.datamodels.MIRIRampModel`
for MIRI exposures). This type of data model is serialized to a ``ramp`` type FITS
file on disk. The original detector pixel values (in the SCI extension) are converted
from integer to floating-point data type. The same is true for the ZEROFRAME and REFOUT
data extensions, if they are present. An ERR array and two types of data quality arrays are
also added to the product. The FITS file layout is as follows:

+-----+------------+----------+-----------+---------------------------------+
| HDU | EXTNAME    | HDU Type | Data Type | Dimensions                      |
+=====+============+==========+===========+=================================+
|  0  | N/A        | primary  | N/A       | N/A                             |
+-----+------------+----------+-----------+---------------------------------+
|  1  | SCI        | IMAGE    | float32   | ncols x nrows x ngroups x nints |
+-----+------------+----------+-----------+---------------------------------+
|  2  | PIXELDQ    | IMAGE    | uint32    | ncols x nrows                   |
+-----+------------+----------+-----------+---------------------------------+
|  3  | GROUPDQ    | IMAGE    | uint8     | ncols x nrows x ngroups x nints |
+-----+------------+----------+-----------+---------------------------------+
|  4  | ERR        | IMAGE    | float32   | ncols x nrows x ngroups x nints |
+-----+------------+----------+-----------+---------------------------------+
|     | ZEROFRAME* | IMAGE    | float32   | ncols x nrows x nints           |
+-----+------------+----------+-----------+---------------------------------+
|     | GROUP      | BINTABLE | N/A       | variable                        |
+-----+------------+----------+-----------+---------------------------------+
|     | INT_TIMES  | BINTABLE | N/A       | nints (rows) x 7 cols           |
+-----+------------+----------+-----------+---------------------------------+
|     | REFOUT*    | IMAGE    | uint16    | ncols x 256 x ngroups x nints   |
+-----+------------+----------+-----------+---------------------------------+
|     | ASDF       | BINTABLE | N/A       | variable                        |
+-----+------------+----------+-----------+---------------------------------+

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
 - INT_TIMES: A table of begining, middle, and end time stamps for each integration in the
   exposure.
 - REFOUT: The MIRI detector reference output values. Only appears in MIRI exposures.
 - ADSF: The data model meta data.

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
image for the entire exposure. These resuls are stored in a ``rate`` product.

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
 - INT_TIMES: A table of begining, middle, and end time stamps for each integration in the
   exposure.
 - VAR_POISSON: 3-D data array containing the per-integration variance estimates for each pixel,
   based on Poisson noise only.
 - VAR_RNOISE: 3-D data array containing the per-integration variance estimates for each pixel,
   based on read noise only.
 - ADSF: The data model meta data.

These FITS files are compatitable with the `~jwst.datamodels.CubeModel` data model.

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
|     | AREA*       | IMAGE    |           | ncols x nrows         |
+-----+-------------+----------+-----------+-----------------------+
|     | RELSENS*    | BINTABLE | N/A       | variable              |
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
 - INT_TIMES: A table of begining, middle, and end time stamps for each integration in the
   exposure.
 - VAR_POISSON: 3-D data array containing the variance estimate for each pixel,
   based on Poisson noise only, for each integration.
 - VAR_RNOISE: 3-D data array containing the variance estimate for each pixel,
   based on read noise only, for each integration.
 - AREA: 2-D data array containing pixel area values, added by the :ref:`photom <photom_step>`
   step, for imaging modes.
 - RELSENS: A table of sensitivity values as a function of wavelength, added by the
   :ref:`photom <photom_step>` step, for some spectroscopic modes.
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
|     | VAR_POISSON               | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | VAR_RNOISE                | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | AREA*                     | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | RELSENS*                  | BINTABLE | N/A       | variable      |
+-----+---------------------------+----------+-----------+---------------+
|     | RELSENS2D*                | BINTABLE | N/A       | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | PATHLOSS_POINTSOURCE*     | IMAGE    | float32   | ncols         |
+-----+---------------------------+----------+-----------+---------------+
|     | WAVELENGTH_POINTSOURCE*   | IMAGE    | float32   | ncols         |
+-----+---------------------------+----------+-----------+---------------+
|     | PATHLOSS_UNIFORMSOURCE*   | IMAGE    | float32   | ncols         |
+-----+---------------------------+----------+-----------+---------------+
|     | WAVELENGTH_UNIFORMSOURCE* | IMAGE    | float32   | ncols         |
+-----+---------------------------+----------+-----------+---------------+
|     | BARSHADOW*                | IMAGE    | float32   | ncols x nrows |
+-----+---------------------------+----------+-----------+---------------+
|     | WAVELENGTH*               | IMAGE    | float32   | ncols x nrows |
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
 - AREA: 2-D data array containing pixel area values, added by the :ref:`photom <photom_step>`
   step, for imaging modes.
 - RELSENS: A table of sensitivity values as a function of wavelength, added by the
   :ref:`photom <photom_step>` step, for some spectroscopic modes.
 - RELSENS2D: 2-D data array of sensitivity values per pixel, added by the
   :ref:`photom <photom_step>` step, for IFU spectroscopic modes.
 - PATHLOSS_POINTSOURCE: 1-D data array of point-source pathloss correction factors, added by
   the :ref:`pathloss <pathloss_step>` step, for some spectroscopic modes.
 - WAVELENGTH_POINTSOURCE: 1-D data array of wavelength values associated with the
   PATHLOSS_POINTSOURCE correction factors, added by the :ref:`pathloss <pathloss_step>` step,
   for some spectroscopic modes.
 - PATHLOSS_UNIFORMSOURCE: 1-D data array of uniform-source pathloss correction factors, added by
   the :ref:`pathloss <pathloss_step>` step, for some spectroscopic modes.
 - WAVELENGTH_UNIFORMSOURCE: 1-D data array of wavelength values associated with the
   PATHLOSS_UNIFORMSOURCE correction factors, added by the :ref:`pathloss <pathloss_step>` step,
   for some spectroscopic modes.
 - BARSHADOW: 2-D data array of NIRSpec MSA bar shadow correction factors, added by the
   :ref:`barshadow <barshadow_step>` step, for NIRSpec MSA exposures only.
 - WAVELENGTH: 2-D data array of wavelength values for each pixel, for some spectroscopic modes.
 - ADSF: The data model meta data.

For spectroscopic modes that contain data for multiple sources, such as NIRSpec MOS,
NIRCam WFSS, and NIRISS WFSS, there will be multiple tuples of the SCI, ERR, DQ, VAR_POISSON,
VAR_RNOISE, RELSENS, etc. extensions, where each tuple contains the data for a given source or
slit, as created by the :ref:`extract_2d <extract_2d_step>` step. FITS "EXTVER" keywords are
used in each extension header to segregate the multiple instances of each extension type by
source.

Resampled 2-D data: ``i2d`` and ``s2d``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Images and spectra that have been resampled by the :ref:`resample <resample_step>` step use a
different set of data arrays than other science products. Resampled 2-D images are stored in
``i2d`` products and resampled 2-D spectra are stored in ``s2d`` products.
The FITS file structure for ``i2d`` products is as follows:

+-----+-------------+----------+-----------+---------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions    |
+=====+=============+==========+===========+===============+
|  0  | N/A         | primary  | N/A       | N/A           |
+-----+-------------+----------+-----------+---------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows |
+-----+-------------+----------+-----------+---------------+
|  2  | CON         | IMAGE    | int32     | ncols x nrows |
+-----+-------------+----------+-----------+---------------+
|  3  | WHT         | IMAGE    | float32   | ncols x nrows |
+-----+-------------+----------+-----------+---------------+
|     | HDRTAB*     | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+
|     | ASDF        | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+

 - SCI: 2-D data array containing the pixel values, in units of surface brightness.
 - CON: 2-D context image, which encodes information about which input images contribute
   to a specific output pixel.
 - WHT: 2-D weight image giving the relative weight of the output pixels (effectively a
   relative exposure time map).
 - HDRTAB: A table containing meta data (FITS keyword values) for all of the input images
   that were combined to produce the output image. Only appears when multiple inputs are used.
 - ADSF: The data model meta data.
 
The FITS file structure for ``s2d`` products is as follows:

+-----+-------------+----------+-----------+---------------+
| HDU | EXTNAME     | HDU Type | Data Type | Dimensions    |
+=====+=============+==========+===========+===============+
|  0  | N/A         | primary  | N/A       | N/A           |
+-----+-------------+----------+-----------+---------------+
|  1  | SCI         | IMAGE    | float32   | ncols x nrows |
+-----+-------------+----------+-----------+---------------+
|  2  | CON         | IMAGE    | int32     | ncols x nrows |
+-----+-------------+----------+-----------+---------------+
|  3  | WHT         | IMAGE    | float32   | ncols x nrows |
+-----+-------------+----------+-----------+---------------+
|  4  | RELSENS     | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+
|     | HDRTAB*     | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+
|     | ASDF        | BINTABLE | N/A       | variable      |
+-----+-------------+----------+-----------+---------------+

 - SCI: 2-D data array containing the pixel values, in units of surface brightness.
 - CON: 2-D context image, which encodes information about which input images contribute
   to a specific output pixel.
 - WHT: 2-D weight image giving the relative weight of the output pixels (effectively a
   relative exposure time map).
 - RELSENS: A table of sensitivity values as a function of wavelength, carried along from the
   input calibrated product.
 - HDRTAB: A table containing meta data (FITS keyword values) for all of the input images
   that were combined to produce the output image. Only appears when multiple inputs are used.
 - ADSF: The data model meta data.

For exposure-based products that contain spectra for more than one source or slit
(e.g. NIRSpec MOS) there will be multiple tuples of the SCI, CON, WHT, and RELSENS
extensions, one set for each source or slit. FITS "EXTVER" keywords are used in each
extension header to segregate the multiple instances of each extension type by
source.

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
|  3  | DQ          | IMAGE    | int32     | ncols x nrows x nwaves |
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
 - DQ: 3-D data array containing DQ flags for each spaxel.
 - ERR: 3-D data array containing uncertainty estimates for each spaxel.
 - WMAP: 3-D weight image giving the relative weight of the output spaxels.
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

+-------------+-----------+-------------------+----------------+
| Column Name | Data Type | Contents          | Units          |
+=============+===========+===================+================+
| WAVELENGTH  | float64   | Wavelength values | :math:`\mu` m  |
+-------------+-----------+-------------------+----------------+
| FLUX        | float64   | Flux values       | mJy            |
+-------------+-----------+-------------------+----------------+
| ERROR       | float64   | Error values      | mJy            |
+-------------+-----------+-------------------+----------------+
| DQ          | int32     | DQ flags          | N/A            |
+-------------+-----------+-------------------+----------------+
| NET         | float64   | Net flux          | DN/s           |
+-------------+-----------+-------------------+----------------+
| NERROR      | float64   | Net error         | DN/s           |
+-------------+-----------+-------------------+----------------+
| BACKGROUND  | float64   | Background signal | DN/s           |
+-------------+-----------+-------------------+----------------+
| BERROR      | float64   | Background error  | DN/s           |
+-------------+-----------+-------------------+----------------+

The table is constructed using a simple 2-D layout, using one row per extracted spectral
element in the dispersion direction of the data (i.e. one row per wavelength bin).
For some spectroscopic modes, such as MIRI MRS and NIRSpec IFU, the data that are used
as input to the :ref:`extract_1d` step are already in calibrated units of surface
brightness and therefore it's not possible to present 1-D extracted results for the net
and background spectral data in units of DN/s. In these cases the NET, NERROR, BACKGROUND,
and BERROR table columns will be zero-filled and only the WAVELENGTH, FLUX, ERROR, and DQ
columns will be populated.
