FITS file structures and contents
=================================

Here we describe the structure and content of the most frequently used
forms of FITS files for JWST science data products. Each type of FITS
file is the result of serialization of a corresponding data model. All
JWST pipeline input and output images, with the exception of a few
reference files, are serialized as FITS files. The asdf representation
of the datamodel is serialized as the final extension of the FITS
file, in the ASDF extension, as eight bit unsigned ints, essentialy a
text charavter serialization in yaml format.  The asdf representation
is read from the extension when the FITS file is read to provide the
initial datamodel. Values in the other FiTS extensions then either
override this initial model or are adde to it.


Common Features
---------------

All FITS science products have a few common features to their structure
and organization:

1. The primary Header-Data Unit (HDU) only contains header information,
   in the form of keyword records, with an empty data array, which is
   indicated by the occurence of NAXIS=0 in the primary header. Meta
   data that pertains to the entire product is stored in keywords in the
   primary header. Meta data related to specific extensions (see below)
   should be stored in keywords in the headers of those extensions.

2. All data related to the product are contained in one or more FITS
   Image or Table extensions. The header of each extension may contain
   keywords that pertain uniquely to that extension.

 3.
 
Level-1 and Level-2 exposure-based products, which contain the data
from an individual exposure on an individual detector, use the
following file naming scheme::

    jw{ppppp}{ooo}{vvv}_{gg}{s}{aa}_{eeeee}_{detector}_{suffix}.fits

where:

   - ppppp: program ID number
   - ooo: observation number
   - vvv: visit number
   - gg: visit group
   - s: parallel sequence ID (1=prime, 2-5=parallel)
   - aa: activity number (base 36)
   - eeeee: exposure number
   - detector: detector name (e.g. 'nrca1', 'nrcblong', 'mirimage')
   - suffix: product type identifier (e.g. 'uncal', 'rate', 'cal')

An example Level-2a product FITS file name is::

   jw93065002001_02101_00001_nrca1_rate.fits

Specific products
-----------------

This section lists the organization and contents of each type of
science product in FITS form.

Raw Level-1b (suffix = `uncal`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exposure raw data (level-1b) products are designated with a file name
suffix of "uncal." These files usually contain only the raw pixel values
from an exposure, with the addition of a table extension that contains
some downlinked meta data pertaining to individual groups. Additional
extensions can be included for certain instruments and readout types.
If the zero-frame was requested to be downlinked, an additional image
extension is included that contains those data.
MIRI exposures also contain an additional image extension with the values
from the reference output. The FITS file structure is as follows.

+-----+-------------------+-----------+----------+-----------+---------------------------------+
| HDU | Content           | EXTNAME   | HDU Type | Data Type | Dimensions                      |
+=====+===================+===========+==========+===========+=================================+
|  0  | Primary header    | N/A       | N/A      | N/A       | N/A                             |
+-----+-------------------+-----------+----------+-----------+---------------------------------+
|  1  | Pixel values      | SCI       | IMAGE    | uint16    | ncols x nrows x ngroups x nints |
+-----+-------------------+-----------+----------+-----------+---------------------------------+
|  2  | Group meta        | GROUP     | BINTABLE | N/A       | variable                        |
+-----+-------------------+-----------+----------+-----------+---------------------------------+
|  3  | Zero frame images | ZEROFRAME | IMAGE    | uint16    | ncols x nrows x nints           |
+-----+-------------------+-----------+----------+-----------+---------------------------------+
|  4  | Reference output  | REFOUT    | IMAGE    | uint16    | ncols x 256 x ngroups x nints   |
+-----+-------------------+-----------+----------+-----------+---------------------------------+

The raw pixel values in the SCI extension are stored as a 4-D data array,
having dimensions equal to the 2-D size of the detector readout, with
the data from the multiple groups (ngroups) within each integration stored
along the 3rd axis, and the multiple integrations (nints) stored along the
4th axis.

If zero-frame data are downlinked, there will be one zero-frame image
for each integration, stored as a 3-D cube (each cube plane corresponds
to an integration).

Level-2 ramp data (suffix = `ramp`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As soon as raw level-1b products are loaded into the calibration
pipeline the contents of the product is modified to include
additional data extensions, as well as converting the raw SCI
(and ZEROFRAME and REFOUT, if present) array values from integer to
floating-point data type. New data arrays that are added include an
ERR extension and two types of data quality flag extensions. There
is a 2-D `PIXELDQ` extension that will contain flags that pertain
to all groups and all integrations, and there is also a 4-D
`GROUPDQ` extension for containing flags that pertain to individual
groups within individual integrations. The FITS file layout is as
follows:

+-----+-------------------+-----------+----------+-----------+---------------------------------+
| HDU | Content           | EXTNAME   | HDU Type | Data Type | Dimensions                      |
+=====+===================+===========+==========+===========+=================================+
|  0  | Primary header    | N/A       | N/A      | N/A       | N/A                             |
+-----+-------------------+-----------+----------+-----------+---------------------------------+
|  1  | Pixel values      | SCI       | IMAGE    | float32   | ncols x nrows x ngroups x nints |
+-----+-------------------+-----------+----------+-----------+---------------------------------+
|  2  | 2-D data quality  | PIXELDQ   | IMAGE    | uint32    | ncols x nrows                   |
+-----+-------------------+-----------+----------+-----------+---------------------------------+
|  3  | 4-D data quality  | GROUPDQ   | IMAGE    | uint8     | ncols x nrows x ngroups x nints |
+-----+-------------------+-----------+----------+-----------+---------------------------------+
|  4  | Error values      | ERR       | IMAGE    | float32   | ncols x nrows x ngroups x nints |
+-----+-------------------+-----------+----------+-----------+---------------------------------+

Any additional extensions that were present in the raw level-1b
file (e.g. GROUP, ZEROFRAME, REFOUT) will be carried along and
will also appear in the level-2 ramp product.

Level-2a countrate products (suffix = `rate` and `rateints`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Countrate products are produced by applying ramp-fitting to the
integrations within an exposure, in order to compute count rates
from the original accumulating signal. For exposures that
contain multiple integrations (nints > 1) this is done in two ways,
which results in two separate products that are produced. First,
countrates are computed for each integration within the exposure,
the resuls of which are stored in a `rateints` product. These
products will contain 3-D science data arrays, where each plane
of the data cube contains the countrate image for an integration.

The results for each integration are also averaged together to
form a single 2-D countrate image for the entire exposure. These
resuls are stored in a `rate` product.

The FITS file structure for a `rateints` product is as follows:

+-----+----------------+---------+----------+-----------+-----------------------+
| HDU | Content        | EXTNAME | HDU Type | Data Type | Dimensions            |
+=====+================+=========+==========+===========+=======================+
|  0  | Primary header | N/A     | N/A      | N/A       | N/A                   |
+-----+----------------+---------+----------+-----------+-----------------------+
|  1  | Pixel values   | SCI     | IMAGE    | float32   | ncols x nrows x nints |
+-----+----------------+---------+----------+-----------+-----------------------+
|  2  | Data quality   | DQ      | IMAGE    | uint32    | ncols x nrows x nints |
+-----+----------------+---------+----------+-----------+-----------------------+
|  3  | Error values   | ERR     | IMAGE    | float32   | ncols x nrows x nints |
+-----+----------------+---------+----------+-----------+-----------------------+

The FITS file structure for a `rate` product is as follows:

+-----+----------------+---------+----------+-----------+---------------+
| HDU | Content        | EXTNAME | HDU Type | Data Type | Dimensions    |
+=====+================+=========+==========+===========+===============+
|  0  | Primary header | N/A     | N/A      | N/A       | N/A           |
+-----+----------------+---------+----------+-----------+---------------+
|  1  | Pixel values   | SCI     | IMAGE    | float32   | ncols x nrows |
+-----+----------------+---------+----------+-----------+---------------+
|  2  | Data quality   | DQ      | IMAGE    | uint32    | ncols x nrows |
+-----+----------------+---------+----------+-----------+---------------+
|  3  | Error values   | ERR     | IMAGE    | float32   | ncols x nrows |
+-----+----------------+---------+----------+-----------+---------------+

Note that the two separate forms of PIXELDQ and GROUPDQ flags from the
previous types of products have been combined into a single DQ extension
with the same dimensions as the SCI and ERR components.

Level-2b calibrated products (suffix = `cal` and `calints`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Single exposure calibrated products duplicate the format and content of
level-2a products. As with level-2a, there are two different forms of
calibrated products: one containing results for individual integrations
(`calints`) and one for exposure-wide results (`cal`).

The FITS file structure for a `calints` product is as follows:

+-----+----------------+---------+----------+-----------+-----------------------+
| HDU | Content        | EXTNAME | HDU Type | Data Type | Dimensions            |
+=====+================+=========+==========+===========+=======================+
|  0  | Primary header | N/A     | N/A      | N/A       | N/A                   |
+-----+----------------+---------+----------+-----------+-----------------------+
|  1  | Pixel values   | SCI     | IMAGE    | float32   | ncols x nrows x nints |
+-----+----------------+---------+----------+-----------+-----------------------+
|  2  | Data quality   | DQ      | IMAGE    | uint32    | ncols x nrows x nints |
+-----+----------------+---------+----------+-----------+-----------------------+
|  3  | Error values   | ERR     | IMAGE    | float32   | ncols x nrows x nints |
+-----+----------------+---------+----------+-----------+-----------------------+

The FITS file structure for a `cal` product is as follows:

+-----+----------------+---------+----------+-----------+---------------+
| HDU | Content        | EXTNAME | HDU Type | Data Type | Dimensions    |
+=====+================+=========+==========+===========+===============+
|  0  | Primary header | N/A     | N/A      | N/A       | N/A           |
+-----+----------------+---------+----------+-----------+---------------+
|  1  | Pixel values   | SCI     | IMAGE    | float32   | ncols x nrows |
+-----+----------------+---------+----------+-----------+---------------+
|  2  | Data quality   | DQ      | IMAGE    | uint32    | ncols x nrows |
+-----+----------------+---------+----------+-----------+---------------+
|  3  | Error values   | ERR     | IMAGE    | float32   | ncols x nrows |
+-----+----------------+---------+----------+-----------+---------------+
