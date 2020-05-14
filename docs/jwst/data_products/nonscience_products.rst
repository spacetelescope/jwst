Non-science products
--------------------

.. _dark:

Dark exposure: ``dark``
^^^^^^^^^^^^^^^^^^^^^^^
Dark exposures processed by the :ref:`calwebb_dark <calwebb_dark>` pipeline result in a
product that has the same structure and content as the :ref:`ramp <ramp>` product described above.
The details are as follows:

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
|  5  | GROUP      | BINTABLE | N/A       | variable                          |
+-----+------------+----------+-----------+-----------------------------------+
|  6  | INT_TIMES  | BINTABLE | N/A       | nints (rows) x 7 cols             |
+-----+------------+----------+-----------+-----------------------------------+
|     | ZEROFRAME* | IMAGE    | float32   | ncols x nrows x nints             |
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
 - GROUP: A table of meta data for some (or all) of the data groups.
 - INT_TIMES: A table of begining, middle, and end time stamps for each integration in the
   exposure.
 - ZEROFRAME: 3-D data array containing the pixel values of the zero-frame for each
   integration in the exposure, where each plane of the cube corresponds to a given integration.
   Only appears if the zero-frame data were requested to be downlinked separately.
 - REFOUT: The MIRI detector reference output values. Only appears in MIRI exposures.
 - ADSF: The data model meta data.

.. _trfld:

Charge trap state data: ``trapsfilled``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`persistence <persistence_step>` step in the :ref:`calwebb_detector1 <calwebb_detector1>`
pipeline produces an image containing information on the number of filled charge traps in each
pixel at the end of an exposure. Internally these data exist as a `~jwst.datamodels.TrapsFilledModel`
data model, which is saved to a ``trapsfilled`` FITS product. The FITS file has the following
format:

+-----+---------+----------+-----------+-------------------+
| HDU | EXTNAME | HDU Type | Data Type | Dimensions        |
+=====+=========+==========+===========+===================+
|  0  | N/A     | primary  | N/A       | N/A               |
+-----+---------+----------+-----------+-------------------+
|  1  | SCI     | IMAGE    | float32   | ncols x nrows x 3 |
+-----+---------+----------+-----------+-------------------+
|  2  | ASDF    | BINTABLE | N/A       | variable          |
+-----+---------+----------+-----------+-------------------+

 - SCI: 3-D data array giving the number of charge traps per pixel, with each plane
   corresponding to a different trap family.
 - ADSF: The data model meta data.

.. _wfscmb:

WFS&C combined image: ``wfscmb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`wfs_combine <wfs_combine_step>` step in the :ref:`calwebb_wfs-image3 <calwebb_wfs-image3>`
pipeline combines dithered pairs of Wavefront Sensing and Control (WFS&C) images, with the
result being stored in a ``wfscmb`` product. Unlike the drizzle methods used to combine and resample
science images, resulting in an :ref:`i2d <i2d>` product, the WFS&C combination is a simple shift and add
technique that results in a standard imaging FITS file structure, as shown below.

+-----+---------+----------+-----------+---------------+
| HDU | EXTNAME | HDU Type | Data Type | Dimensions    |
+=====+=========+==========+===========+===============+
|  0  | N/A     | primary  | N/A       | N/A           |
+-----+---------+----------+-----------+---------------+
|  1  | SCI     | IMAGE    | float32   | ncols x nrows |
+-----+---------+----------+-----------+---------------+
|  2  | ERR     | IMAGE    | float32   | ncols x nrows |
+-----+---------+----------+-----------+---------------+
|  3  | DQ      | IMAGE    | uint32    | ncols x nrows |
+-----+---------+----------+-----------+---------------+
|  4  | ASDF    | BINTABLE | N/A       | variable      |
+-----+---------+----------+-----------+---------------+

 - SCI: 2-D data array containing the pixel values, in units of surface brightness.
 - ERR: 2-D data array containing uncertainty estimates for each pixel.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - ADSF: The data model meta data.

