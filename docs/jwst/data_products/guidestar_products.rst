.. _guidestar_products:

Guide star data products
------------------------
The FGS guiding capabilities are provided by 4 guide star functions: Identification, Acquisition,
Track, and Fine Guide. Data downlinked by these functions is processed by DMS to provide uncalibrated
and calibrated guide star data products. The uncalibrated products consist of raw pixel value data
arrays, as well as different kinds of tabular information related to the guide stars and centroid
locations of the guide star as computed by the on-board FGS Flight Software (FSW). Calibrated
guide star products are created by the :ref:`calwebb_guider <calwebb_guider>` pipeline. Briefly, the
processing performed applies bad pixel masks and flat-fields the science data, as well as computing
countrate images from the multiple groups within each integration contained in a given product. The
countrate images, which are computed by the :ref:`guider_cds <guider_cds_step>` step, are computed
for most modes by simply differencing groups 1 and 2 in each integration and dividing by the group
time.

ID mode
^^^^^^^
The "Identification" guiding function images the field of view by reading the detector in a series
subarray "strips" that, collectively, cover most of the field. A total of 36 subarray strips are
read out, each of which is 64 x 2048 pixels in size. Each strip has 8 pixels of overlap with its
adjoining strips, resulting in a total of 2024 unique detector rows that've been read out. ID mode
uses 2 groups per integration and 2 integrations, resulting in a total of 4 reads. Each subarray
strip has its 4 reads performed before moving on to the next subarray.

DMS creates 2 different forms of products for ID mode data: one in which an image is constructed
by simply stacking or butting the data from adjacent subarray strips against one another and the
other in which the overlap regions of the strips are taken into account by averaging the pixel
values. The first form is referred to as a "stacked" product and the second as an "image" product.

The file name syntax for ID mode products is as follows::

 jw<pppppooovvv>_gs-id_<m>_<stacked | image>-<uncal | cal>.fits

where:

 - pppppooovvv: the visit id, consisting of program id, obervation number, and visit number
 - gs-id: guide star identification product
 - m: the number of the identification attempt (0-8) within a visit
 - stacked: indicates that the subarray strips are stacked (butted) against one another in the image
 - image: indicates that the subarray strips have been merged in their overlap regions
 - uncal: indicates the uncalibrated (raw) product
 - cal: indicates the calibrated product

The FITS file structure for uncalibrated ID "image" products is as follows:

+-----+-------------------------+----------+-----------+---------------------+
| HDU | EXTNAME                 | HDU Type | Data Type | Dimensions          |
+=====+=========================+==========+===========+=====================+
|  0  | N/A                     | primary  | N/A       | N/A                 |
+-----+-------------------------+----------+-----------+---------------------+
|  1  | SCI                     | IMAGE    | uint16    | 2024 x 2048 x 2 x 2 |
+-----+-------------------------+----------+-----------+---------------------+
|  2  | Flight Reference Stars  | BINTABLE | N/A       | 4 cols x nstars     |
+-----+-------------------------+----------+-----------+---------------------+
|  3  | Planned Reference Stars | BINTABLE | N/A       | 10 cols x nstars    |
+-----+-------------------------+----------+-----------+---------------------+

 - SCI: 4-D data array containing the raw pixel values. The subarray overlaps have been accounted for,
   resulting in image dimensions of 2024 x 2048 pixels, with the 2 groups and 2 integrations stacked
   along the 3rd and 4th array axes.
 - Flight Reference Stars: A table containing information on the actual reference stars
   used by the FSW. Detailed contents are listed :ref:`below <flight_ref_stars>`.
 - Planned Reference Stars: A table containing information on the planned reference stars.
   Detailed contents are listed :ref:`below <planned_ref_stars>`.

The FITS file structure for uncalibrated ID "stacked" products is as follows:

+-----+-------------------------+----------+-----------+---------------------+
| HDU | EXTNAME                 | HDU Type | Data Type | Dimensions          |
+=====+=========================+==========+===========+=====================+
|  0  | N/A                     | primary  | N/A       | N/A                 |
+-----+-------------------------+----------+-----------+---------------------+
|  1  | SCI                     | IMAGE    | uint16    | 2304 x 2048 x 2 x 2 |
+-----+-------------------------+----------+-----------+---------------------+
|  2  | Flight Reference Stars  | BINTABLE | N/A       | 4 cols x nstars     |
+-----+-------------------------+----------+-----------+---------------------+
|  3  | Planned Reference Stars | BINTABLE | N/A       | 10 cols x nstars    |
+-----+-------------------------+----------+-----------+---------------------+

 - SCI: 4-D data array containing the raw pixel values. The subarray data are butted against one
   another, resulting in image dimensions of 2304 x 2048 pixels, with the 2 groups and 2 integrations
   stacked along the 3rd and 4th array axes.
 - Flight Reference Stars: A table containing information on the actual reference stars
   used by the FSW. Detailed contents are listed :ref:`below <flight_ref_stars>`.
 - Planned Reference Stars: A table containing information on the planned reference stars.
   Detailed contents are listed :ref:`below <planned_ref_stars>`.

The FITS file structure for calibrated ID "image" products is as follows:

+-----+-------------------------+----------+-----------+------------------+
| HDU | EXTNAME                 | HDU Type | Data Type | Dimensions       |
+=====+=========================+==========+===========+==================+
|  0  | N/A                     | primary  | N/A       | N/A              |
+-----+-------------------------+----------+-----------+------------------+
|  1  | SCI                     | IMAGE    | float32   | 2024 x 2048 x 1  |
+-----+-------------------------+----------+-----------+------------------+
|  2  | ERR                     | IMAGE    | float32   | 2024 x 2048 x 1  |
+-----+-------------------------+----------+-----------+------------------+
|  3  | DQ                      | IMAGE    | uint32    | 2024 x 2048      |
+-----+-------------------------+----------+-----------+------------------+
|  4  | Flight Reference Stars  | BINTABLE | N/A       | 4 cols x nstars  |
+-----+-------------------------+----------+-----------+------------------+
|  5  | Planned Reference Stars | BINTABLE | N/A       | 10 cols x nstars |
+-----+-------------------------+----------+-----------+------------------+
|  6  | ASDF                    | BINTABLE | N/A       | variable         |
+-----+-------------------------+----------+-----------+------------------+

 - SCI: 3-D data array containing the pixel values, in units of DN/s. The data for the 2 integrations
   has been combined into a single image, as is done by the on-board FSW, resulting in a data array
   with NAXIS3 = 1.
 - ERR: 3-D data array containing uncertainty estimates for each pixel.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - Flight Reference Stars: A table containing information on the actual reference stars
   used by the FSW. Detailed contents are listed :ref:`below <flight_ref_stars>`.
 - Planned Reference Stars: A table containing information on the planned reference stars.
   Detailed contents are listed :ref:`below <planned_ref_stars>`.
 - ADSF: The data model meta data.

The FITS file structure for calibrated ID "stacked" products is as follows:

+-----+-------------------------+----------+-----------+------------------+
| HDU | EXTNAME                 | HDU Type | Data Type | Dimensions       |
+=====+=========================+==========+===========+==================+
|  0  | N/A                     | primary  | N/A       | N/A              |
+-----+-------------------------+----------+-----------+------------------+
|  1  | SCI                     | IMAGE    | float32   | 2304 x 2048 x 1  |
+-----+-------------------------+----------+-----------+------------------+
|  2  | ERR                     | IMAGE    | float32   | 2304 x 2048 x 1  |
+-----+-------------------------+----------+-----------+------------------+
|  3  | DQ                      | IMAGE    | uint32    | 2304 x 2048      |
+-----+-------------------------+----------+-----------+------------------+
|  4  | Flight Reference Stars  | BINTABLE | N/A       | 4 cols x nstars  |
+-----+-------------------------+----------+-----------+------------------+
|  5  | Planned Reference Stars | BINTABLE | N/A       | 10 cols x nstars |
+-----+-------------------------+----------+-----------+------------------+
|  6  | ASDF                    | BINTABLE | N/A       | variable         |
+-----+-------------------------+----------+-----------+------------------+

 - SCI: 3-D data array containing the pixel values, in units of DN/s. The data for the 2 integrations
   has been combined into a single image, as is done by the on-board FSW, resulting in a data array
   with NAXIS3=1.
 - ERR: 3-D data array containing uncertainty estimates for each pixel.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - Flight Reference Stars: A table containing information on the actual reference stars
   used by the FSW. Detailed contents are listed :ref:`below <flight_ref_stars>`.
 - Planned Reference Stars: A table containing information on the planned reference stars.
   Detailed contents are listed :ref:`below <planned_ref_stars>`.
 - ADSF: The data model meta data.

.. _flight_ref_stars:

Flight reference stars table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The structure and content of the Flight Reference Stars table is as follows.

+-------------------+-----------+-------------------------------+
| Column Name       | Data Type | Description                   |
+===================+===========+===============================+
| reference_star_id | char*2    | Reference star index          |
+-------------------+-----------+-------------------------------+
| id_x              | float64   | x position in FGS Ideal frame |
+-------------------+-----------+-------------------------------+
| id_y              | float64   | y position in FGS Ideal frame |
+-------------------+-----------+-------------------------------+
| count_rate        | float64   | count rate                    |
+-------------------+-----------+-------------------------------+

.. _planned_ref_stars:

Planned reference stars table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The structure and content of the Planned Reference Stars table is as follows.

+-------------------+-----------+-------------------------------+
| Column Name       | Data Type | Description                   |
+===================+===========+===============================+
| guide_star_order  | int32     | Guide star index within list  |
+-------------------+-----------+-------------------------------+
| reference_star_id | char*12   | GSC II identifier             |
+-------------------+-----------+-------------------------------+
| ra                | float64   | ICRS RA of the star           |
+-------------------+-----------+-------------------------------+
| dec               | float64   | ICRS Dec of the star          |
+-------------------+-----------+-------------------------------+
| id_x              | float64   | x position in FGS Ideal frame |
+-------------------+-----------+-------------------------------+
| id_y              | float64   | y position in FGS Ideal frame |
+-------------------+-----------+-------------------------------+
| fgs_mag           | float64   | magnitude                     |
+-------------------+-----------+-------------------------------+
| fgs_mag_uncert    | float64   | magnitude uncertainty         |
+-------------------+-----------+-------------------------------+
| count_rate        | float64   | count rate                    |
+-------------------+-----------+-------------------------------+
| count_rate_uncert | float64   | count rate uncertainty        |
+-------------------+-----------+-------------------------------+

ACQ1 mode
^^^^^^^^^
The "Acquistion" guiding function ACQ1 performs 128 x 128 pixel subarray readouts of the
detector, using 2 groups per integration and a total of 6 integrations.

The file name syntax for ACQ1 mode products is as follows::

 jw<pppppooovvv>_gs-acq1_<yyyydddhhmmss>-<uncal | cal>.fits

where:

 - pppppooovvv: the visit id, consisting of program id, obervation number, and visit number
 - gs-acq1: guide star acquisition1 product
 - yyyydddhhmmss: the time stamp at the *end* of the data contained in the file
 - uncal: indicates the uncalibrated (raw) product
 - cal: indicates the calibrated product

The FITS file structure for ACQ1 uncalibrated products is as follows:

+-----+---------+----------+-----------+-------------------+
| HDU | EXTNAME | HDU Type | Data Type | Dimensions        |
+=====+=========+==========+===========+===================+
|  0  | N/A     | primary  | N/A       | N/A               |
+-----+---------+----------+-----------+-------------------+
|  1  | SCI     | IMAGE    | uint16    | 128 x 128 x 2 x 6 |
+-----+---------+----------+-----------+-------------------+

 - SCI: 4-D data array containing the raw pixel values.

The FITS file structure for ACQ1 calibrated products is as follows:

+-----+---------+----------+-----------+---------------+
| HDU | EXTNAME | HDU Type | Data Type | Dimensions    |
+=====+=========+==========+===========+===============+
|  0  | N/A     | primary  | N/A       | N/A           |
+-----+---------+----------+-----------+---------------+
|  1  | SCI     | IMAGE    | float32   | 128 x 128 x 6 |
+-----+---------+----------+-----------+---------------+
|  2  | ERR     | IMAGE    | float32   | 128 x 128 x 6 |
+-----+---------+----------+-----------+---------------+
|  3  | DQ      | IMAGE    | uint32    | 128 x 128     |
+-----+---------+----------+-----------+---------------+
|  4  | ASDF    | BINTABLE | N/A       | variable      |
+-----+---------+----------+-----------+---------------+

 - SCI: 3-D data array containing the pixel values, in units of DN/s. Count rate images have been
   computed for each of the 6 integrations by differencing the 2 groups of each integration.
 - ERR: 3-D data array containing uncertainty estimates for each pixel.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - ADSF: The data model meta data.

ACQ2 mode
^^^^^^^^^
The "Acquisition" guiding function ACQ2 performs 32 x 32 pixel subarray readouts of the detector,
using 2 groups per integration and a total of 5 integrations.

The file name syntax for ACQ2 mode products is as follows::

 jw<pppppooovvv>_gs-acq2_<yyyydddhhmmss>-<uncal | cal>.fits

where:

 - pppppooovvv: the visit id, consisting of program id, obervation number, and visit number
 - gs-acq2: guide star acquisition2 product
 - yyyydddhhmmss: the time stamp at the *end* of the data contained in the file
 - uncal: indicates the uncalibrated (raw) product
 - cal: indicates the calibrated product

The FITS file structure for ACQ2 uncalibrated products is as follows:

+-----+---------+----------+-----------+-----------------+
| HDU | EXTNAME | HDU Type | Data Type | Dimensions      |
+=====+=========+==========+===========+=================+
|  0  | N/A     | primary  | N/A       | N/A             |
+-----+---------+----------+-----------+-----------------+
|  1  | SCI     | IMAGE    | uint16    | 32 x 32 x 2 x 5 |
+-----+---------+----------+-----------+-----------------+

 - SCI: 4-D data array containing the raw pixel values.

The FITS file structure for ACQ2 calibrated products is as follows:

+-----+---------+----------+-----------+-------------+
| HDU | EXTNAME | HDU Type | Data Type | Dimensions  |
+=====+=========+==========+===========+=============+
|  0  | N/A     | primary  | N/A       | N/A         |
+-----+---------+----------+-----------+-------------+
|  1  | SCI     | IMAGE    | float32   | 32 x 32 x 5 |
+-----+---------+----------+-----------+-------------+
|  2  | ERR     | IMAGE    | float32   | 32 x 32 x 5 |
+-----+---------+----------+-----------+-------------+
|  3  | DQ      | IMAGE    | uint32    | 32 x 32     |
+-----+---------+----------+-----------+-------------+
|  4  | ASDF    | BINTABLE | N/A       | variable    |
+-----+---------+----------+-----------+-------------+

 - SCI: 3-D data array containing the pixel values, in units of DN/s. Count rate images have been
   computed for each of the 5 integrations by differencing the 2 groups of each integration.
 - ERR: 3-D data array containing uncertainty estimates for each pixel.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - ADSF: The data model meta data.

Track mode
^^^^^^^^^^
The "Track" guiding function performs 32 x 32 pixel subarray readouts, the location of which
can move on the detector as the FGS FSW tracks the position of the guide star. The subarray
readouts are performed with a cadence of 16 Hz. Each integration consists of 2 groups, and the
total number of integrations (NINTS) can be very large (in the thousands).

The file name syntax for TRACK mode products is as follows::

 jw<pppppooovvv>_gs-track_<yyyydddhhmmss>-<uncal | cal>.fits

where:

 - pppppooovvv: the visit id, consisting of program id, obervation number, and visit number
 - gs-track: guide star track product
 - yyyydddhhmmss: the time stamp at the *end* of the data contained in the file
 - uncal: indicates the uncalibrated (raw) product
 - cal: indicates the calibrated product

The FITS file structure for TRACK uncalibrated products is as follows:

+-----+----------------------+----------+-----------+---------------------+
| HDU | EXTNAME              | HDU Type | Data Type | Dimensions          |
+=====+======================+==========+===========+=====================+
|  0  | N/A                  | primary  | N/A       | N/A                 |
+-----+----------------------+----------+-----------+---------------------+
|  1  | SCI                  | IMAGE    | uint16    | 32 x 32 x 2 x nints |
+-----+----------------------+----------+-----------+---------------------+
|  2  | Pointing             | BINTABLE | N/A       | 12 cols x nrows     |
+-----+----------------------+----------+-----------+---------------------+
|  3  | FGS Centroid Packet  | BINTABLE | N/A       | 17 cols x nrows     |
+-----+----------------------+----------+-----------+---------------------+
|  4  | Track subarray table | BINTABLE | N/A       | 5 cols x nrows      |
+-----+----------------------+----------+-----------+---------------------+

 - SCI: 4-D data array containing the raw pixel values.
 - Pointing: A table containing guide star position and jitter information.
   See :ref:`below <pointing_table>` for details of the contents.
 - FGS Centroid Packet: A table containing guide star centroiding information.
   See :ref:`below <centroid_table>` for details of the contents.
 - Track subarray table: A table containing subarray information over the duration of the product.
   See :ref:`below <subarray_table>` for details of the contents.

The FITS file structure for TRACK calibrated products is as follows:

+-----+----------------------+----------+-----------+-----------------+
| HDU | EXTNAME              | HDU Type | Data Type | Dimensions      |
+=====+======================+==========+===========+=================+
|  0  | N/A                  | primary  | N/A       | N/A             |
+-----+----------------------+----------+-----------+-----------------+
|  1  | SCI                  | IMAGE    | float32   | 32 x 32 x nints |
+-----+----------------------+----------+-----------+-----------------+
|  2  | ERR                  | IMAGE    | float32   | 32 x 32 x nints |
+-----+----------------------+----------+-----------+-----------------+
|  3  | DQ                   | IMAGE    | uint32    | 32 x 32         |
+-----+----------------------+----------+-----------+-----------------+
|  4  | POINTING             | BINTABLE | N/A       | 12 cols x nrows |
+-----+----------------------+----------+-----------+-----------------+
|  5  | FGS CENTROID PACKET  | BINTABLE | N/A       | 17 cols x nrows |
+-----+----------------------+----------+-----------+-----------------+
|  6  | TRACK SUBARRAY TABLE | BINTABLE | N/A       | 5 cols x nrows  |
+-----+----------------------+----------+-----------+-----------------+
|  7  | ASDF                 | BINTABLE | N/A       | variable        |
+-----+----------------------+----------+-----------+-----------------+

 - SCI: 3-D data array containing the pixel values, in units of DN/s. Count rate images for each
   integration have been computed by differencing the 2 groups in each integration.
 - ERR: 3-D data array containing uncertainty estimates for each pixel.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - Pointing: A table containing guide star position and jitter information.
   See :ref:`below <pointing_table>` for details of the contents.
 - FGS Centroid Packet: A table containing guide star centroiding information.
   See :ref:`below <centroid_table>` for details of the contents.
 - Track subarray table: A table containing subarray information over the duration of the product.
   See :ref:`below <subarray_table>` for details of the contents.
 - ADSF: The data model meta data.

.. _pointing_table:

Pointing table
~~~~~~~~~~~~~~
The structure and content of the Pointing table is as follows.

+-------------------+-----------+--------------+----------------------------------------------------+
| Column Name       | Data Type | Units        | Description                                        |
+===================+===========+==============+====================================================+
| time              | float64   | milli-sec    | Time since start of data file                      |
+-------------------+-----------+--------------+----------------------------------------------------+
| jitter            | float64   | milli-arcsec | :math:`sqrt(delta\_ddc\_ra^2 + delta\_ddc\_dec^2)` |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_ddc_ra      | float64   | milli-arcsec | Initial DDC RA - Current                           |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_ddc_dec     | float64   | milli-arcsec | Initial DDC Dec - Current                          |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_aperture_pa | float64   | milli-arcsec | Initial PA - Current                               |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_v1_ra       | float64   | milli-arcsec | Initial V frame RA - Current                       |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_v1_dec      | float64   | milli-arcsec | Initial V frame Dec - Current                      |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_v3_pa       | float64   | milli-arcsec | Initial V frame PA - Current                       |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_j1_ra       | float64   | milli-arcsec | Initial J frame RA - Current                       |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_j1_dec      | float64   | milli-arcsec | Initial J frame Dec - Current                      |
+-------------------+-----------+--------------+----------------------------------------------------+
| delta_j3_pa       | float64   | milli-arcsec | Initial J frame PA - Current                       |
+-------------------+-----------+--------------+----------------------------------------------------+
| HGA_motion        | int32     | N/A          | | HGA state: 0 = moving,                           |
|                   |           |              | | 1 = finished, 2 = offline                        |
+-------------------+-----------+--------------+----------------------------------------------------+

.. _centroid_table:

FGS Centroid Packet table
~~~~~~~~~~~~~~~~~~~~~~~~~
The structure and content of the Centroid Packet table is as follows.

+------------------------------------------------+-----------+----------------------------------------------+
| Column Name                                    | Data Type | Description                                  |
+================================================+===========+==============================================+
| observatory_time                               | char*23   | UTC time when packet was generated           |
+------------------------------------------------+-----------+----------------------------------------------+
| centroid_time                                  | char*23   | Fine guidance centroid time                  |
+------------------------------------------------+-----------+----------------------------------------------+
| guide_star_position_x                          | float64   | FGS Ideal Frame (arcsec)                     |
+------------------------------------------------+-----------+----------------------------------------------+
| guide_star_position_y                          | float64   | FGS Ideal Frame (arcsec)                     |
+------------------------------------------------+-----------+----------------------------------------------+
| guide_star_instrument_counts_per_sec           | float64   | Instrument counts/sec                        |
+------------------------------------------------+-----------+----------------------------------------------+
| signal_to_noise_current_frame                  | float64   | For current image frame                      |
+------------------------------------------------+-----------+----------------------------------------------+
| delta_signal                                   | float64   | Between current and previous frame           |
+------------------------------------------------+-----------+----------------------------------------------+
| delta_noise                                    | float64   | Between current and previous frame           |
+------------------------------------------------+-----------+----------------------------------------------+
| psf_width_x                                    | int32     | Bias from ideal guide star position (pixels) |
+------------------------------------------------+-----------+----------------------------------------------+
| psf_width_y                                    | int32     | Bias from ideal guide star position (pixels) |
+------------------------------------------------+-----------+----------------------------------------------+
| data_quality                                   | int32     | Centroid data quality                        |
+------------------------------------------------+-----------+----------------------------------------------+
| bad_pixel_flag                                 | char*4    | Bad pixel status for current subwindow (0/1) |
+------------------------------------------------+-----------+----------------------------------------------+
| bad_centroid_dq_flag                           | char*50   | Bad centroid for current subwindow (0/1)     |
+------------------------------------------------+-----------+----------------------------------------------+
| cosmic_ray_hit_flag                            | char*5    | NO/YES                                       |
+------------------------------------------------+-----------+----------------------------------------------+
| sw_subwindow_loc_change_flag                   | char*5    | NO/YES                                       |
+------------------------------------------------+-----------+----------------------------------------------+
| guide_star_at_detector_subwindow_boundary_flag | char*5    | NO/YES                                       |
+------------------------------------------------+-----------+----------------------------------------------+
| subwindow_out_of_FOV_flag                      | char*5    | NO/YES                                       |
+------------------------------------------------+-----------+----------------------------------------------+

.. _subarray_table:

Track Subarray table
~~~~~~~~~~~~~~~~~~~~
The Track Subarray table contains location and size information for the detector subarray window
that is used during the track function to follow the guide star.
The structure and content of the Track Subarray table is as follows.

+------------------+-----------+------------------------------------+
| Column Name      | Data Type | Description                        |
+==================+===========+====================================+
| observatory_time | char*23   | UTC time when packet was generated |
+------------------+-----------+------------------------------------+
| x_corner         | float64   | Subarray x corner (pixels)         |
+------------------+-----------+------------------------------------+
| y_corner         | float64   | Subarray y corner (pixels)         |
+------------------+-----------+------------------------------------+
| x_size           | int16     | Subarray x size (pixels)           |
+------------------+-----------+------------------------------------+
| y_size           | int16     | Subarray y size (pixels)           |
+------------------+-----------+------------------------------------+

FineGuide mode
^^^^^^^^^^^^^^
The "FineGuide" guiding function performs 8 x 8 pixel subarray readouts, at a fixed location
on the detector, and with a cadence of 16 Hz, from which the FGS FSW computes centroids for the
guide star. To reduce readout noise contribution to the centroid calculation, "Fowler" sampling
of the readouts is employed. Each integration consists of 4 readouts at the beginning, a
signal accumulation period, and 4 readouts at the end. The detector is then reset and the
readout cycle repeats for the next integration. The 4 readouts at the beginning are averaged
together, the 4 readouts at the end are averaged together, and then the difference of the 2
averages is computed to form a final countrate image for each integration. This approach to
creating the countrate images is used both on-board and in the :ref:`calwebb_guider <calwebb_guider>`
pipeline when the raw data are processed on the ground.

The file name syntax for FineGuide mode products is as follows::

 jw<pppppooovvv>_gs-fg_<yyyydddhhmmss>-<uncal | cal>.fits

where:

 - pppppooovvv: the visit id, consisting of program id, obervation number, and visit number
 - gs-fg: guide star FineGuide product
 - yyyydddhhmmss: the time stamp at the *end* of the data contained in the file
 - uncal: indicates the uncalibrated (raw) product
 - cal: indicates the calibrated product

The FITS file structure for FineGuide uncalibrated products is as follows:

+-----+----------------------+----------+-----------+-------------------+
| HDU | EXTNAME              | HDU Type | Data Type | Dimensions        |
+=====+======================+==========+===========+===================+
|  0  | N/A                  | primary  | N/A       | N/A               |
+-----+----------------------+----------+-----------+-------------------+
|  1  | SCI                  | IMAGE    | uint16    | 8 x 8 x 8 x nints |
+-----+----------------------+----------+-----------+-------------------+
|  2  | Pointing             | BINTABLE | N/A       | 12 cols x nrows   |
+-----+----------------------+----------+-----------+-------------------+
|  3  | FGS Centroid Packet  | BINTABLE | N/A       | 17 cols x nrows   |
+-----+----------------------+----------+-----------+-------------------+

 - SCI: 4-D data array containing the raw pixel values.
 - Pointing: A table containing guide star position and jitter information.
   See :ref:`above <pointing_table>` for details of the contents.
 - FGS Centroid Packet: A table containing guide star centroiding information.
   See :ref:`above <centroid_table>` for details of the contents.

The FITS file structure for FineGuide calibrated products is as follows:

+-----+----------------------+----------+-----------+-----------------+
| HDU | EXTNAME              | HDU Type | Data Type | Dimensions      |
+=====+======================+==========+===========+=================+
|  0  | N/A                  | primary  | N/A       | N/A             |
+-----+----------------------+----------+-----------+-----------------+
|  1  | SCI                  | IMAGE    | float32   | 8 x 8 x nints   |
+-----+----------------------+----------+-----------+-----------------+
|  2  | ERR                  | IMAGE    | float32   | 8 x 8 x nints   |
+-----+----------------------+----------+-----------+-----------------+
|  3  | DQ                   | IMAGE    | uint32    | 8 x 8           |
+-----+----------------------+----------+-----------+-----------------+
|  4  | POINTING             | BINTABLE | N/A       | 12 cols x nrows |
+-----+----------------------+----------+-----------+-----------------+
|  5  | FGS CENTROID PACKET  | BINTABLE | N/A       | 17 cols x nrows |
+-----+----------------------+----------+-----------+-----------------+
|  6  | ASDF                 | BINTABLE | N/A       | variable        |
+-----+----------------------+----------+-----------+-----------------+

 - SCI: 3-D data array containing the pixel values, in units of DN/s. Count rate images for each
   integration have been computed using the Fowler sampling scheme described above.
 - ERR: 3-D data array containing uncertainty estimates for each pixel.
 - DQ: 2-D data array containing DQ flags for each pixel.
 - Pointing: A table containing guide star position and jitter information.
   See :ref:`above <pointing_table>` for details of the contents.
 - FGS Centroid Packet: A table containing guide star centroiding information.
   See :ref:`above <centroid_table>` for details of the contents.
 - ADSF: The data model meta data.

