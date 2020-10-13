Description
===========

This step creates a catalog of source photometry and morphologies.
Both aperture and isophotal (segment-based) photometry are calculated.
Source morphologies are based on 2D image moments within the source
segment.


Source Detection
^^^^^^^^^^^^^^^^

Sources are detected using `image segmentation
<https://en.wikipedia.org/wiki/Image_segmentation>`_, which is a
process of assigning a label to every pixel in an image such that
pixels with the same label are part of the same source.  The
segmentation procedure used is from `Photutils source extraction
<https://photutils.readthedocs.io/en/latest/segmentation.html>`_.
Detected sources must have a minimum number of connected pixels that
are each greater than a specified threshold value in an image.  The
threshold level is usually defined at some multiple of the background
standard deviation above the background.  The image can also be
filtered before thresholding to smooth the noise and maximize the
detectability of objects with a shape similar to the filter kernel.

Source Deblending
^^^^^^^^^^^^^^^^^

Overlapping sources are detected as single sources.  Separating those
sources requires a deblending procedure, such as a multi-thresholding
technique used by `SExtractor
<https://www.astromatic.net/software/sextractor>`_.  Here we use the
`Photutils deblender
<https://photutils.readthedocs.io/en/latest/segmentation.html#source-deblending>`_,
which is an algorithm that deblends sources using a combination of
multi-thresholding and `watershed segmentation
<https://en.wikipedia.org/wiki/Watershed_(image_processing)>`_.  In
order to deblend sources, they must be separated enough such that
there is a saddle between them.

Source Photometry and Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After detecting sources using image segmentation, we can measure their
photometry, centroids, and morphological properties.  The aperture
photometry is measured in three apertures, based on the input
encircled energy values.  The total aperture-corrected flux and
magnitudes are also calculated, based on the largest aperture.  Both
AB and Vega magnitudes are calculated.

The isophotal photometry is based on `photutils segmentation
<https://photutils.readthedocs.org/en/latest/segmentation.html>`_.
The properties that are currently calculated for each source include
source centroids (both in pixel and sky coordinates), isophotal fluxes
(and errors), AB and Vega magnitudes (and errors), isophotal area,
semimajor and semiminor axis lengths, orientation of the major axis,
and sky coordinates at corners of the minimal bounding box enclosing
the source.

.. Note::

   Errors are only created when an image has an error extension.  Products
   created from the resampling step currently do not have an error extension 
   and the error columns are currently filled with a value of nan. 

Source Catalog Table
^^^^^^^^^^^^^^^^^^^^

The output source catalog table is saved in `ECSV format
<https://docs.astropy.org/en/stable/io/ascii/write.html#ecsv-format>`_.

The table contains a row for each source, with the following columns:

+------------------------+----------------------------------------------------+
| Column                 | Description                                        |
+========================+====================================================+
| id                     | Unique source identification number                |
+------------------------+----------------------------------------------------+
| xcentroid              | X pixel value of the source centroid               |
+------------------------+----------------------------------------------------+
| ycentroid              | Y pixel value of the source centroid               |
+------------------------+----------------------------------------------------+
| sky_centroid           | Sky coordinate of the source centroid              |
+------------------------+----------------------------------------------------+
| aper_bkg_flux          | The local background value calculated as the       |
|                        | sigma-clipped median value in the background       |
|                        | annulus aperture                                   |
+------------------------+----------------------------------------------------+
| aper_bkg_flux_err      | The standard error of the sigma-clipped median     |
|                        | background value                                   |
+------------------------+----------------------------------------------------+
| aper30_flux            | Flux within the 30% encircled energy circular      |
|                        | aperture                                           |
+------------------------+----------------------------------------------------+
| aper30_flux_err        | Flux error within the 30% encircled energy         |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper50_flux            | Flux within the 50% encircled energy circular      |
|                        | aperture                                           |
+------------------------+----------------------------------------------------+
| aper50_flux_err        | Flux error within the 50% encircled energy         |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper70_flux            | Flux within the 70% encircled energy circular      |
|                        | aperture                                           |
+------------------------+----------------------------------------------------+
| aper70_flux_err        | Flux error within the 70% encircled energy         |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper_total_flux        | Total aperture-corrected flux based on the 70%     |
|                        | encircled energy circular aperture; calculated     |
|                        | only for stars                                     |
+------------------------+----------------------------------------------------+
| aper_total_flux_err    | Total aperture-corrected flux error based on the   |
|                        | 70% encircled energy circular aperture; calculated |
|                        | only for stars                                     |
+------------------------+----------------------------------------------------+
| aper30_abmag           | AB magnitude within the 30% encircled energy       |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper30_abmag_err       | AB magnitude error within the 30% encircled energy |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper50_abmag           | AB magnitude within the 50% encircled energy       |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper50_abmag_err       | AB magnitude error within the 50% encircled energy |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper70_abmag           | AB magnitude within the 70% encircled energy       |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper70_abmag_err       | AB magnitude error within the 70% encircled energy |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper_total_abmag       | Total aperture-corrected AB magnitude based on the |
|                        | 70% encircled energy circular aperture; calculated |
|                        | only for stars                                     |
+------------------------+----------------------------------------------------+
| aper_total_abmag_err   | Total aperture-corrected AB magnitude error based  |
|                        | on the 70% encircled energy circular aperture;     |
|                        | calculated only for stars                          |
+------------------------+----------------------------------------------------+
| aper30_vegamag         | Vega magnitude within the 30% encircled energy     |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper30_vegamag_err     | Vega magnitude error within the 30% encircled      |
|                        | energy circular aperture                           |
+------------------------+----------------------------------------------------+
| aper50_vegamag         | Vega magnitude within the 50% encircled energy     |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper50_vegamag_err     | Vega magnitude error within the 50% encircled      |
|                        | energy circular aperture                           |
+------------------------+----------------------------------------------------+
| aper70_vegamag         | Vega magnitude within the 70% encircled energy     |
|                        | circular aperture                                  |
+------------------------+----------------------------------------------------+
| aper70_vegamag_err     | Vega magnitude error within the 70% encircled      |
|                        | energy circular aperture                           |
+------------------------+----------------------------------------------------+
| aper_total_vegamag     | Total aperture-corrected Vega magnitude based on   |
|                        | the 70% encircled energy circular aperture;        |
|                        | calculated only for stars                          |
+------------------------+----------------------------------------------------+
| aper_total_vegamag_err | Total aperture-corrected Vega magnitude error      |
|                        | based on the 70% encircled energy circular         |
|                        | aperture; calculated only for stars                |
+------------------------+----------------------------------------------------+
| CI_30_50               | Concentration index calculated as aper30_abmag -   |
|                        | aper50_abmag                                       |
+------------------------+----------------------------------------------------+
| CI_50_70               | Concentration index calculated as aper50_abmag -   |
|                        | aper70_abmag                                       |
+------------------------+----------------------------------------------------+
| CI_30_70               | Concentration index calculated as aper30_abmag -   |
|                        | aper70_abmag                                       |
+------------------------+----------------------------------------------------+
| is_star                | Flag indicating whether the source is a star       |
+------------------------+----------------------------------------------------+
| sharpness              | The DAOFind source sharpness statistic             |
+------------------------+----------------------------------------------------+
| roundness              | The DAOFind source roundness statistic             |
+------------------------+----------------------------------------------------+
| nn_dist                | The distance in pixels to the nearest neighbor     |
+------------------------+----------------------------------------------------+
| nn_abmag               | The AB magnitude of the nearest neighbor.  If the  |
|                        | object is a star it is the total aperture-         |
|                        | corrected AB magnitude, otherwise it is the        |
|                        | isophotal AB magnitude.                            |
+------------------------+----------------------------------------------------+
| isophotal_flux         | Isophotal flux                                     |
+------------------------+----------------------------------------------------+
| isophotal_flux_err     | Isophotal flux error                               |
+------------------------+----------------------------------------------------+
| isophotal_abmag        | Isophotal AB magnitude                             |
+------------------------+----------------------------------------------------+
| isophotal_abmag_err    | Isophotal AB magnitude error                       |
+------------------------+----------------------------------------------------+
| isophotal_vegamag      | Isophotal Vega magnitude                           |
+------------------------+----------------------------------------------------+
| isophotal_vegamag_err  | Isophotal Vega magnitude error                     |
+------------------------+----------------------------------------------------+
| isophotal_area         | Isophotal area                                     |
+------------------------+----------------------------------------------------+
| semimajor_sigma        | 1-sigma standard deviation along the semimajor     |
|                        | axis of the 2D Gaussian function that has the same |
|                        | second-order central moments as the source         |
+------------------------+----------------------------------------------------+
| semiminor_sigma        | 1-sigma standard deviation along the semiminor     |
|                        | axis of the 2D Gaussian function that has the same |
|                        | second-order central moments as the source         |
+------------------------+----------------------------------------------------+
| ellipticity            | 1 minus the ratio of the 1-sigma lengths of the    |
|                        | semimajor and semiminor axes                       |
+------------------------+----------------------------------------------------+
| orientation            | The angle (degrees) between the positive X axis    |
|                        | and the major axis (increases counter-clockwise)   |
+------------------------+----------------------------------------------------+
| sky_orientation        | The position angle (degrees) from North of the     |
|                        | major axis                                         |
+------------------------+----------------------------------------------------+
| sky_bbox_ll            | Sky coordinate of the lower-left vertex of the     |
|                        | minimal bounding box of the source                 |
+------------------------+----------------------------------------------------+
| sky_bbox_ul            | Sky coordinate of the upper-left vertex of the     |
|                        | minimal bounding box of the source                 |
+------------------------+----------------------------------------------------+
| sky_bbox_lr            | Sky coordinate of the lower-right vertex of the    |
|                        | minimal bounding box of the source                 |
+------------------------+----------------------------------------------------+
| sky_bbox_ur            | Sky coordinate of the upper-right vertex of the    |
|                        | minimal bounding box of the source                 |
+------------------------+----------------------------------------------------+
