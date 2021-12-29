Introduction
------------
Subtraction of background signal can take several different forms depending on the
observing mode and the available data. Here we give an overview of the different
methods that are available, when they can be used, and where they occur in the
processing flow. Imaging and spectroscopic observations share one method for background
subtraction, while others are unique to spectroscopic data only. See the documentation
for the individual steps mentioned here for complete details on how each of them
function.

Imaging Mode
------------
Background subtraction for imaging data is currently available in several places
within the calibration pipeline stages.

1. Image-from-image subtraction can be performed by the
   :ref:`background <background_step>` step during :ref:`calwebb_image2 <calwebb_image2>`
   processing. The background images come from observations of a dedicated
   background target.
2. Background matching and subtraction can be performed within an ensemble of
   images by the :ref:`skymatch <skymatch_step>` step during
   :ref:`calwebb_image3 <calwebb_image3>` processing.
3. Local background subtraction for individual sources can be performed by the
   :ref:`source_catalog <source_catalog_step>` step within the
   :ref:`calwebb_image3 <calwebb_image3>` pipeline.

Spectroscopic Modes
-------------------
Spectroscopic observations allow for some additional ways of performing
background subtraction. The list of options includes:

1. Image-from-image subtraction can be performed by the
   :ref:`background <background_step>` step during :ref:`calwebb_spec2 <calwebb_spec2>`
   processing. The background images can come from:

   a) Observations of a dedicated background target
   b) Nodded observations of a point-like science target

2. Subtraction of a "master" background spectrum, where the master background
   spectrum can come from:

   a) Observations of a dedicated background target
   b) Nodded observations of a point-like science target
   c) Dedicated background slitlets in a NIRSpec MOS exposure
   d) A user-supplied spectrum

3. Local background subtraction for individual spectral can be performed by
   the :ref:`extract_1d <extract_1d_step>` step when doing 1D spectral
   extraction.

The following table shows the list of image-from-image and master background
subtraction methods available for various spectroscopic observation modes, and
indicates the pipeline and step in which the subtraction operation occurs.
The table also shows which method is applied by default in the operational pipeline
when the available data support multiple methods.

.. Note:: Master background subtraction is applied in the
          :ref:`calwebb_spec3 <calwebb_spec3>` pipeline for most spectroscopic modes,
          but for **NIRSpec MOS** mode it is applied during
          :ref:`calwebb_spec2 <calwebb_spec2>` processing.

+--------------------------+---------------+-------------------+-----------------------------+
|                          | calwebb_spec2 | calwebb_spec3     | calwebb_spec2               |
|                          |               |                   |                             |
| Mode                     | background    | master_background | master_background_nrs_slits |
+==========================+===============+===================+=============================+
| **NIRSpec Fixed Slit:**  |               |                   |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Dedicated background     | Default       | Optional          |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Nodded point source      | Default       | Optional          |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| User supplied            |               | Default           |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| **NIRSpec IFU:**         |               |                   |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Dedicated background     | Default       | Optional          |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Nodded point source      | Default       | Optional          |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| User supplied            |               | Default           |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| **NIRSpec MOS:**         |               |                   |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Background slitlets      |               |                   | Default                     |
+--------------------------+---------------+-------------------+-----------------------------+
| Nodded point source      | Default       |                   |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| User supplied            |               |                   | Default                     |
+--------------------------+---------------+-------------------+-----------------------------+
| **MIRI LRS Fixed Slit:** |               |                   |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Dedicated background     | Default       | Optional          |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Nodded point source      | Default       | Optional          |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| User supplied            |               | Default           |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| **MIRI MRS:**            |               |                   |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Dedicated background     | Default       | Optional          |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| Nodded point source      | Default       | Optional          |                             |
+--------------------------+---------------+-------------------+-----------------------------+
| User supplied            |               | Default           |                             |
+--------------------------+---------------+-------------------+-----------------------------+

These background subtraction methods are only available for the observing modes
listed in the table. Other spectroscopic modes, including NIRCam and NIRISS Wide Field
Slitless Spectroscopy (WFSS), NIRCam Time Series Grism, NIRISS Single Object Slitless
Spectroscopy (SOSS), and MIRI LRS slitless, use other ways of handling background.

Image-from-Image Subtraction
----------------------------
As explained in the documentation for the :ref:`background <background_step>` step,
this process combines one or more exposures to be treated as backgrounds into a
sigma-clipped mean background image, which is then directly subtracted, in
detector space, from an exposure being processed in the :ref:`calwebb_image2 <calwebb_image2>`
or :ref:`calwebb_spec2 <calwebb_spec2>` pipelines for imaging or spectroscopic
data, respectively. For imaging mode observations this is only possible when
observations of a designated background target have been obtained. For spectroscopic
modes this is possible either through observations of a designated background target
or when nodded exposures of a point-like target are obtained (e.g. using the MIRI LRS
"ALONG-SLIT-NOD" dither pattern for an LRS fixed slit observation). Exposures from
one nod position can be used as background for exposures at the other nod position,
assuming the source is point-like.

In either instance, the exposures to be used as background are included in the
`image2` or `spec2` ASN file used to process the science target exposures, where
the background exposures are labeled with an ASN member type of "background".

Spectroscopic observations that have designated background target exposures or
nodded exposures can use either the image-from-image or master background subtraction
methods. In the operational pipeline the image-from-image subtraction method is applied
by default and the master background subtraction is skipped. A user has the option to
reprocess the data and apply the other method, if desired.

Master Background Subtraction
-----------------------------
In general, the master background subtraction method works by taking a 1D
background spectrum, interpolating it back into the 2D space of a science image,
and then subtracting it. This allows for higher SNR background data to be used
than what might be obtainable by doing direct image-from-image subtraction using
only one or a few background images. The 1D master background spectrum can either
be constructed on-the-fly by the calibration pipeline from available background
data or supplied by the user. See the documentation for the
:ref:`master background subtraction <master_background_step>` step for full details.

As with image-from-image subtraction, there are different ways of obtaining the
data necessary for constructing a master background spectrum, depending on the
observing mode:

1. Observations of a designated background target
2. Nodded observations of a point-like source
3. Dedicated background slitlets in a NIRSpec MOS exposure
4. User-supplied master background spectrum

All of these scenarios apply the master background subtraction during
:ref:`calwebb_spec3 <calwebb_spec3>` processing, except for NIRSpec MOS observations.
Master background subtraction for NIRSpec MOS, using either data from background
slitlets contained in each MOS exposure or a user-supplied master background spectrum,
is applied during :ref:`calwebb_spec2 <calwebb_spec2>`, due to unique methods that
must be used for MOS exposures.

For scenarios that apply master background subtraction during
:ref:`calwebb_spec3 <calwebb_spec3>` processing, the fully-calibrated 1D spectra
("x1d" products) from either dedicated background target exposures or nodded
science exposures are used by the :ref:`master_background <master_background_step>`
step to construct the 1D master background spectrum. These are the x1d products created
during the last step of the preceding :ref:`calwebb_spec2 <calwebb_spec2>` pipeline
when it is used to process each exposure. Again, see the documentation for the
:ref:`master background subtraction <master_background_step>` step for full
details of the source of the background data for these scenarios.

If the user supplies a 1D master background spectrum, the construction of the
master background spectrum in the pipeline is skipped and the user-supplied
spectrum is used in its place. This applies to all modes, including NIRSpec MOS.

As mentioned above, NIRSpec MOS observations require special handling to correctly
apply master background subtraction. If a MOS observation uses an MSA configuration
that includes one or more slitlets containing only background signal, the background
slitlets are fully calibrated and extracted to produce one or more 1D background
spectra. The background spectra are combined into a 1D master background spectrum,
which is then interpolated back into the 2D space of all slitlets and subtracted.
If the user supplies a master background spectrum for a MOS observation,
that spectrum is used to do the subtraction. Again note that for NIRSpec MOS mode
these operations take place during :ref:`calwebb_spec2 <calwebb_spec2>` pipeline
processing, not :ref:`calwebb_spec3 <calwebb_spec3>` like all other modes.
