Description
===========

:Class: `jwst.picture_frame.PictureFrameStep`
:Alias: picture_frame

Overview
--------
The ``picture_frame`` step removes thermal artifacts (the "picture frame"
effect) from calibrated ramp images, after the :ref:`jump <jump_step>` step
and prior to performing the :ref:`clean_flicker_noise <clean_flicker_noise_step>`
or the :ref:`ramp_fitting <ramp_fitting_step>` step.

The correction is available for NIRSpec FULL frame exposures only.

The picture frame artifacts are corrected by scaling and subtracting a reference rate
image, containing only the thermal background, stored in a PICTUREFRAME reference file.

Correction Algorithm
--------------------
The ``picture_frame`` step proceeds similarly to the
:ref:`clean_flicker_noise <clean_flicker_noise_step>` step: it generates a draft
rate image from the input ramp, creates a background mask from the rate image,
fits the correction to the background data in each group in the ramp, and subtracts the
correction from the group image.

To create the background mask, the step
uses the process outlined in :ref:`scene_mask`, with the difference that science
regions are masked by default for picture frame processing.

After background data is identified, median values are computed over the
edges (:math:`m_{eref}`) and center (:math:`m_{cref}`) of the reference picture
frame image, at the valid background locations only.  Then, for each group image
in each integration, the cleaning process is:

#. Make a difference image (current group – first group) to correct.  The first
   group is left uncorrected.

#. Calculate the median values at the edges (:math:`m_e`) and center (:math:`m_c`)
   of the difference image at valid background locations.

#. Use the median values for the difference image to scale and offset
   the reference picture frame rate image (:math:`p`):

    .. math::

       scale &= (m_c - m_e) / (m_{cref} - m_{eref})

       correction &= scale * (p - m_{eref}) + m_e

#. Subtract the correction image from the group image.


References
==========

The picture frame correction algorithm is based on work by M. Regan and
E. Bergeron (in prep).
