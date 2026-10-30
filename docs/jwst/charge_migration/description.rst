Description
===========

:Class: `jwst.charge_migration.charge_migration_step.ChargeMigrationStep`
:Alias: charge_migration


Overview
--------
This step addresses artifacts seen in undersampled JWST exposures due to
:ref:`charge migration <charge_migration>`: the spilling of charge from a saturated
pixel into its neighboring pixels. This charge migration causes the full width half
maximum of the point spread function (PSF) to be larger for bright point sources
compared to faint point sources, also known as the "brighter-fatter effect" (BFE).

This charge migration causes group-to-group differences
to decrease significantly once the signal level is greater than ~25,000 ADU. As a result,
the last several groups of these ramps get flagged by the :ref:`jump <jump_step>`
step.  The smaller number of groups used for these pixels in the :ref:`ramp_fitting <ramp_fitting_step>`
step results in them having  larger read noise variances, which in turn may lead to lower weights used
during resampling steps in later pipeline processing.

This step is primarily intended to address charge migration by flagging pixels with
signal over a specified threshold, along with their immediate neighbors, so that they
may be handled appropriately in later processing.  For this flagging algorithm, changes
are expected only to the GROUPDQ extension.

For correcting the effects of charge migration to recover signal from affected pixels,
this step also implements an experimental method to fit and correct for the excess or
missing flux on a per-pixel basis.  This algorithm is off by default for all modes,
but is available for testing.  If applied, it will directly modify the ramp values for
affected pixels, in the SCI extension.

The input must be in the form of a `~stdatamodels.jwst.datamodels.RampModel`.
The output will have the same format, with the updates to the internal data and group DQ
arrays applied.


Algorithm
---------

Flagging Charge Loss
^^^^^^^^^^^^^^^^^^^^

For each pixel in each integration, the first group, and all subsequent groups, exceeding
the value of the ``signal_threshold`` parameter is flagged as CHARGELOSS. ``signal_threshold``
is in units of ADUs. These groups will also be flagged as DO_NOT_USE, and will not
be included in the slope calculation during the ``ramp_fitting`` step. Despite being flagged
as DO_NOT_USE, these CHARGELOSS groups are still included in the calculation of the
variance due to readnoise.

This results in a readnoise variance for affected pixels that is similar to that of
pixels unaffected by charge migration. For the Poisson noise variance calculation in
:ref:`ramp_fitting <ramp_fitting_step>`, the CHARGELOSS/DO_NOT_USE groups are not included.

For integrations having only 1 or 2 groups, no flagging will be performed.

Correcting for Charge Migration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   This correction algorithm is currently an experimental prototype.
   It is not yet recommended for general use.

To correct for the brighter-fatter effect, this step fits a linear trend to the group
differences for any pixel with maximum signal greater than a specified threshold, along
with its neighbors.  For pixels affected by charge loss, the linear trend will have a negative slope,
as the rate decreases over time.  For the neighbors receiving excess charge, the linear
trend will have a positive slope. In either case, correcting the ramp values for the
fit slope will mitigate the effects of the charge migration, while tending to conserve
flux over the PSF.

The process is:

1. Identify pixels with maximum counts over the correction threshold, along with
   neighbors within a specified padding window (parameter ``mask_pad``).
2. Compute differences of all groups in the ramp, masking any bad pixels, and
   median-combine across all integrations.
3. If more than 5 diffs are available, reject outliers with sigma clipping.
4. For each identified pixel, fit a line to the clipped median diffs.
5. Reject the fit if it does not meet these criteria:

   a. The value change due to the slope over the range of the fit values is
      less than a specified factor times the read noise.  The factor is set by the
      ``snr_limit`` parameter.
   b. The fit value for the first group difference is positive.
   c. If there are more than 5 diffs, the RMS of the corrected data is reduced
      by at least a specified factor of the original RMS. The factor is set
      by the ``rms_limit`` parameter.

6. For pixels with valid corrections, subtract the fit slope from the group
   differences.
7. Reconstruct the ramps from the modified differences, adding each new diff
   into the previous corrected group.


References
----------

The flagging algorithm is based on work by P. Goudfrooij et. al.,
"An Algorithm to Mitigate Charge Migration Effects in Data from the Near
Infrared Imager and Slitless Spectrograph on the James Webb Space Telescope."
(`2024, PASP, 136, 4503 <https://ui.adsabs.harvard.edu/abs/2024PASP..136a4503G/abstract>`__).

The correction algorithm is based on work by D. Law, in prep.
