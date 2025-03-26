Description
===========

:Class: `jwst.jump.JumpStep`
:Alias: jump

This step finds and flags outliers (usually caused by cosmic-ray hits) in
each pixel of an "up the ramp" IR exposure.

Assumptions
-----------
We assume that the :ref:`saturation <saturation_step>` step has already been applied to
the input exposure, so that saturated ramp groups are appropriately flagged in the
input GROUPDQ array. We also assume that steps such as
:ref:`reference pixel correction <refpix_step>` and
:ref:`non-linearity correction <linearity_step>` have been applied,
so that the input data ramps do not have any non-linearities or noise above the modeled Poisson
and read noise due to instrumental effects. The absence of any of these preceding corrections
or the presence of residual non-linearities and noise can lead to false detection of jumps
in the ramps, due to departure from linearity.

The ``jump`` step will automatically skip execution if the input data contain fewer
than 3 groups per integration, because the baseline algorithm requires at least
two first differences to work.

Note that the core algorithms for this step are called from the external package
``stcal``, an STScI effort to unify common calibration processing algorithms
for use by multiple observatories.

The keywords ``PRIMECRS`` and ``EXTNCRS`` are set in this step.  The ``PRIMECRS`` keyword
is the number of primary cosmic rays found per thousand pixels per second.  The ``EXTNCRS``
keyword is the number of extended events (snowball and shower) per million pixels per
second.


:ref:`Algorithm <stcal:jump_algorithm>`
---------------------------------------

Large Events (Snowballs and Showers)
------------------------------------
All the detectors on JWST are affected by large cosmic ray
events. While these events, in general, affect a large number of
pixels, the more distinguishing characteristic is that they are
surrounded by a halo of pixels that have a low level of excess
counts. These excess counts are, in general, below the detection
threshold of normal cosmic rays.

To constrain the effect of this halo, the jump step will fit ellipses or circles that
enclose the large events and expand the ellipses and circles by the input
`expand_factor` and mark them as jump (see :ref:`jump step arguments <jump_arguments>`
for details).

The two different types of JWST detectors respond differently. The large events in the near-infrared
detectors are almost always circles with a central region that is saturated.
The saturated core allows the search for smaller events without false positives.
The mid-IR (MIRI) detectors do not, in general, have a saturated center and are only rarely circular.
Thus, we fit the minimum enclosing ellipse and do not require that there are saturated pixels
within the ellipse.  Likewise, MIRI showers are only flagged when detected features are consistent
with the maximum known amplitude (in DN/s) of shower artifacts.

Multiprocessing
---------------
This step has the option of running in multiprocessing mode. In that mode it will
split the input data cube into a number of row slices based on the number of available
cores on the host computer and the value of the ``max_cores`` input parameter. By
default the step runs on a single processor. At the other extreme, if ``max_cores`` is
set to "all", it will use all available cores (real and virtual). Testing has shown
a reduction in the elapsed time for the step proportional to the number of real
cores used. Using the virtual cores also reduces the elapsed time, but at a slightly
lower rate than the real cores.

If multiprocessing is requested, the input cube will be divided into a number of
slices in the row dimension (with the last slice being slightly larger, if needed),
and sent for processing in parallel.
In the event the number of cores (and hence slices) selected exceeds the number of
available image rows, the number of slices will be reduced to match the number of rows.
After all the slices have finished processing, the output GROUPDQ cube - containing
the DQ flags for detected jumps - is reassembled from the slices.

Subarrays
---------
Full-frame reference files can be used for all science exposures even if the
science exposure was taken in a subarray mode. If so, subarrays will be
extracted from the reference file data to match the science exposure.
Alternatively, subarray-specific reference files, which match the science
exposure, may be used.
