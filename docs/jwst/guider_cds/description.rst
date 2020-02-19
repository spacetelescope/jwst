Description
============

The ``guider_cds`` step computes countrate images from the
Correlated Double Sampling (CDS) detector readouts used in FGS
guiding mode data. The exact way in which the countrate images
are computed depends on the guiding mode (ID, ACQ1, ACQ2,
TRACK, FineGuide) in use.

ID mode
-------
The ID mode uses 2 integrations (NINTS=2) with 2 groups per
integration (NGROUPS=2). For this mode the ``guider_cds`` step
first computes a difference image for each integration by
subtracting group 1 from group 2. A final difference image is
then computed by taking the minimum value at each pixel from
the 2 integrations. The minimum difference image is then divided
by the group time to produce a countrate image. The output
data array is 3D, with dimensions of (ncols x nrows x 1).

ACQ1, ACQ2, and TRACK modes
---------------------------
These modes use multiple integrations (NINTS>1) with 2 groups
per integration (NGROUPS=2). For these modes the ``guider_cds``
step computes a countrate image for each integration, by
subtracting group 1 from group 2 and dividing by the group time.
The output data array is 3D, with dimensions of
(ncols x nrows x nints).

FineGuide mode
--------------
The FineGuide mode uses many integrations (NINTS>>1) with 4
groups at the beginning and 4 groups at the end of each
integration. The ``guider_cds`` step computes a countrate
image for each integration by subtracting the average of the
first 4 groups from the average of the last 4 groups and
dividing by the group time. The output data array is
3D, with dimensions of (ncols x nrows x nints).

After successful completion of the step, the "BUNIT" keyword in
the output data is updated to "DN/s" and the "S_GUICDS"
keyword is set to "COMPLETE".
