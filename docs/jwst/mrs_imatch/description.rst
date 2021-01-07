.. _mrs_imatch-description-label:

Description
============

Overview
--------
The ``mrs_imatch`` step "matches" image intensities of several input
2D MIRI MRS images by fitting polynomials to cube intensities (cubes built
from the input 2D images), in such a way as to minimize - in the least squares
sense - inter-image mismatches in intensity. The "background matching" polynomials
are defined in the frame of world coordinates (e.g. RA, DEC, lambda).


If any of background polynomial coefficients are a nan then the step is skipped and
S_MRSMAT is set to SKIPPED.

Any sources in the scene are identified via sigma clipping and removed from the
matching region.

Assumptions
-----------
Because the fitted polynomials are defined in terms of world coordinates, and because
the algorithm needs to build 3D cubes for each input image, all input images need
to have a valid WCS defined.

Algorithm
---------
This step builds a system of linear equations

.. math::
    a \cdot c = b

whose solution :math:`c` is a set of coefficients of (multivariate)
polynomials that represent the "background" in each input image (these are
polynomials that are "corrections" to the intensities in the input images) such
that the following sum is minimized:

.. math::
    L = \sum^N_{n,m=1,n \neq m} \sum_k \frac{\left[I_n(k) - I_m(k) - P_n(k) + P_m(k)\right]^2}{\sigma^2_n(k) + \sigma^2_m(k)}.

In the above equation, index :math:`k=(k_1,k_2,...)` labels a position
in an input image's pixel grid [NOTE: all input images share a common
pixel grid].

"Background" polynomials :math:`P_n(k)` are defined through the
corresponding coefficients as:

.. math::
    P_n(k_1,k_2,...) = \sum_{d_1=0,d_2=0,...}^{D_1,D_2,...} c_{d_1,d_2,...}^n \cdot k_1^{d_1} \cdot k_2^{d_2}  \cdot \ldots .

Step Arguments
==============
The ``mrs_imatch`` step has two optional arguments:

``bkg_degree``
  The background polynomial degree (int; default=1)

``subtract``
  Indicates whether the computed matching "backgrounds" should be subtracted
  from the image data (bool; default=False)

Reference Files
===============
This step does not require any reference files.

Also See
========
See :doc:`wiimatch package documentation <../wiimatch/index>` for more details.
