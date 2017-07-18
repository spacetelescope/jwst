.. _mrs_imatch-description-label:

Description
============

Overview
--------
The ``mrs_imatch`` step that "matches" image intensities of several input
2D MIRI MRS images by fitting polynomials to cube intensities (cubes built
from input 2D images) in such a way as to minimize inter-image mismatch
in the least squares sense. These "background matching" polynomials
are defined in terms of world coordinates (e.g., ``RA``, ``DEC``, ``lambda``).

Assumptions
-----------
Because polynomials are defined in terms of world coordinates, and because
the algorithm needs to build 3D cubes for each input image, input images need
to have valid WCS.

Algorithm
---------
This step builds a system of linear equations

.. math::
    a \cdot c = b

whose solution :math:`c` is a set of coefficients of (multivariate)
polynomials that represent the "background" in each input image (these are
polynomials that are "corrections" to intensities of input images) such
that the following sum is minimized:

.. math::
    L = \sum^N_{n,m=1,n \neq m} \sum_k \frac{\left[I_n(k) - I_m(k) - P_n(k) + P_m(k)\right]^2}{\sigma^2_n(k) + \sigma^2_m(k)}.

In the above equation, index :math:`k=(k_1,k_2,...)` labels a position
in input image's pixel grid [NOTE: all input images share a common
pixel grid].

"Background" polynomials :math:`P_n(k)` are defined through the
corresponding coefficients as:

.. math::
    P_n(k_1,k_2,...) = \sum_{d_1=0,d_2=0,...}^{D_1,D_2,...} c_{d_1,d_2,...}^n \cdot k_1^{d_1} \cdot k_2^{d_2}  \cdot \ldots .

Step Arguments
==============
The mrs_imatch step has two optional argument:

* ``bkg_degree``: An integer background polynomial degree (Default: 1)

* ``subtract``: A boolean value indicating whether computed matching
  "backgrounds" should subtracted from image data (Default: `False`).

Reference Files
===============
This step does not require any reference files.

Also See
========
See :doc:`wiimatch package documentation <../wiimatch/index>` for more details.
