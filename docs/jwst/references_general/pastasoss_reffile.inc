.. _pastasoss_reffile:

PASTASOSS Reference File
------------------------

:REFTYPE: PASTASOSS
:Data model: `~stdatamodels.jwst.datamodels.PastasossModel`

The PASTASOSS (Predicting Accurate Spectral Traces in Astrophysical SOSS Spectra)
reference file contains the spectral trace positions and wavelength calibration models for
NIRISS SOSS (Single Object Slitless Spectroscopy) mode observations. It is used by the
:ref:`extract_1d <extract_1d_step>` step to predict the locations of spectral orders 1, 2,
and 3 on the detector, accounting for small rotations of the trace introduced by
visit-to-visit offsets of the GR700XD grism. It also provides the spectral throughput
for each order.

Reference Selection Keywords for PASTASOSS
+++++++++++++++++++++++++++++++++++++++++++
CRDS selects the appropriate PASTASOSS reference file based on the following keywords.
PASTASOSS is only applicable for NIRISS SOSS observations.
All keywords used for file selection are *required*.

========== =====================================================================
Instrument Keywords
========== =====================================================================
NIRISS     INSTRUME, DETECTOR, EXP_TYPE, FILTER, DATE-OBS, TIME-OBS
========== =====================================================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The PASTASOSS reference file is in ASDF format. It contains trace position arrays,
polynomial wavelength calibration models, and throughput curves for
each of the NIRISS SOSS spectral orders.

The top-level structure of the reference file is as follows:

:meta:

    :pwcpos_cmd: Commanded pupil wheel position angle (degrees).
    :pwcpos_bounds: Two-element array giving the valid range of the pupil wheel position angle (degrees).
    :order2_separation_cutoff: Wavelengths short of which to consider order 2 to be
        well-separated from order 1.
        Two values are defined. The first marks the shortest wavelength to include
        in the combined orders 1 and 2 extraction; the second marks the longest wavelength
        to include in the "well-separated" order 2 extraction. The values are different so
        that the extractions overlap in wavelength a bit, ameliorating numerical issues.
    :soss_xdim: Number of columns in the SOSS subarray.
    :soss_ydim: Number of rows in the SOSS subarray.
    :subarray_ymin: Starting y-pixel row of the SOSS subarray relative to the full detector frame.
    :wavemap_nwl: Number of wavelength grid points used to construct the 2D wavemaps.
    :wavemap_wlmin: Minimum wavelength (microns) of the wavemap wavelength grid.
    :wavemap_wlmax: Maximum wavelength (microns) of the wavemap wavelength grid.

:traces: List of trace tables, one per spectral order, each containing:

    :spectral_order: Spectral order number.
    :trace: 2D array of shape (N, 2) giving the (x, y) centroid positions of the trace.
    :cutoff: Maximum x pixel column to consider for this spectral order.
    :short_cutoff: Minimum wavelength (microns) cutoff for trace validity.
    :pivot_x: X pixel coordinate of the pivot point used for the rotation transformation.
    :pivot_y: Y pixel coordinate of the pivot point used for the rotation transformation.
    :oversampling: Oversampling factor.

:wavecal_models: List of wavelength calibration model tables, one per spectral order, defining 2-D polynomials
    that map (scaled) x pixel and (scaled) pupil wheel offset to wavelength. Each table contains:

    :spectral_order: Spectral order number.
    :polynomial_degree: Degree of the polynomial wavelength calibration model.
    :coefficients: 1D array of polynomial coefficients.
    :coefficient_labels: 1D arrays of strings describing the coefficient ordering convention.
    :scale_extents: 2D array ``[[x_min, offset_min], [x_max, offset_max]]`` used to normalize inputs to [0, 1] before evaluating the polynomial.

:throughputs: List of throughput tables, one per spectral order, each containing:

    :spectral_order: Spectral order number.
    :wavelength: 1D array of wavelength values (microns).
    :throughput: 1D array of fractional throughput values corresponding to each wavelength.
