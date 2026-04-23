Description
===========

:Class: `~jwst.wfss_contam.wfss_contam_step.WfssContamStep`
:Alias: wfss_contam

The Wide Field Slitless Spectroscopy (WFSS) contamination correction
(``wfss_contam``) step is applied to grism exposures in an
attempt to correct effects due to overlapping spectral traces, which often
happens in observations of crowded fields. It is to be applied to individual
grism exposures in the latter stages of the :ref:`calwebb_spec2 <calwebb_spec2>`
pipeline.

Briefly, source fluxes from a direct image of the field are used
to simulate grism spectra for each source. Each source spectrum is then
corrected for contamination by subtracting the simulated spectra of nearby
sources. Details of the procedures and all input/output products are given
in the following sections.

Inputs
------

The method utilized to perform the correction requires several input data
products, including:

1. The grism data to be corrected. The step is applied near the end of the
   :ref:`calwebb_spec2 <calwebb_spec2>` pipeline, after the application of
   the :ref:`extract_2d <extract_2d_step>` and :ref:`srctype <srctype_step>`
   steps, but before the :ref:`photom <photom_step>` step. Thus individual
   2D cutouts exist for each identified source in the grism image, and the
   data are still in units of countrate.
2. The resampled direct image (:ref:`i2d <i2d>` product) of the field,
   usually obtained from the same WFSS observation as the grism image. The
   name of the direct image to use is retrieved from the "DIRIMAGE" keyword
   in the input grism image, which should've been populated at the
   beginning of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline from an
   entry in the "spec2" input ASN file.
3. The segmentation map (:ref:`segm <segm>` product) created from the direct image
   during :ref:`calwebb_image3 <calwebb_image3>` processing. The name of
   the segmentation map to use is retrieved from the "SEGMFILE" keyword in
   the input grism image, which should've been populated at the beginning
   of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline from an entry in
   the "spec2" input ASN file.

The Method
----------

Here we describe the steps used to perform the contamination correction:

1. First, a full-frame intermediate image, matching the size and shape of the
   grism image to be corrected, is created and populated with simulated spectra of
   all known sources in the field. The simulated spectra are created as follows:

   a. The segmentation (:ref:`segm <segm>`) file is searched for pixels with
      non-zero values and lists of pixels belonging to each source are created.
   b. The fluxes of each pixel in the lists are loaded from the direct image
      (:ref:`i2d <i2d>`), creating a list of per-pixel flux values for each source.
   c. A list of wavelength values is created for each source, which will be used to
      create the simulated spectra. The wavelength values span the range given by
      minimum and maximum wavelengths read from the WAVELENGTHRANGE reference file
      and are order-dependent.
   d. The direct image pixel locations and wavelengths for each source are transformed
      into dispersed pixel locations within the grism image using the WCS transforms
      of the input grism image.
   e. The flux of each direct image pixel belonging to each source is
      "dispersed" into the list of grism image pixel locations, thus creating a
      simulated spectrum.
   f. The initial simulated spectra are in flux-calibrated units, so each spectrum
      is divided by the sensitivity curve from the PHOTOM reference file, to convert
      the simulated spectra to units of countrates, thus matching the units of the
      observed grism data.
   g. The simulated spectrum for each source is stored in the full-frame image.
   h. Steps c-g are repeated for all spectral orders defined in the WAVELENGTHRANGE
      reference file.

2. 2D cutouts are created from the full-frame simulated grism image, matching the
   cutouts of each source in the input grism data.
3. For each source cutout, the simulated spectrum of the primary source is removed
   from the simulated cutout, leaving only the simulated spectra of any nearby
   contaminating sources.
4. The simulated contamination cutout is subtracted from the observed source cutout,
   thereby removing the signal from contaminating spectra.

Outputs
-------

There is one primary output and two optional outputs from the step:

1. The primary output is the contamination-corrected grism data, in the form of a
   `~stdatamodels.jwst.datamodels.MultiSlitModel` data model. In the :ref:`calwebb_spec2 <calwebb_spec2>`
   pipeline flow, this data model is passed along to the :ref:`photom <photom_step>` step
   for further processing.
2. If the step argument ``--save_simulated_image`` is set to `True`, the full-frame
   image containing all simulated spectra (the result of step 1 above) is saved to
   a file. See :ref:`wfss_contam_step_args`.
3. If the step argument ``--save_contam_images`` is set to `True`, the simulated
   contamination cutouts (the result of step 3 above) are saved to a file.
   See :ref:`wfss_contam_step_args`.

Polynomial Flux Modeling
------------------------

By default, each source is simulated with a spectrally flat flux model - that is, the
flux at every wavelength is taken directly from the direct image pixel values. 
When the step argument ``--polyfit_degree`` is set to an integer *N*, the step
fits a more flexible spectral model to each source. The procedure is:

1. In addition to the standard flat-spectrum simulation (the constant, degree-0 term),
   *N* additional grism-frame images are simulated for each source, one for each
   polynomial basis function :math:`\lambda^k` (:math:`k = 1, 2, \ldots, N`), where
   :math:`\lambda` is the wavelength. Recall that each dispersed-image pixel represents
   a linear combination of the contribution of several direct-image pixels at different
   wavelengths. These basis functions therefore must be computed before the dispersed image
   is discretized onto a pixel grid, i.e., just after the dispersion calculation.

2. For each source, the observed 2D spectrum is fit as a linear combination of
   these :math:`N+1` basis images, i.e.

   .. math::

      \text{observed} \approx c_0 \cdot B_0 + c_1 \cdot B_1 + \cdots + c_N \cdot B_N

   where :math:`B_0` is the flat-spectrum simulation and :math:`B_k` is the simulation
   driven by the :math:`\lambda^k` flux model.  The coefficients :math:`c_k` are
   determined by linear least-squares over all valid (finite, non-flagged) pixels.
3. The best-fit linear combination replaces the original simulation for that source, and
   this spectrally corrected simulation is used as the contamination model.

The ``--n_iterations`` argument controls how many times the polynomial fit is repeated.
On the first iteration the observed spectrum still contains contamination from neighboring
sources; on subsequent iterations it is replaced by the contamination-corrected spectrum
from the previous pass, so the polynomial flux fit is less biased by contamination
from other sources. Iteration has no effect when ``--polyfit_degree`` is not set.

Multiprocessing
---------------
The step can make use of multiple CPU cores to speed up the simulation of the
dispersed spectra. In short, the direct image pixels to be processed are divided into
chunks that are distributed to the available CPU cores. Each core processes its
assigned chunk of pixels, and the results are combined into the final full-frame
simulated grism image.
The number of cores to use can be set using the step argument ``--maximum_cores``,
and the maximum number of direct image pixels to be processed at once can be set using the
step argument ``--max_pixels_per_chunk``; see :ref:`wfss_contam_step_args`.
See :ref:`multiprocessing` for more details and examples of how to run a pipeline step
with multiprocessing enabled.