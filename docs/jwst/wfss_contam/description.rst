Description
============

:Class: `jwst.wfss_contam.WfssContamStep`
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

 1) The grism data to be corrected. The step is applied near the end of the
    :ref:`calwebb_spec2 <calwebb_spec2>` pipeline, after the application of
    the :ref:`extract_2d <extract_2d_step>` and :ref:`srctype <srctype_step>`
    steps, but before the :ref:`photom <photom_step>` step. Thus individual
    2D cutouts exist for each identified source in the grism image, and the
    data are still in units of countrate.

 2) The resampled direct image (:ref:`i2d <i2d>` product) of the field,
    usually obtained from the same WFSS observation as the grism image. The
    name of the direct image to use is retrieved from the "DIRIMAGE" keyword
    in the input grism image, which should've been populated at the
    beginning of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline from an
    entry in the "spec2" input ASN file.

 3) The segmentation map (:ref:`segm <segm>` product) created from the direct image
    during :ref:`calwebb_image3 <calwebb_image3>` processing. The name of
    the segmentation map to use is retrieved from the "SEGMFILE" keyword in
    the input grism image, which should've been populated at the beginning
    of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline from an entry in
    the "spec2" input ASN file.

The Method
----------
Here we describe the steps used to perform the contamination correction.

 1) First, a full-frame intermediate image, matching the size and shape of the
    grism image to be corrected, is created and populated with simulated spectra of
    all known sources in the field. The simulated spectra are created as follows:

    a) The segmentation (:ref:`segm <segm>`) file is searched for pixels with
       non-zero values and lists of pixels belonging to each source are created.
    b) The fluxes of each pixel in the lists are loaded from the direct image
       (:ref:`i2d <i2d>`), creating a list of per-pixel flux values for each source.
    c) A list of wavelength values is created for each source, which will be used to
       create the simulated spectra. The wavelength values span the range given by
       minimum and maximum wavelengths read from the WAVELENGTHRANGE reference file
       and are order-dependent.
    d) The direct image pixel locations and wavelengths for each source are transformed
       into dispersed pixel locations within the grism image using the WCS transforms
       of the input grism image.
    e) The flux of each direct image pixel belonging to each source is
       "dispersed" into the list of grism image pixel locations, thus creating a
       simulated spectrum.
    f) The initial simulated spectra are in flux-calibrated units, so each spectrum
       is divided by the sensitivity curve from the PHOTOM reference file, to convert
       the simulated spectra to units of countrates, thus matching the units of the
       observed grism data.
    g) The simulated spectrum for each source is stored in the full-frame image.
    h) Steps c-g are repeated for all spectral orders defined in the WAVELENGTHRANGE
       reference file.
 2) 2D cutouts are created from the full-frame simulated grism image, matching the
    cutouts of each source in the input grism data.
 3) For each source cutout, the simulated spectrum of the primary source is removed
    from the simulated cutout, leaving only the simulated spectra of any nearby
    contaminating sources.
 4) The simulated contamination cutout is subtracted from the observed source cutout,
    thereby removing the signal from contaminating spectra.

Outputs
-------
There is one primary output and two optional outputs from the step:

 1) The primary output is the contamination-corrected grism data, in the form of a
    `~jwst.datamodels.MultiSlitModel` data model. In the :ref:`calwebb_spec2 <calwebb_spec2>`
    pipeline flow, this data model is passed along to the :ref:`photom <photom_step>` step
    for further processing.

 2) If the step argument `--save_simulated_image` is set to `True`, the full-frame
    image containing all simulated spectra (the result of step 1 above) is saved to
    a file. See :ref:`wfss_contam_step_args`.

 3) If the step argument `--save_contam_images` is set to `True`, the simulated
    contamination cutouts (the result of step 3 above) are saved to a file.
    See :ref:`wfss_contam_step_args`.
