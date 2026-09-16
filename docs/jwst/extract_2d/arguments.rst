Step Arguments
==============
The ``extract_2d`` step has the following optional arguments.

The following arguments apply to **all modes** unless otherwise specified:

``--slit_names`` (list of strings or integers, default=None)
  Names of specific slits to extract.
  If None (default), all known slits for the instrument mode will be extracted.

``--source_ids`` (list of strings or integers, default=None)
  Source IDs of specific slits to extract.
  If None (default), all known slits for the instrument will be extracted.
  For WFSS data, the selected source IDs correspond to the ``label`` column of the source catalog.
  ``slit_names`` and ``source_ids`` can be used at the same time, and duplicates will be filtered out.
  If either argument is specified, but no valid slits are identified, an error will be
  raised and the step will exit.

The following arguments apply to **TSO and WFSS modes** only:

``--extract_orders`` (list, default=None)
  The list of spectral orders to extract. If not explicitly specified, the default
  is taken from the :ref:`bg_wlrange_reffile`.

The following arguments apply to **TSO mode** only:

``--tsgrism_extract_height`` (integer, default=None)
  The cross-dispersion extraction size, in units of pixels.

The following arguments apply to **WFSS modes** only:

``--source_ra`` (list of floats, default=None)
  The RA coordinates (in decimal degrees) of specific sources to extract from
  the source catalog. Must match the length of ``source_dec``.
  ``source_ids`` can be used at the same time as ``source_ra`` and ``source_dec``;
  duplicates will be filtered out.

``--source_dec`` (list of floats, default=None)
  The Dec coordinates (in decimal degrees) of specific sources to extract from
  the source catalog. Must match the length of ``source_ra``.
  ``source_ids`` can be used at the same time as ``source_ra`` and ``source_dec``;
  duplicates will be filtered out.

``--source_max_sep`` (float, default=2.0)
  The maximum separation in arcseconds within which ``source_ra`` and ``source_dec``
  will be matched to sources in the catalog. If no source is found within this radius, a warning
  will be emitted and no source will be extracted corresponding to that RA, Dec pair.

``--grism_objects`` (list, default=None)
  A list of `~stdatamodels.jwst.transforms.GrismObject` that override
  the default extraction boxes.

``--wfss_extract_half_height`` (int, default=5)
  The cross-dispersion half size of the extraction region, in pixels, applied to
  point sources.

``--wfss_mmag_extract`` (float, default=None)
  The minimum (faintest) magnitude object to extract, based on
  the value of ``isophotal_abmag`` in the source catalog.
  If None (default), sources of any magnitude will be extracted.

``--wfss_nbright`` (int, default=1000)
  The number of brightest source catalog objects to extract.
  Can be used in conjunction with ``wfss_mmag_extract``.
