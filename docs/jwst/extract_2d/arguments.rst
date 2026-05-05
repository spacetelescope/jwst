Step Arguments
==============

The ``extract_2d`` step has the following optional arguments:

``--slit_names`` (list of strings or integers)
  Names of specific slits to extract.
  If None, all known slits for the instrument mode to be extracted.
  Applies to all modes.

``--source_ids`` (list of strings or integers)
  Source IDs of specific slits to extract.
  If None, all known slits for the instrument to be extracted.
  For WFSS data, the selected source IDs correspond to the ``label`` column of the source catalog.
  ``slit_names`` and ``source_ids`` can be used at the same time, and duplicates will be filtered out.
  If either argument is specified, but no valid slits are identified, an error will be
  raised and the step will exit.  Applies to all modes.

``--tsgrism_extract_height`` (int)
  The cross-dispersion extraction size, in units of pixels. Only applies to TSO mode.

``--source_ra`` (list of floats)
  The RA coordinates (in decimal degrees) of specific sources to extract from
  the source catalog. Must match the length of ``source_dec``.
  ``source_ids`` can be used at the same time as ``source_ra`` and ``source_dec``;
  duplicates will be filtered out. Only applies to WFSS mode.

``--source_dec`` (list of floats)
  The Dec coordinates (in decimal degrees) of specific sources to extract from
  the source catalog. Must match the length of ``source_ra``.
  ``source_ids`` can be used at the same time as ``source_ra`` and ``source_dec``;
  duplicates will be filtered out. Only applies to WFSS mode.

``--source_max_sep`` (float)
  The maximum separation in arcseconds within which ``source_ra`` and ``source_dec``
  will be matched to sources in the catalog. If no source is found within this radius, a warning
  will be emitted and no source will be extracted corresponding to that ra, dec pair.
  Only applies to WFSS mode.

``--wfss_extract_half_height`` (int)
  The cross-dispersion half size of the extraction region, in pixels, applied to
  point sources. Only applies to WFSS mode.

``--wfss_mmag_extract`` (float)
  The minimum (faintest) magnitude object to extract, based on
  the value of ``isophotal_abmag`` in the source catalog.
  If None, sources of any magnitude will be extracted.
  Only applies to WFSS mode.

``--wfss_nbright`` (int)
  The number of brightest source catalog objects to extract.
  Can be used in conjunction with ``wfss_mmag_extract``. Only applies to WFSS mode.

``--extract_orders`` (list)
  The list of spectral orders to extract. If not explicitly specified, the default
  is taken from the ``wavelengthrange`` reference file. Applies to both WFSS and TSO modes.

``--grism_objects`` (list)
  A list of `~stdatamodels.jwst.transforms.GrismObject` that override
  the default extraction boxes. Only applies to WFSS mode.
