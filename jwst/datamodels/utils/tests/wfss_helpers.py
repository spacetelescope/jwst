import numpy as np
import stdatamodels.jwst.datamodels as dm
from jwst.datamodels.utils.wfss_multispec import make_wfss_multiexposure


N_SOURCES = 5
N_EXPOSURES = 4
N_ROWS = 3


def example_spec():
    """
    Create a mock SpecModel with spec_table generated from the spec schema.

    Returns
    -------
    spec : dm.SpecModel
        A SpecModel with a mock WCS and a spec_table with N_ROWS rows.
        The spec_table has columns "WAVELENGTH" and "FLUX".
    """

    def mock_wcs(*args, **kwargs):  # Noqa: ARG001
        return 0.0, 0.0, 0.0

    spec = dm.SpecModel()
    spectable_dtype = spec.schema["properties"]["spec_table"]["datatype"]
    recarray_dtype = [(d["name"], d["datatype"]) for d in spectable_dtype]
    spec.meta.wcs = mock_wcs
    spec_table = np.recarray((N_ROWS,), dtype=recarray_dtype)
    spec_table["WAVELENGTH"] = np.linspace(1.0, 10.0, N_ROWS)
    spec_table["FLUX"] = np.ones(N_ROWS)
    spec.spec_table = spec_table
    spec.spec_table.columns["wavelength"].unit = "um"
    return spec


def _add_multispec_meta(spec):
    """
    Add specmeta attributes to a spec-like ObjectNode inside a MultiSpecModel.

    This only includes attributes that do NOT need to change for each source/exposure
    for the purposes of this test.

    input spec is updated in place.
    """
    spec.source_type = "POINT"
    spec.source_ra = 0.0
    spec.source_dec = 0.0
    spec.extract2d_xstart = 0.0
    spec.extract2d_ystart = 0.0
    spec.extract2d_xstop = 0.0
    spec.extract2d_ystop = 0.0
    spec.extraction_xstart = 0.0
    spec.extraction_ystart = 0.0
    spec.extraction_xstop = 0.0
    spec.extraction_ystop = 0.0


def wfss_spec2_multi():
    """
    Set up a MultiSpecModel object that looks like outputs from extract_1d during calwebb_spec2.

    Each of the spectra is from the SAME exposure, but DIFFERENT sources.

    Returns
    -------
    multi : dm.MultiSpecModel
        A MultiSpecModel with N_SOURCES in its spec list.
    """
    spec0 = example_spec()
    multi = dm.MultiSpecModel()
    for i in range(N_SOURCES):
        # create a new SpecModel for each source
        spec = spec0.copy()
        spec.meta.filename = "exposure_1.fits"  # all sources in the same exposure
        spec.meta.group_id = "1"  # all sources in the same exposure
        spec.source_id = N_SOURCES - i  # reverse the order to test sorting
        spec.name = str(spec.source_id)
        _add_multispec_meta(spec)
        multi.spec.append(spec)

    return multi


def wfss_spec3_multi():
    """
    Set up a MultiSpecModel object that looks like outputs from extract_1d during calwebb_spec3.

    Each of the spectra is from the SAME source, but DIFFERENT exposures.

    Returns
    -------
    multi : dm.MultiSpecModel
        A MultiSpecModel with N_EXPOSURES in its spec list.
    """
    spec0 = example_spec()
    multi = dm.MultiSpecModel()
    multi.meta.exposure.exposure_time = 7.0
    multi.meta.exposure.integration_time = 7.0
    for j in range(N_EXPOSURES):
        # create a new SpecModel for each exposure
        spec = spec0.copy()
        spec.meta.filename = f"exposure_{j}.fits"
        spec.meta.group_id = str(j + 1)
        spec.source_id = 999  # all sources are the same
        spec.dispersion_direction = 3
        spec.name = str(spec.source_id)
        _add_multispec_meta(spec)
        multi.spec.append(spec)

    return multi


def wfss_multi():
    """
    Make a MultiExposureSpecModel object with N_EXPOSURES exposures and N_SOURCES sources.

    Returns
    -------
    dm.MultiExposureSpecModel
        A MultiExposureSpecModel with N_EXPOSURES in its spec list,
        each containing N_SOURCES sources.
    """
    inputs_list = []
    source0 = wfss_spec3_multi()
    for i in range(N_SOURCES):
        this_source = source0.copy()
        for spec in this_source.spec:
            spec.source_id = N_SOURCES - i
        inputs_list.append(this_source)
    output_model = make_wfss_multiexposure(inputs_list)
    return output_model


def wfss_comb():
    """
    Make a MultiCombinedSpecModel object with N_SOURCES sources, N_ROWS rows, and 2 orders.

    This looks like the output of combine_1d.

    Returns
    -------
    dm.MultiCombinedSpecModel
        A MultiCombinedSpecModel with 2 SpecModels from different spectral orders,
        each containing a spec_table with N_ROWS rows.
    """
    multi = dm.MultiCombinedSpecModel()
    spec = dm.CombinedSpecModel()
    _add_multispec_meta(spec)
    spectable_dtype = spec.schema["properties"]["spec_table"]["datatype"]
    recarray_dtype = [(d["name"], d["datatype"]) for d in spectable_dtype]
    spec_table = np.recarray((N_ROWS,), dtype=recarray_dtype)
    spec_table["WAVELENGTH"] = np.linspace(1.0, 10.0, N_ROWS)
    spec_table["FLUX"] = np.ones(N_ROWS)
    spec.spec_table = spec_table
    spec.spec_table.columns["wavelength"].unit = "um"
    spec.dispersion_direction = 3

    spec.spectral_order = 1
    spec2 = spec.copy()
    spec2.spectral_order = 2
    multi.spec.append(spec)
    multi.spec.append(spec2)
    return multi
