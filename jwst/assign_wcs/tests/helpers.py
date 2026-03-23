import numpy as np
import stdatamodels.jwst.datamodels as dm

from jwst.lib.pipe_utils import generate_substripe_ranges

# Ordered from lowest to highest subarray row (= order in which stripes are
# packed into the packed subarray by generate_stripe_array).  The regions
# reference file is built with the same ordering, so the two stay in sync.
NRCA1_DHS_STRIPE_IDS = [10, 9, 8, 7]

__all__ = [
    "NRCA1_DHS_STRIPE_IDS",
    "make_mock_dhs_nrca1_rate",
    "make_mock_dhs_nrca1_regions",
    "make_mock_dhs_nrcalong_rate",
    "make_mock_dhs_nrcalong_regions",
]


def _populate_dhs_shared_metadata(model):
    """
    Populate metadata that is identical between NRCA1 and NRCALONG DHS modes.

    Updates *model* in place.
    """
    # Aperture
    model.meta.aperture.pps_name = "NRCA5_GRISM_F444W"

    # Exposure
    model.meta.exposure.type = "NRC_TSGRISM"
    model.meta.exposure.ngroups = 5

    # Observation
    model.meta.observation.date = "2026-05-03"
    model.meta.observation.time = "00:00:00.000"

    # Subarray
    model.meta.subarray.name = "SUB164STRIPE4_DHS"
    model.meta.subarray.fastaxis = -1
    model.meta.subarray.slowaxis = 2
    model.meta.subarray.num_superstripe = 0
    model.meta.subarray.repeat_stripe = 1
    model.meta.subarray.xstart = 1
    model.meta.subarray.xsize = 2048
    model.meta.subarray.ystart = 1
    model.meta.subarray.ysize = 164

    # WCS info
    model.meta.wcsinfo.ra_ref = 80.0
    model.meta.wcsinfo.dec_ref = -69.5
    model.meta.wcsinfo.v2_ref = 120.576793
    model.meta.wcsinfo.v3_ref = -527.501431
    model.meta.wcsinfo.roll_ref = 305.03225951982046
    model.meta.wcsinfo.velosys = -9378.83

    # Velocity aberration
    model.meta.velocity_aberration.scale_factor = 1.0


def _populate_dhs_regions_metadata(model):
    """
    Populate metadata that is identical between NRCA1 and NRCALONG regions files.

    Updates *model* in place.
    """
    model.meta.description = "Mock DHS regions for testing"
    model.meta.author = "test"
    model.meta.pedigree = "GROUND"
    model.meta.useafter = "2000-01-01T00:00:00"
    model.meta.instrument.name = "NIRCAM"
    model.meta.subarray.name = "SUB164STRIPE4_DHS"


def make_mock_dhs_nrca1_rate():
    """
    Return a mock DHS NRCA1 rate CubeModel.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.CubeModel`
        Mock NRCA1 DHS rate model.
    """
    model = dm.CubeModel((5, 164, 2048))
    _populate_dhs_shared_metadata(model)

    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.channel = "SHORT"
    model.meta.instrument.detector = "NRCA1"
    model.meta.instrument.filter = "F150W2"
    model.meta.instrument.pupil = "GDHS0"
    model.meta.instrument.module = "A"

    model.meta.subarray.interleave_reads1 = 1
    model.meta.subarray.multistripe_reads1 = 1
    model.meta.subarray.multistripe_reads2 = 40
    model.meta.subarray.multistripe_skips1 = 1526
    model.meta.subarray.multistripe_skips2 = 85

    model.meta.wcsinfo.siaf_xref_sci = 1024.5
    model.meta.wcsinfo.siaf_yref_sci = 24.5

    # Fill each packed stripe row range with the stripe's ID value so tests
    # can verify that the correct detector region ends up in each output slit.
    _, sub_ranges = generate_substripe_ranges(model, subarray_ranges=True)
    for i, stripe_id in enumerate(NRCA1_DHS_STRIPE_IDS):
        y0, y1 = sub_ranges[i]
        model.data[:, y0:y1, :] = stripe_id

    return model


def make_mock_dhs_nrca1_regions(sci_model, tmp_path):
    """
    Write a mock NRCA1 DHS regions reference file and return its path.

    Stripe IDs run from highest to lowest (10 - 7) as row index increases.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.CubeModel`
        The NRCA1 rate model whose multistripe parameters define the layout.
    tmp_path : `pathlib.Path`
        Writable temporary directory.

    Returns
    -------
    str
        Absolute path to the saved ASDF regions file.
    """
    reads1 = sci_model.meta.subarray.multistripe_reads1
    skips1 = sci_model.meta.subarray.multistripe_skips1
    reads2 = sci_model.meta.subarray.multistripe_reads2
    skips2 = sci_model.meta.subarray.multistripe_skips2

    regions = np.zeros((2048, 2048), dtype=np.float64)
    row_start = reads1 + skips1
    for stripe_id in NRCA1_DHS_STRIPE_IDS:
        regions[row_start : row_start + reads2, :] = stripe_id
        row_start += reads2 + skips2

    regions_path = tmp_path / "mock_nrca1_regions.asdf"
    model = dm.RegionsModel()
    model.regions = regions
    _populate_dhs_regions_metadata(model)
    model.save(str(regions_path))
    model.close()

    return str(regions_path)


def make_mock_dhs_nrcalong_rate():
    """
    Return a mock DHS NRCALONG rate CubeModel.

    Returns
    -------
    model : `~stdatamodels.jwst.datamodels.CubeModel`
        Mock NRCALONG DHS rate model.
    """
    model = dm.CubeModel((5, 164, 2048))
    _populate_dhs_shared_metadata(model)

    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.channel = "LONG"
    model.meta.instrument.detector = "NRCALONG"
    model.meta.instrument.filter = "F444W"
    model.meta.instrument.pupil = "GRISMR"
    model.meta.instrument.module = "A"

    model.meta.subarray.interleave_reads1 = 0
    model.meta.subarray.multistripe_reads1 = 1
    model.meta.subarray.multistripe_reads2 = 40
    model.meta.subarray.multistripe_skips1 = 971
    model.meta.subarray.multistripe_skips2 = 0

    model.meta.wcsinfo.siaf_xref_sci = 862
    model.meta.wcsinfo.siaf_yref_sci = 20.5

    # Fill each packed stripe row range with the same pattern
    # since for NRCALONG all readouts should be from the same physical detector region
    _, sub_ranges = generate_substripe_ranges(model, subarray_ranges=True)
    for sub_range in sub_ranges.values():
        y0, y1 = sub_range
        model.data[:, y0 + 5 : y1 - 5, :] = 1.0

    return model


def make_mock_dhs_nrcalong_regions(sci_model, tmp_path):
    """
    Write a mock NRCALONG DHS regions reference file and return its path.

    For NRCALONG the same detector band is read in every readout, so the
    regions map is a single nonzero stripe.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.CubeModel`
        The NRCALONG rate model whose multistripe parameters define the layout.
    tmp_path : `pathlib.Path`
        Writable temporary directory.

    Returns
    -------
    str
        Absolute path to the saved ASDF regions file.
    """
    reads1 = sci_model.meta.subarray.multistripe_reads1
    skips1 = sci_model.meta.subarray.multistripe_skips1
    reads2 = sci_model.meta.subarray.multistripe_reads2

    regions = np.zeros((2048, 2048), dtype=np.float64)
    row_start = reads1 + skips1
    regions[row_start : row_start + reads2, :] = 1

    regions_path = tmp_path / "mock_nrcalong_regions.asdf"
    model = dm.RegionsModel()
    model.regions = regions
    _populate_dhs_regions_metadata(model)
    model.save(str(regions_path))
    model.close()

    return str(regions_path)
