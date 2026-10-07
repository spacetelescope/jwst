import numpy as np
import stdatamodels.jwst.datamodels as dm

from jwst.lib.stripe_utils import generate_substripe_ranges
from jwst.tests.nircam_rate_helpers import NRCA1_DHS_STRIPE_IDS

__all__ = ["make_mock_dhs_nrca1_regions"]


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


def make_mock_dhs_nrca1_regions(sci_model, tmp_path):
    """
    Write a mock subarray-size NRCA1 DHS regions reference file and return its path.

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
    # make a regions file matching the subarray size for the input data
    sub_ranges = generate_substripe_ranges(sci_model, science_frame=True)["subarray"]
    regions = np.zeros(sci_model.data.shape[-2:], dtype=np.float64)
    for i, stripe_id in enumerate(NRCA1_DHS_STRIPE_IDS):
        row_start, row_stop = sub_ranges[i]
        regions[row_start:row_stop, :] = stripe_id

    regions_path = tmp_path / "mock_nrca1_regions.asdf"
    model = dm.RegionsModel()
    model.regions = regions
    _populate_dhs_regions_metadata(model)
    model.save(str(regions_path))
    model.close()

    return str(regions_path)
