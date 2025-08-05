import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.extract_1d.tests.helpers import mock_nirspec_fs_one_slit_func


@pytest.fixture
def mock_nirspec_fs_one_slit():
    """
    Make mock_nirspec_fs_one_slit_function available as a fixture.

    Yields
    ------
    SlitModel
        The mock model.
    """
    model = mock_nirspec_fs_one_slit_func()
    # mock WCS was originally a function which is not serializable
    model.meta.wcs = None
    yield model
    model.close()


@pytest.fixture
def mock_nirspec_mos(mock_nirspec_fs_one_slit):
    """
    Mock three slits in NIRSpec MOS mode.

    Yields
    ------
    MultiSlitModel
        The mock model.
    """
    model = dm.MultiSlitModel()
    model.meta.instrument.name = "NIRSPEC"
    model.meta.instrument.detector = "NRS1"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.exposure.type = "NRS_MSASPEC"
    model.meta.exposure.nints = 1

    nslit = 3
    for i in range(nslit):
        slit = mock_nirspec_fs_one_slit.copy()
        slit.name = str(i + 1)
        slit.source_id = i + 1
        model.slits.append(slit)

    yield model
    model.close()


@pytest.fixture
def mock_input(mock_nirspec_mos):
    """
    Give the slits different source IDs.

    Result is 5 sources, where the middle one is present in all three exposures.

    Returns
    -------
    list[MultiSlitModel]
        List of mock MultiSlitModel objects with different source IDs.
    """
    inputs = []
    for i in range(3):
        this_model = mock_nirspec_mos.copy()
        this_model.meta.filename = f"mock_nirspec_mos_{i}.fits"
        for slit in this_model.slits:
            slit.source_id += i
        inputs.append(this_model)
    return inputs
