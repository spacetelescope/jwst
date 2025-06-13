"""
Unit tests for straylight step configuration
"""

from stdatamodels.jwst.datamodels import CubeModel, IFUImageModel

from jwst.straylight import StraylightStep
import numpy as np
import pytest


@pytest.fixture(scope="module")
def miri_mrs_short_tso():
    """Set up MIRI MRS SHORT TSO data"""

    image = CubeModel((5, 1024, 1032))
    image.data = np.random.random((5, 1024, 1032))
    image.meta.instrument.name = "MIRI"
    image.meta.instrument.detector = "MIRIFUSHORT"
    image.meta.exposure.type = "MIR_MRS"
    image.meta.instrument.channel = "12"
    image.meta.instrument.band = "SHORT"
    image.meta.filename = "test_miri_short_tso.fits"
    image.meta.observation.date = "2019-01-01"
    image.meta.observation.time = "10:10:10"
    return image


@pytest.fixture(scope="module")
def miri_mrs_short():
    """Set up MIRI MRS SHORT data"""

    image = IFUImageModel((1024, 1032))
    image.data = np.random.random((1024, 1032))
    image.meta.instrument.name = "MIRI"
    image.meta.instrument.detector = "MIRIFUSHORT"
    image.meta.exposure.type = "MIR_MRS"
    image.meta.instrument.channel = "12"
    image.meta.instrument.band = "SHORT"
    image.meta.filename = "test_miri_short.fits"
    image.meta.observation.date = "2019-01-01"
    image.meta.observation.time = "10:10:10"
    return image


def test_call_straylight_mrsshort_tso(tmp_cwd, miri_mrs_short_tso):
    """Test step is skipped for MRS IFUSHORT TSO data"""
    result = StraylightStep.call(miri_mrs_short_tso)
    assert result.meta.cal_step.straylight == "SKIPPED"


def test_step_save_shower_model(tmp_path, miri_mrs_short):
    model = miri_mrs_short

    StraylightStep.call(
        model,
        skip=False,
        clean_showers=True,
        save_results=True,
        save_shower_model=True,
        output_dir=str(tmp_path),
    )

    output_files = ["test_miri_short_straylightstep.fits", "test_miri_short_shower_model.fits"]

    for output_file in output_files:
        assert (tmp_path / output_file).exists()
