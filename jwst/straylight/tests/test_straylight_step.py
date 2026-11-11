"""
Unit tests for straylight step configuration
"""

import numpy as np
import pytest
from stdatamodels.jwst.datamodels import CubeModel, IFUImageModel
from gwcs import wcs
from jwst.straylight import StraylightStep
from jwst.tests.wcs_helpers import get_reference_files
from jwst.assign_wcs import miri

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
    image.meta.wcsinfo.v2_ref = -503.3406
    image.meta.wcsinfo.v3_ref = -318.864
    image.meta.wcsinfo.roll_ref = 234.28804423615944
    image.meta.wcsinfo.ra_ref = 206.88338859678234
    image.meta.wcsinfo.dec_ref = 49.3133194779185
    image.meta.wcsinfo.dispersion_direction = 2

    ref = get_reference_files(image)
    pipeline = miri.create_pipeline(image, ref)
    image.meta.wcs = wcs.WCS(pipeline)

    return image


def test_step_complete(miri_mrs_short):
    result = StraylightStep.call(miri_mrs_short)

    # Step is complete
    assert result.meta.cal_step.straylight == "COMPLETE"

    # Input is not modified
    assert result is not miri_mrs_short
    assert miri_mrs_short.meta.cal_step.straylight is None


def test_skip_missing_reffile(miri_mrs_short):
    result = StraylightStep.call(miri_mrs_short, override_mrsxartcorr="N/A")

    # Step is skipped
    assert result.meta.cal_step.straylight == "SKIPPED"

    # Input is not modified
    assert result is not miri_mrs_short
    assert miri_mrs_short.meta.cal_step.straylight is None


def test_call_straylight_mrsshort_tso(tmp_cwd, miri_mrs_short_tso):
    """Test step is skipped for MRS IFUSHORT TSO data"""
    result = StraylightStep.call(miri_mrs_short_tso)

    # Step is skipped
    assert result.meta.cal_step.straylight == "SKIPPED"

    # Input is not modified
    assert result is not miri_mrs_short_tso
    assert miri_mrs_short_tso.meta.cal_step.straylight is None


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
