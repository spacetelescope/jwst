import pytest
import os
import numpy as np

import stdatamodels.jwst.datamodels as dm
import jwst
from jwst.datamodels import SourceModelContainer
from jwst.datamodels.utils.tests.wfss_helpers import wfss_multi

from jwst.stpipe import Step
from jwst.extract_1d.tests.conftest import mock_nis_wfss_l2


INPUT_WFSS = "mock_wfss_cal.fits"
INPUT_WFSS_2 = "mock_wfss_2_cal.fits"
INPUT_ASN = "mock_wfss_asn.json"


@pytest.fixture
def wfss_multiexposure():
    return wfss_multi()


@pytest.fixture
def mock_niriss_wfss_l2():
    model = mock_nis_wfss_l2()
    yield model
    model.close()


@pytest.fixture
def spec3_wfss_asn(mock_niriss_wfss_l2, tmp_cwd):
    model = mock_niriss_wfss_l2
    for slit in model.slits:
        slit.meta.wcs = None  # mock WCS coming in from fixture is not serializable
    model.save(INPUT_WFSS)
    model2 = model.copy()
    model2.meta.group_id = "8"
    model2.save(INPUT_WFSS_2)
    os.system(f"asn_from_list -o {INPUT_ASN} --product-name test {INPUT_WFSS} {INPUT_WFSS_2}")


@pytest.fixture
def run_spec3_wfss(spec3_wfss_asn, monkeypatch, wfss_multiexposure):
    """
    Run the spec3 pipeline on a WFSS association.

    Extract_1d and combine_1d steps are mocked.
    Pixel_replace is skipped.
    This only tests the pipeline logic.
    """

    def mock_extract1d(self, input_model, *args, **kwargs):
        """
        Mock the Extract1dStep process() method.

        Ensure it receives the right input type for a WFSS association,
        and outputs the correct type.
        """
        if not isinstance(input_model, SourceModelContainer):
            raise TypeError("Input to extract_1d is not a SourceModelContainer")

        return wfss_multiexposure

    monkeypatch.setattr("jwst.extract_1d.Extract1dStep.process", mock_extract1d)

    def mock_combine1d(self, input_model, *args, **kwargs):
        """
        Mock the Combine1dStep process() method.

        Ensure it receives the right input type for a WFSS association.
        """
        if not isinstance(input_model, dm.WFSSMultiSpecModel):
            raise TypeError("Input to combine_1d is not a WFSSMultiSpecModel")
        output_model = dm.MultiCombinedSpecModel()
        spec = dm.SpecModel()
        output_model.spec.append(spec)
        output_model.meta.cal_step.combine_1d = "COMPLETE"
        return output_model

    monkeypatch.setattr("jwst.combine_1d.Combine1dStep.process", mock_combine1d)

    def mock_wfss_multiexposure(input_model):
        """
        Bypass reorganizing the output list.

        Ensure the input type is correct, which is equivalent to ensuring the result of
        the pipeline has the correct type.
        """
        if not isinstance(input_model, list):
            raise TypeError("Input to make_wfss_multiexposure is not a list")
        if not all(isinstance(m, dm.WFSSMultiSpecModel) for m in input_model):
            raise TypeError("Input to make_wfss_multiexposure is not a list of WFSSMultiSpecModel")
        output_model = dm.WFSSMultiSpecModel()
        output_model.spec.append(dm.WFSSSpecModel())
        return output_model

    monkeypatch.setattr(
        jwst.pipeline.calwebb_spec3, "make_wfss_multiexposure", mock_wfss_multiexposure
    )

    def mock_wfss_multicombined(input_model):
        """
        Bypass reorganizing the output list.

        Ensure the input type is correct, which is equivalent to ensuring the result of
        the pipeline has the correct type.
        """
        if not isinstance(input_model, list):
            raise TypeError("Input to make_wfss_multicombined is not a list")
        if not all(isinstance(m, dm.MultiCombinedSpecModel) for m in input_model):
            raise TypeError("Input to make_wfss_multicombined is not a list of MultiSpecModel")
        output_model = dm.WFSSMultiCombinedSpecModel()
        return output_model

    monkeypatch.setattr(
        jwst.pipeline.calwebb_spec3, "make_wfss_multicombined", mock_wfss_multicombined
    )

    args = [
        "calwebb_spec3",
        INPUT_ASN,
    ]
    Step.from_cmdline(args)


def test_spec3_wfss(run_spec3_wfss):
    """Smoke test to ensure pipeline runs for WFSS input."""
    files_created = os.listdir(".")
    assert "test_x1d.fits" in files_created
    assert "test_c1d.fits" in files_created
    x1d = dm.open("test_x1d.fits")
    assert len(x1d.spec[0].s_region) > 0
