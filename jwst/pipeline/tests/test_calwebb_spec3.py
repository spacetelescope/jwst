import pytest
import os
from pathlib import Path

import stdatamodels.jwst.datamodels as dm
import jwst
from jwst.datamodels import SourceModelContainer
from jwst.stpipe import Step
from jwst.combine_1d._fileio import save_wfss_c1d
from jwst.extract_1d.tests.conftest import mock_niriss_wfss_l2, mock_nirspec_fs_one_slit, simple_wcs


INPUT_WFSS = "mock_wfss_cal.fits"
INPUT_WFSS_2 = "mock_wfss_2_cal.fits"
INPUT_ASN = "mock_wfss_asn.json"

@pytest.fixture
def spec3_wfss_asn(mock_niriss_wfss_l2, tmp_cwd):

    model = mock_niriss_wfss_l2
    for slit in model.slits:
        slit.meta.wcs = None # mock WCS coming in from fixture is not serializable
    model.save(INPUT_WFSS)
    model2 = model.copy()
    model2.meta.observation.exposure_number = "8"
    model2.save(INPUT_WFSS_2)
    os.system(f"asn_from_list -o {INPUT_ASN} --product-name test {INPUT_WFSS} {INPUT_WFSS_2}")


@pytest.fixture
def run_spec3_wfss(spec3_wfss_asn, monkeypatch):
    """
    Run the spec3 pipeline on a WFSS association.
    
    Extract_1d and combine_1d steps are mocked.
    Pixel_replace is skipped.
    This only tests the pipeline logic.
    """

    def mock_extract1d(self, input_model, *args, **kwargs):
        """
        Mock the Extract1dStep process() method.
        
        Ensure it receives the right input type for a WFSS association.
        """
        if not isinstance(input_model, SourceModelContainer):
            raise TypeError("Input to extract_1d is not a SourceModelContainer")
        output_model = dm.MultiSpecModel()
        spec = dm.SpecModel()
        spec.source_ra = 0.0
        spec.source_dec = 0.0
        output_model.spec.append(spec)
        output_model.meta.cal_step.extract_1d = "COMPLETE"
        return output_model
    monkeypatch.setattr("jwst.extract_1d.Extract1dStep.process", mock_extract1d)


    def mock_combine1d(self, input_model, *args, **kwargs):
        """
        Mock the Combine1dStep process() method.
        
        Ensure it receives the right input type for a WFSS association.
        """
        if not isinstance(input_model, dm.MultiSpecModel):
            raise TypeError("Input to extract_1d is not a SourceModelContainer")
        output_model = dm.MultiCombinedSpecModel()
        spec = dm.SpecModel()
        output_model.spec.append(spec)
        output_model.meta.cal_step.combine_1d = "COMPLETE"
        return output_model
    monkeypatch.setattr("jwst.combine_1d.Combine1dStep.process", mock_combine1d)


    def mock_save_wfss(input_model, filename):
        """
        Bypass saving the files.
        
        Ensure the input type is correct, which is equivalent to ensuring the result of
        the pipeline has the correct type.
        """
        if not isinstance(input_model, list):
            raise TypeError("Input to save_wfss is not a list")
        if not all(isinstance(model, (dm.MultiSpecModel, dm.MultiCombinedSpecModel)) for model in input_model):
            raise TypeError("Input to save_wfss is not a list of MultiSpecModel")
        Path(filename).write_text("Mocked save_wfss output")
    monkeypatch.setattr(jwst.pipeline.calwebb_spec3, "save_wfss_x1d", mock_save_wfss)
    monkeypatch.setattr(jwst.pipeline.calwebb_spec3, "save_wfss_c1d", mock_save_wfss)

    args = ["calwebb_spec3", INPUT_ASN,]
    Step.from_cmdline(args)


def test_spec3_wfss(run_spec3_wfss):
    """Smoke test to ensure pipeline runs for WFSS input."""
    files_created = os.listdir(".")
    assert "test_x1d.fits" in files_created
    assert "test_c1d.fits" in files_created