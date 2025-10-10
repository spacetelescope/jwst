import os

import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.modeling.models import Identity, Multiply
from gwcs import coordinate_frames as cf
from gwcs import wcs
from stcal.alignment.util import compute_s_region_keyword, sregion_to_footprint

import jwst
from jwst.datamodels import SourceModelContainer
from jwst.datamodels.utils.tests.wfss_helpers import wfss_multi
from jwst.extract_1d.tests.helpers import mock_nis_wfss_l2
from jwst.stpipe import Step

INPUT_WFSS = "mock_wfss_cal.fits"
INPUT_WFSS_2 = "mock_wfss_2_cal.fits"
INPUT_ASN = "mock_wfss_asn.json"


def simple_wcs_wfss():
    """
    Create a simple WCS whose transform is just to multiply by enough to make all
    points sit in different pixels.

    4 inputs and 4 outputs, to simulate WFSS (x, y, wavelength, order).
    """
    input_frame = cf.Frame2D(name="detector", axes_order=(0, 1))
    output_frame = cf.Frame2D(name="world", axes_order=(0, 1))
    transform = Multiply(0.005) & Multiply(0.005) & Identity(1) & Identity(1)
    return wcs.WCS(forward_transform=transform, input_frame=input_frame, output_frame=output_frame)


@pytest.fixture
def wfss_multiexposure():
    return wfss_multi()


@pytest.fixture
def mock_niriss_wfss_l2():
    model = mock_nis_wfss_l2()
    yield model
    model.close()


def _offset_sregion(s_region, dx, dy):
    """
    Offset an S_REGION string by dx, dy.

    Parameters
    ----------
    s_region : str
        The S_REGION string to offset.
    dx : float
        Offset in x direction.
    dy : float
        Offset in y direction.

    Returns
    -------
    str
        The offset S_REGION string.
    """
    footprint = sregion_to_footprint(s_region)
    footprint[:, 0] += dx
    footprint[:, 1] += dy
    return compute_s_region_keyword(footprint)


@pytest.fixture
def spec3_wfss_asn(mock_niriss_wfss_l2, tmp_cwd):
    model = mock_niriss_wfss_l2
    for slit in model.slits:
        slit.meta.wcs = simple_wcs_wfss()
    model.save(INPUT_WFSS)
    model2 = model.copy()
    model2.meta.group_id = "8"
    model2.meta.wcsinfo.s_region = _offset_sregion(model.meta.wcsinfo.s_region, 0.005, 0.005)
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
    assert (
        x1d.spec[0].s_region
        == "POLYGON ICRS  247.901987783 30.174116268 247.864126916 30.158804440 247.846405241 30.190721550 247.852683427 30.193419510 247.851405241 30.195721550 247.888569817 30.211692493 247.906987783 30.179116268 247.900617472 30.176539964"
    )
