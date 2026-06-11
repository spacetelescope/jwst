import os
from pathlib import Path

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.modeling.models import Identity, Multiply
from gwcs import coordinate_frames as cf
from gwcs import wcs
from stcal.alignment.util import compute_s_region_keyword, sregion_to_footprint

import jwst
from jwst.datamodels.utils.wfss_multispec import make_wfss_multiexposure_spec3
from jwst.extract_1d.tests.helpers import mock_nirspec_fs_one_slit_func
from jwst.pipeline import Spec3Pipeline
from jwst.stpipe import Step

N_SOURCES = 5
N_EXPOSURES = 4
INPUT_WFSS = "mock_wfss_x1d.fits"
INPUT_WFSS_2 = "mock_wfss_2_x1d.fits"
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


def wfss_one_spec(source_id):
    spec = dm.WFSSSpecModel()
    spec.source_id = source_id
    spec.exposure_time = 7.0
    spec.integration_time = 7.0
    spec.dispersion_direction = 3
    spec.spectral_order = 1
    spec.s_region = (
        "POLYGON ICRS  247.883569817 30.206692493 247.901987783 "
        "30.174116268 247.864126916 30.158804440 247.846405241 30.190721550"
    )
    # fmt: off
    spec_table = np.array(
        [(source_id, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, "POINT", 0, 0, 0, 0, 0, 0, 0, 0),
         (source_id, 5.5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, "POINT", 0, 0, 0, 0, 0, 0, 0, 0),
         (source_id, 10, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, "POINT", 0, 0, 0, 0, 0, 0, 0, 0),
         ],
        dtype=[
            ('SOURCE_ID', 'int32'),
            ('WAVELENGTH', 'float64'),
            ('FLUX', 'float64'),
            ('FLUX_ERROR', 'float64'),
            ('FLUX_VAR_POISSON', 'float64'),
            ('FLUX_VAR_RNOISE', 'float64'),
            ('FLUX_VAR_FLAT', 'float64'),
            ('SURF_BRIGHT', 'float64'),
            ('SB_ERROR', 'float64'),
            ('SB_VAR_POISSON', 'float64'),
            ('SB_VAR_RNOISE', 'float64'),
            ('SB_VAR_FLAT', 'float64'),
            ('DQ', 'uint32'),
            ('BACKGROUND', 'float64'),
            ('BKGD_ERROR', 'float64'),
            ('BKGD_VAR_POISSON', 'float64'),
            ('BKGD_VAR_RNOISE', 'float64'),
            ('BKGD_VAR_FLAT', 'float64'),
            ('NPIXELS', 'float64'),
            ('N_ALONGDISP', 'uint32'),
            ('SOURCE_TYPE', 'str'),
            ('SOURCE_XPOS', 'float64'),
            ('SOURCE_YPOS', 'float64'),
            ('SOURCE_RA', 'float64'),
            ('SOURCE_DEC', 'float64'),
            ('EXTRACT2D_XSTART', 'uint32'),
            ('EXTRACT2D_YSTART', 'uint32'),
            ('EXTRACT2D_XSTOP', 'uint32'),
            ('EXTRACT2D_YSTOP', 'uint32')])
    # fmt: on
    spec.spec_table = spec_table
    spec.spec_table.columns["WAVELENGTH"].unit = "um"
    return spec


def wfss_exp_spec3(source_id):
    spec0 = wfss_one_spec(source_id)
    multi = dm.WFSSMultiSpecModel()
    multi.meta.instrument.name = "NIRISS"
    multi.meta.instrument.detector = "ANY"
    multi.meta.exposure.type = "NIS_WFSS"
    multi.meta.observation.date = "2024-01-01"
    multi.meta.observation.time = "00:00:00.000"
    multi.meta.wcs = simple_wcs_wfss()
    for j in range(N_EXPOSURES):
        # create a new SpecModel for each exposure
        spec = spec0.copy()
        spec.filename = f"exposure_{j}.fits"
        spec.group_id = str(j + 1)
        multi.spec.append(spec)
    return multi


@pytest.fixture
def spec3_wfss_asn(tmp_cwd):
    model = make_wfss_multiexposure_spec3([wfss_exp_spec3(N_SOURCES - i) for i in range(N_SOURCES)])
    model.save(INPUT_WFSS)
    model2 = model.copy()
    new_s_region = _offset_sregion(model.spec[0].s_region, 0.005, 0.005)
    for spec in model2.spec:
        spec.group_id = "8"
        spec.s_region = new_s_region
    model2.save(INPUT_WFSS_2)
    os.system(f"asn_from_list -o {INPUT_ASN} --product-name test {INPUT_WFSS} {INPUT_WFSS_2}")


@pytest.fixture
def run_spec3_wfss(spec3_wfss_asn, monkeypatch):
    """
    Run the spec3 pipeline on a WFSS association.

    combine_1d step is mocked.
    pixel_replace is skipped.
    This only tests the pipeline logic.
    """

    def mock_combine1d(self, input_model, *args, **kwargs):
        """
        Mock the Combine1dStep process() method.

        Ensure it receives the right input type for a WFSS association.
        """
        if not isinstance(input_model, dm.MultiSpecModel):
            raise TypeError(f"Input to combine_1d is not a MultiSpecModel: {input_model}")
        output_model = dm.MultiCombinedSpecModel()
        spec = dm.SpecModel()
        output_model.spec.append(spec)
        output_model.meta.cal_step.combine_1d = "COMPLETE"
        return output_model

    monkeypatch.setattr("jwst.combine_1d.Combine1dStep.process", mock_combine1d)

    def mock_wfss_multiexposure_spec3(input_model):
        """
        Bypass reorganizing the output list.

        Ensure the input type is correct, which is equivalent to ensuring the result of
        the pipeline has the correct type.
        """
        if not isinstance(input_model, list):
            raise TypeError("Input to make_wfss_multiexposure_spec3 is not a list")
        if not all(isinstance(m, dm.MultiSpecModel) for m in input_model):
            raise TypeError(
                f"Input to make_wfss_multiexposure_spec3 is not a list of MultiSpecModel: {input_model}"
            )
        output_model = dm.WFSSMultiSpecModel()
        output_model.spec.append(dm.WFSSSpecModel())
        return output_model

    monkeypatch.setattr(
        jwst.pipeline.calwebb_spec3, "make_wfss_multiexposure_spec3", mock_wfss_multiexposure_spec3
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

    def mock_wfss_multispec_to_source(input_model):
        output_model = dm.MultiSpecModel()
        return [output_model]

    monkeypatch.setattr(
        jwst.pipeline.calwebb_spec3,
        "wfss_multispec_to_source",
        mock_wfss_multispec_to_source,
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
    with dm.open("test_x1d.fits") as x1d:
        assert (
            x1d.spec[0].s_region
            == "POLYGON ICRS  247.901987783 30.174116268 247.864126916 30.158804440 247.846405241 30.190721550 247.852683427 30.193419510 247.851405241 30.195721550 247.888569817 30.211692493 247.906987783 30.179116268 247.900617472 30.176539964"
        )


def test_spec3_nrs_fs(tmp_cwd):
    input_model = mock_nirspec_fs_one_slit_func()
    input_model.meta.filename = "test_spec3_x1d.fits"
    model_copy = input_model.copy()
    steps = {"outlier_detection": {"skip": True}, "resample_spec": {"skip": True}}
    Spec3Pipeline.call([input_model], steps=steps, save_results=True)

    # check for expected output
    assert Path("test_spec3_x1d.fits").exists()

    # make sure input model was not modified
    assert input_model.meta.cal_step.instance == model_copy.meta.cal_step.instance
    np.testing.assert_allclose(input_model.data, model_copy.data)
