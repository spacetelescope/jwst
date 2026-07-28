from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.wfss_contam.wfss_contam_step import WfssContamStep


@pytest.fixture(scope="module")
def multislitmodel(
    tmp_path_factory, direct_image_with_gradient, segmentation_map, source_catalog, grism_wcs
):
    model = dm.MultiSlitModel()
    # add metadata
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.instrument.filter = "GR150C"
    model.meta.instrument.pupil = "F200W"
    model.meta.observation.date = "2023-05-30"
    model.meta.observation.time = "12:00:00"

    model.meta.exposure.type = "NIS_WFSS"
    model.meta.subarray.xsize = 2048
    model.meta.subarray.ysize = 2048

    # save direct image and segmentation map to file, then point model to those
    tmp_path = tmp_path_factory.mktemp("data")
    dim = str(tmp_path / "direct_image.fits")
    direct_image_with_gradient.save(dim)
    seg = str(tmp_path / "segmentation_map.fits")
    segmentation_map.save(seg)
    srccat = str(tmp_path / "source_catalog.ecsv")
    source_catalog.write(srccat, format="ascii.ecsv")
    model.meta.direct_image = dim
    model.meta.segmentation_map = seg
    model.meta.source_catalog = srccat

    # add a slit model with the grism WCS
    for i in range(len(source_catalog)):
        this_source = source_catalog[i]
        slit = dm.SlitModel()
        slit.meta.wcs = grism_wcs
        slit.meta.wcsinfo.spectral_order = 1
        slit.source_id = this_source["label"]
        slit.xstart = int(this_source["xcentroid"] - 10)
        slit.ystart = int(this_source["ycentroid"] - 10)
        slit.data = np.ones((20, 20))
        slit.xsize = slit.data.shape[1]
        slit.ysize = slit.data.shape[0]
        model.slits.append(slit)

    return model


def test_wfss_contam_step(tmp_cwd, multislitmodel):
    """
    Smoke test that the step runs with some user-defined options enabled.

    Right now none of the slits overlap with the simulated slits because of the incompatibility
    between a WCS taken from a random real image and the mock data.
    This could be fixed in the future by mocking the WCS object.
    """
    result = WfssContamStep.call(
        multislitmodel,
        output_file="multislit_model",
        save_simulated_image=True,
        save_contam_images=True,
        magnitude_limit=25,
        orders=[1],
    )
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "COMPLETE"
    assert Path("multislit_model_simul.fits").exists()
    assert Path("multislit_model_simul_slits.fits").exists()
    assert Path("multislit_model_contam.fits").exists()
    result.close()


def test_wfss_contam_step_defaults(tmp_cwd, multislitmodel):
    """
    Smoke test that the step runs with all default options.
    Also check that input is not modified by the step.
    """
    input_copy = multislitmodel.copy()

    result = WfssContamStep.call(multislitmodel)
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "COMPLETE"

    # Input is not modified
    assert result is not multislitmodel
    result.close()

    # Input data is not modified
    assert multislitmodel.meta.cal_step.wfss_contam is None
    i_modified = [
        i
        for i in range(len(multislitmodel.slits))
        if (not np.allclose(multislitmodel.slits[i].data, input_copy.slits[i].data))
    ]
    if len(i_modified) > 0:
        raise AssertionError(f"Slits modified: {i_modified}")


def test_wfss_contam_skip_maglimit(tmp_cwd, multislitmodel):
    """
    Test that the step is skipped if no sources meet the magnitude limit.
    """
    result = WfssContamStep.call(
        multislitmodel,
        save_simulated_image=True,
        save_contam_images=True,
        magnitude_limit=0,  # very bright, so no sources will meet this
        orders=[1],
    )
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "SKIPPED"
    result.close()


def test_wfss_contam_skip_bad_order(tmp_cwd, multislitmodel):
    """
    Test that the step is skipped if no valid spectral orders are found.
    """
    result = WfssContamStep.call(
        multislitmodel,
        save_simulated_image=True,
        save_contam_images=True,
        magnitude_limit=25,
        orders=[99],
    )
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "SKIPPED"
    result.close()


def test_wfss_contam_step_cube_direct_image(
    tmp_cwd, multislitmodel, direct_image_cube_with_gradient
):
    """
    Smoke test that the step completes when the direct image is a WFSSCubeModel.

    Reuses the multislitmodel fixture (slits, WCS, segmentation map, source catalog)
    but just swaps in the cube as the direct image.
    """
    direct_image_cube_with_gradient.save("direct_image_cube.fits")
    model = deepcopy(multislitmodel)
    model.meta.direct_image = str(Path("direct_image_cube.fits").resolve())
    result = WfssContamStep.call(model, magnitude_limit=25, orders=[1])
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "COMPLETE"
    result.close()


def test_wfss_contam_step_with_polyfit(tmp_cwd, multislitmodel):
    """Smoke test that the step completes when polyfit_degree and n_iterations are set."""
    result = WfssContamStep.call(multislitmodel, orders=[1], polyfit_degree=2, n_iterations=2)
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "COMPLETE"
    result.close()
