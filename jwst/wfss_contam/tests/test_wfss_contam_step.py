import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.wfss_contam.wfss_contam_step import WfssContamStep


@pytest.fixture(scope="module")
def multislitmodel(
    tmp_cwd_module, direct_image_with_gradient, segmentation_map, source_catalog, grism_wcs
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

    # save direct image and segmentation map to file, then point model to those
    dim = "direct_image.fits"
    direct_image_with_gradient.save(dim)
    seg = "segmentation_map.fits"
    segmentation_map.save(seg)
    srccat = "source_catalog.ecsv"
    source_catalog.write(srccat, format="ascii.ecsv")  # , overwrite=True)
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

    fname = "multislit_model.fits"
    model.save(fname)
    return fname


def test_wfss_contam_step(multislitmodel, tmp_cwd_module):
    """
    Smoke test that the step runs with some user-defined options enabled.

    Right now none of the slits overlap with the simulated slits because of the incompatibility
    between a WCS taken from a random real image and the mock data.
    This could be fixed in the future by mocking the WCS object.
    """
    result = WfssContamStep.call(
        multislitmodel,
        save_simulated_image=True,
        save_contam_images=True,
        magnitude_limit=25,
        orders=[1],
    )
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "COMPLETE"
    assert (tmp_cwd_module / "multislit_model_simul.fits").exists()
    assert (tmp_cwd_module / "multislit_model_simul_slits.fits").exists()
    assert (tmp_cwd_module / "multislit_model_contam.fits").exists()


def test_wfss_contam_step_defaults(multislitmodel, tmp_cwd_module):
    """
    Smoke test that the step runs with all default options.
    """
    result = WfssContamStep.call(multislitmodel)
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "COMPLETE"


def test_wfss_contam_skip_maglimit(multislitmodel, tmp_cwd_module):
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


def test_wfss_contam_skip_bad_order(multislitmodel, tmp_cwd_module):
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


def test_output_is_not_input(multislitmodel, tmp_cwd_module):
    """Check that input is not modified by the step."""
    datamodel = dm.open(multislitmodel)
    input_copy = datamodel.copy()

    result = WfssContamStep.call(datamodel)
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "COMPLETE"

    # input is not modified
    assert result is not datamodel
    assert datamodel.meta.cal_step.wfss_contam is None
    any_modified = False
    for i in range(len(datamodel.slits)):
        # Input data is not modified
        np.testing.assert_allclose(datamodel.slits[i].data, input_copy.slits[i].data)

        # Output data may have been modified
        if not np.allclose(result.slits[i].data, datamodel.slits[i].data):
            any_modified = True

    # There was at least one slit updated in the output and not modified in the input
    assert any_modified
