import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.wfss_contam.wfss_contam_step import WfssContamStep


@pytest.fixture
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
        slit.source_id = this_source["source_id"]
        slit.xstart = int(this_source["xcentroid"] - 10)
        slit.ystart = int(this_source["ycentroid"] - 10)
        slit.data = np.ones((20, 20))
        slit.xsize = slit.data.shape[1]
        slit.ysize = slit.data.shape[0]
        model.slits.append(slit)

    # manually change x,y offset because took transform from a real direct image, with different
    # pixel 0,0 than the mock data. This puts i=1, order 1 onto the real grism image
    model.slits[0].xstart = 2200
    model.slits[0].ystart = 1000

    fname = "multislit_model.fits"
    model.save(fname)
    return fname


def test_wfss_contam_step(multislitmodel, tmp_cwd_module):
    """
    Smoke test that the step runs.

    Right now none of the slits overlap with the simulated slits because of the incompatibility
    between a WCS taken from a random real image and the mock data.
    This could be fixed in the future by mocking the WCS object.
    """
    result = WfssContamStep.call(multislitmodel, save_simulated_image=True, save_contam_images=True)
    assert isinstance(result, dm.MultiSlitModel)
    assert result.meta.cal_step.wfss_contam == "COMPLETE"
