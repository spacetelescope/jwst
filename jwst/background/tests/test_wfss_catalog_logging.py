import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.background import background_sub_wfss


def test_catalog_mask_does_not_log_missing_catalog(monkeypatch, caplog):
    model = datamodels.ImageModel(
        data=np.ones((4, 4)),
        err=np.ones((4, 4)),
        dq=np.zeros((4, 4), dtype=np.uint32),
    )
    model.meta.source_catalog = "catalog.ecsv"
    model.meta.wcsinfo.dispersion_direction = 1
    bkg = datamodels.ImageModel(
        data=np.ones((4, 4)),
        dq=np.zeros((4, 4), dtype=np.uint32),
    )
    monkeypatch.setattr(background_sub_wfss.datamodels, "open", lambda _: bkg)
    monkeypatch.setattr(
        background_sub_wfss, "get_subarray_model", lambda science, reference: reference
    )
    monkeypatch.setattr(
        background_sub_wfss,
        "_mask_from_source_cat",
        lambda *args, **kwargs: np.zeros((4, 4), dtype=bool),
    )
    with caplog.at_level(logging.WARNING, logger=background_sub_wfss.__name__):
        result = background_sub_wfss.subtract_wfss_bkg(
            model, "background.fits", "wavelengthrange.asdf"
        )
    try:
        assert result.meta.cal_step.bkg_subtract == "FAILED"
        assert "No source_catalog found" not in caplog.text
        assert "No sources will be masked" not in caplog.text
        assert "Not enough background pixels" in caplog.text
    finally:
        result.close()
